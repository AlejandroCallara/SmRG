function [p_tot,a1,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_multimix(Vin,K0_double,vB_double,mu_sk,rk)
%           fits model described in [1] on data input Vin. The model
%           describes a single class k of signal pixels (negative-binomial)
%           and the background (gaussian).For more details refer to [1].
%           For negative binomial parameter update refer to [2].
%
% Syntax:
%           [p_tot,a1]=SmRG_mixtureModelFitting_multimix(vin)
%           [p_tot,a1,K0_double,vB,mu_sk,rk]=SmRG_mixtureModelFitting_multimix(vin)
% Inputs:
%           Vin: 2D image (type: double)*
%
% Outputs:
%           p_tot:  posterior probabilities for each element in Vin(:) of
%                   belonging to neg-bin distribution. The posterior
%                   probabilities of belonging to the gaussian distribution
%                   are simply:
%                               p_gaussian = 1-p_tot
%           a1:     posterior probabilities of voxels belonging to the
%                   negbin distribution
%           K0_double: Gaussian mean
%           vB:     Gaussian variance
%           mu_sk:  negbin mean
%           rk:     negbin dispersion parameter
%
%
%           * note: the function allows any N-dimensional input. When
%           called from SmRG_regionGrowing Vin is a 2D image. However, one
%           can use this function to fit the mixture model of [1] on any
%           kind of input, regardless of dimensionality.
%
%[1]: Calapez,A. and Rosa,A. (2010) A statistical pixel intensity model
%     for segmentation of confocal laser scanning microscopy images.
%     IEEE Trans. Image Process., 19, 2408–2418.
%[2]: Chunmao Huang et al 2019 J. Phys.: Conf. Ser. 1324 012093

% % check inputs
% if nargin<1
%     help SmRG_mixtureModelFitting
%     return
% end
%
% % if nargin<2
% %     % if background not specified never do check on background
% %     background = 0;
% % end
%
% if sum(Vin(:)<0)
%     error ('only images with positive values are allowed')
% end

% if ~isrow(Vin)
%     Vin = Vin';
% end
% if ~isrow(Vin)
%     mu_sk = mu_sk';
% end,mu_sk,rk

nnegbin = 4;
coder.varsize('log_lik_a');
coder.varsize('pa');
coder.varsize('a');
coder.varsize('a1');
coder.varsize('mu_sk');
coder.varsize('v_sk');
coder.varsize('rk');
coder.varsize('lambdak');


% check for nans in data
l = length(Vin);
% pre-allocate p_tot
% p_tot=zeros(1,l);
a  = zeros(l,nnegbin);
a1 = zeros(l,nnegbin);
b  = zeros(l,1);
b1 = zeros(l,1);

i_ok = 1:length(Vin);
V = uint16(Vin);
% if no initial conditions then do it data driven
if nargin <2
    % Gaussian initial condition
    K0_double =2*(min(Vin))+5; % mean
    vB_double = 2*K0_double;    % variance
    
    % Neg-Bin initial condition
    mu_sk = mean(Vin); % mean
    rk = 10*ones(1,nnegbin);
    
    fact = 1;
    for inegbin = 2:nnegbin
        fact = fact+0.1;
        mu_sk(inegbin)=fact*mu_sk(1);
    end
    v_sk = (mu_sk.^2+mu_sk.*rk)./rk;
    pk   = mu_sk./v_sk;
else %else get them
    media_vin = mean(Vin);
    mu_sk_tmp = ones(1,nnegbin)*media_vin; % mean
    mu_sk_tmp(1:length(mu_sk))=mu_sk;
    rk_tmp = 10*ones(1,nnegbin);
    rk_tmp(1:length(mu_sk)) = rk;
    
    
    mu_sk = mu_sk_tmp;
    rk = rk_tmp;
    v_sk = (mu_sk.^2+mu_sk.*rk)./rk;
    pk   = mu_sk./v_sk;
end

% Mixing prior probability (0.5 0.5/nnegbin)
pb   = 0.5;
pa   = (1-pb)/nnegbin*ones(1,nnegbin);

% EM initial conditions
n_of_pixel    = length(Vin);
log_lik_a = zeros(5000,nnegbin);
log_lik_b = zeros(5000,1);
log_lik_b(1)  = 0;
log_lik_b(2)  = 10;
for inegbin = 1:nnegbin
    log_lik_a(1,inegbin)  = 0;
    log_lik_a(2,inegbin)  = 10;
end

tol   = 1e-3;
count = 2;

% reset tol_nbin and tol_gauss
tol_nbin = abs(log_lik_a(count,:)-log_lik_a(count-1,:))./abs(log_lik_a(count,:));
tol_nbin = max(tol_nbin);

% update for tollerance check gaussian
tol_gauss = abs(log_lik_b(count)-log_lik_b(count-1))/abs(log_lik_b(count));

%% mixture, EM
% iterates until convergence
reduced = true;
p_tot =a1;

while reduced
    % initialize log likelihood
    psi_a = zeros(n_of_pixel,nnegbin);
    psi_b = zeros(n_of_pixel,1);
    while (tol_nbin>tol && tol_gauss>tol)
        % E-step
        for ii=1:n_of_pixel
            
            Vtmp = (Vin(ii));
            
            %  Gaussian log likelihood
            psi_b_tmp = (-SmRG_normlike([K0_double sqrt(vB_double)],Vtmp));
            psi_b(ii,1) = psi_b_tmp;
                        
            %  Neg-Bin log likelihood
            for inegbin = 1:nnegbin
                psi_a_tmp = (-nbinlike_mu([mu_sk(inegbin) 1/rk(inegbin)],(Vtmp)));
                psi_a(ii,inegbin) = psi_a_tmp;
            end
            
            % find who's bigger between [log(P(x|a))+log(P(a)] e [log(P(x|b))+log(P(b)]
            % to use the logarithm property of sum: log(sum(ai))=log_b(a0)+log_b(1+sum(b^(log_b(ai)-log_b(a0))
            % with a0>a1>...>aN
            % in this way it is possible to evaluate the log posterior without
            % calculating the exponential: this helps for numerical stability
            % (sometimes exponential blows up)
            
            % sorting
            psort = zeros(1,0);
            for inegbin = 1:nnegbin
                psort = [psort psi_a(ii,inegbin)+log(pa(inegbin))];
            end
            psort = [psort psi_b_tmp+log(pb)];
            [psort isort]=sort(psort,'descend');
            
            % exponential
            sumexp = 0;
            for inegbin = 2:nnegbin
                sumexp = sumexp+exp(psort(inegbin)-psort(1));
            end
            
            % den
            den = psort(1)+log(1+sumexp);
            
            % log posterior gaussian
            b(ii,1) = ((psi_b_tmp)+log(pb))-den;
            
            % log posterior nbin
            for inegbin = 1:nnegbin
                a(ii,inegbin) = psi_a(ii,inegbin)+log(pa(inegbin))-den;
            end
        end
        b1 = exp(b);
        a1 = exp(a);
        
        % M-step
        
        % update Gaussian Parameters
        K0_double = sum(b1'.*Vin)/sum(b1);
        vB_double = sum(b1'.*((Vin-K0_double).^2))/sum(b1);
        
        % update Neg-Bin Parameters
        dk = zeros(n_of_pixel,nnegbin);
        for inegbin = 1:nnegbin
            dk(:,inegbin) = rk(inegbin).*(psi(rk(inegbin)+Vin)-psi(rk(inegbin)));
        end
        tauk    = a1;
        lambdak = sum(tauk.*dk)./sum(tauk);
        betak   = 1- 1./(1-pk)- 1./(log(pk));
        % pk =  (betak.*sum(tauk.*dk))./sum(tauk.*(Vin'-(1-betak).*dk));
        for inegbin = 1:nnegbin
            pk(inegbin) = (betak(inegbin).*sum(tauk(:,inegbin).*dk(:,inegbin)))./sum(tauk(:,inegbin).*(Vin'-(1-betak(inegbin)).*dk(:,inegbin)));
        end
        rk      = -lambdak./log(pk);
        v_sk    = rk.*(1-pk)./(pk.^2);
        mu_sk   = pk.*v_sk;
        
        % today's posterior becomes tomorrow's prior :)
        meana1 = mean(a1);
        pa = meana1;
        meanb1 =mean(b1);
        pb = meanb1;
        
        % USEFUL IN DEBUG MODE: plots model fitting on data hist
        % figure(1),
        % histogram(Vin);hold on
        % plot(50000*normpdf(1:max(Vin(:)),K0_double,sqrt(vB_double)))
        % for inegbin = 1:nnegbin
        %     plot(50000*nbinpdf_mu(1:max(Vin(:)),mu_sk(inegbin),1/rk(inegbin)))
        % end
        % drawnow
        % hold off
        
        % update for tollerance check negbin
        count= count+1;
        tol_nbin = [];
        for inegbin = 1:nnegbin
            log_lik_a(count,inegbin)=mean(psi_a(:,inegbin));
        end
        tol_nbin = abs(log_lik_a(count,:)-log_lik_a(count-1,:));
        tol_nbin = max(tol_nbin(:));
        
        % update for tollerance check gaussian
        log_lik_b(count)=mean(psi_b);
        tol_gauss = (abs(log_lik_b(count)-log_lik_b(count-1)))/abs(log_lik_b(count));
        % prevent infinite loop
        if count > 4999
            break
        end
    end
    p_tot =a1;
    fprintf('convergence reached after %d steps \n', int16(count))
    % end of EM algorithm
    
    % now, check if some clusters are the same
    mus    = mu_sk;
    sigmas = v_sk;
    
    % elbow method kmeans params
    MAX_ELBOW = nnegbin;
    CUTOFF    = 0.95;
    REPEATS   = 3;
    [v_opt,C,SUMD,K]=SmRG_kmeans_opt([mus',sigmas'],MAX_ELBOW,CUTOFF,REPEATS);
    
    % get cluster idx to merge
    idx= struct('to_merge',repmat({false(nnegbin,1)},K,1));
    
    if K<nnegbin
        ck_count = 0;
        while ck_count<K
            ck_count = ck_count+1;
            idx(ck_count).to_merge =v_opt==ck_count;
        end
        
        mus_tmp    = zeros(1,K);
        sigmas_tmp = zeros(1,K);
        p_tmp     = zeros(1,K);
        for ick = 1:length(idx)
            idx_to_merge = idx(ick).to_merge;
            
            % merge
            mus_merge    = mus(idx_to_merge);
            sigmas_merge = sigmas(idx_to_merge);
            p_merge      = pa(idx_to_merge);
            
            mus_merge_mean    = mean(mus_merge,2);
            sigmas_merge_mean = mean(sigmas_merge,2);
            p_merge_sum       = sum(p_merge(:));
            
            mus_tmp(ick)    = mus_merge_mean;
            sigmas_tmp(ick) = sigmas_merge_mean;
            p_tmp(ick)      = p_merge_sum;
        end
        
        mu_sk = mus_tmp;
        v_sk  = sigmas_tmp;
        pa    = p_tmp;
        
        % update params
        rk   = mu_sk.^2./(v_sk-mu_sk);
        pk   = mu_sk./v_sk;
    end
    
    % if merge clusters then run again EM
    if size(pa,2)==nnegbin
        reduced = false
    else
        disp('some clusters were the same...merging...')

        nnegbin = K;
        reduced = true;
        
        % set again initial conditions

        % log_likelihood
        log_lik_a = zeros(5000,nnegbin);
        for inegbin = 1:nnegbin
            log_lik_a(1,inegbin)  = 0;
            log_lik_a(2,inegbin)  = 10;
        end
        log_lik_b = zeros(5000,1);
        log_lik_b(1)  = 0;
        log_lik_b(2)  = 10;
        
        % posteriors
        a  = zeros(l,nnegbin);
        a1 = zeros(l,nnegbin);
        b  = zeros(l,1);
        b1 = zeros(l,1);
        
        
        % tol   = 1e-2;
        count = 2;
        % reset tol_nbin and tol_gauss
        tol_nbin = abs(log_lik_a(count,:)-log_lik_a(count-1,:))./abs(log_lik_a(count,:));
        tol_nbin = max(tol_nbin);
        
        % update for tollerance check gaussian
        tol_gauss = abs(log_lik_b(count)-log_lik_b(count-1))/abs(log_lik_b(count)); 
    end
end
% p_tot(i_nan)=1;
