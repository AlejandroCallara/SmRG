function [p_tot,a1,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_newPost(Vin,K0_double,vB_double,mu_sk,rk)
% SmRG_mixtureModelFitting:
%           fits model described in [1] on data input Vin. The model
%           describes a single class k of signal pixels (negative-binomial)
%           and the background (gaussian).For more details refer to [1].
%           For negative binomial parameter update refer to [2].
%
% Syntax:
%           [p_tot,a1]=SmRG_mixtureModelFitting_newPost(vin)
%           [p_tot,a1,K0_double,vB,mu_sk,rk]=SmRG_mixtureModelFitting_newPost(vin)
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


% check for nans in data
vin = Vin(:);
l = length(vin);
% [Vd,i_nan,i_ok]=SmRG_workWithNans(vin);
Vd=vin;
% pre-allocate p_tot
p_tot=zeros(1,l);
a=zeros(1,l);
a1=zeros(1,l);
b=zeros(1,l);
b1=zeros(1,l);

i_ok = 1:length(vin);
V = uint16(Vd(:));
if nargin <2
    % Gaussian initial condition
    K0_double =2*(min(Vd))+5; % mean
    vB_double = 2*K0_double;    % variance
    
    % Neg-Bin initial condition
    mu_sk = mean(Vd); % mean
    v_sk = 1.2*mu_sk; % variance
    rk = mu_sk^2/(v_sk-mu_sk);
end
v_sk = (mu_sk^2+mu_sk*rk)/rk;
pk   = mu_sk/v_sk;

% Mixing prior probability (0.5 0.5)
pb = 0.5; 
pa = 1-pb; % prior

% EM initial conditions
n_of_pixel    = length(Vd);
log_lik_a = zeros(5000,1);
log_lik_b = zeros(5000,1);
log_lik_a(1)  = 0; log_lik_b(1)  = 0;
log_lik_a(2)  = 10;log_lik_b(2)  = 10;

tol   = 1e-4;
count = 2;
% reset tol_nbin and tol_gauss
tol_nbin = abs(log_lik_a(count)-log_lik_a(count-1))/abs(log_lik_a(count));
tol_nbin = max(tol_nbin);

% update for tollerance check gaussian
tol_gauss = abs(log_lik_b(count)-log_lik_b(count-1))/abs(log_lik_b(count));

% USEFUL IN DEBUG MODE: PLOTS THE MODEL FITTING ON DATA HIST
%     figure
%     h = histogram(Vin)
%     hold on;
%     asse_x =1:1000;% linspace(h.BinLimits(1),h.BinLimits(2),2000);
%     hold on
%     plot(asse_x,90000*normpdf(asse_x,K0_double,sqrt(vB)));
%     plot(asse_x,90000*nbinpdf_mu(asse_x,mu_sk,1/rk));
%     drawnow
%% mixture, EM

% iterates until convergence
while (tol_gauss>tol...
        && tol_nbin>tol )
    
    % initialize log likelihood
    psi_a = zeros(n_of_pixel,1);
    psi_b = psi_a;
    
    % initialize log posterior
    b = psi_a;
    a = psi_a;
    
    % E-step
    for ii=1:n_of_pixel
        
        % cast
        % K0_double = double(K0); vB_double = double(vB);
        Vtmp = (Vd(ii));
        
        %  Gaussian log likelihood
        psi_b_tmp = (-SmRG_normlike([K0_double sqrt(vB_double)],Vtmp));
        
        % USEFUL IN DEBUG MODE
        % if isnan(psi_b_tmp);keyboard;end
        
        %if (Vtmp-K0_double)<0;K0_double=double(min(V));end
        psi_b(ii,1) = psi_b_tmp;
        
        % USEFUL IN DEBUG MODE
        %  Neg-Bin log likelihood
        psi_a_tmp = (-nbinlike_mu([mu_sk 1/rk],(Vtmp)));
        
        % USEFUL IN DEBUG MODE
        %  if isnan(psi_a_tmp);keyboard;end
        
        psi_a(ii,1) = psi_a_tmp;
        % find who's bigger between [log(P(x|a))+log(P(a)] e [log(P(x|b))+log(P(b)]
        % to use the logarithm property of sum: log(sum(ai))=log_b(a0)+log_b(1+sum(b^(log_b(ai)-log_b(a0))
        % with a0>a1>...>aN
        % in this way it is possible to evaluate the log posterior without
        % calculating the exponential: this helps for numerical stability
        % (sometimes exponential blows up)
        
        % find the denominator of log posterior
        if (psi_a_tmp+log(pa))>(psi_b_tmp+log(pb))
            den = psi_a_tmp+log(pa)+log(1+exp(psi_b_tmp+log(pb)-psi_a_tmp+log(pa)));
        else
            den = psi_b_tmp+log(pb)+log(1+exp(psi_a_tmp+log(pa)-psi_b_tmp+log(pb)));
        end
        % old way to evaluate the log posterior
        %  b(ii,1) = ((psi_b_tmp)+log(pb))-log(exp(psi_b_tmp)*pb+exp(psi_a_tmp)*pa);
        %  a(ii,1) = ((psi_a_tmp)+log(pa))-log(exp(psi_b_tmp)*pb+exp(psi_a_tmp)*pa);
        
        % new way; actually much better
        b(ii,1) = ((psi_b_tmp)+log(pb))-den;
        a(ii,1) = ((psi_a_tmp)+log(pa))-den;
    end
    b1 = exp(b);
    a1 = exp(a);
    
    % USEFUL IN DEBUG MODE
    %  if sum(isnan(b1))~=0;keyboard;end
    %  if sum(isnan(a1))~=0;keyboard;end
    
    % M-step
    
    % update Gaussian Parameters
    K0_double = sum(b1.*Vd)/sum(b1);
    vB = sum(b1.*((Vd-K0_double).^2))/sum(b1);
    K0 = int16(K0_double);
    vB_double = vB;
    
    % update Neg-Bin Parameters according to [2]
    dk = rk*(psi(rk+Vd)-psi(rk));
    tauk = a1(:);
    lambdak = sum(tauk.*dk)./sum(tauk);
    betak = 1- 1/(1-pk)- 1/(log(pk));
    pk    = (betak*sum(tauk.*dk))/sum(tauk.*(Vd-(1-betak)*dk));
    rk    = -lambdak/log(pk);
    % USEFUL IN DEBUG MODE
    %  if ~isreal(pk)
    %      keyboard
    %  end
    
    v_sk = rk*(1-pk)/(pk^2);
    mu_sk= pk*v_sk;

    % today's posterior becomes tomorrow's prior :)
    pa = mean(a1(:)');
    pb = mean(b1(:)');
    
    % update for tollerance check
    count= count+1;
    log_lik_a(count)=mean(psi_a);
    log_lik_b(count)=mean(psi_b);
    
    % prevent infinite loop
    if count > 4999
        break
    end
    % reset tol_nbin and tol_gauss
    tol_nbin = abs(log_lik_a(count)-log_lik_a(count-1))/abs(log_lik_a(count));
    tol_nbin = max(tol_nbin);
    
    % update for tollerance check gaussian
    log_lik_b(count)=mean(psi_b);
    tol_gauss = abs(log_lik_b(count)-log_lik_b(count-1))/abs(log_lik_b(count));
        
end
count
p_tot(i_ok) =a1;
% p_tot(i_nan)=1;
