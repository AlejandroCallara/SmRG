function [p_tot,a1,K0_double,vB,mu_sk,rk] = SmRG_mixtureModelFitting(Vin)
% SmRG_mixtureModelFitting: 
%           fits model described in [1] on data input Vin. The model 
%           describes a single class k of signal pixels (negative-binomial)
%           and the background (gaussian).For more details refer to [1].
%
% Syntax:
%           [p_tot,a1]=SmRG_mixtureModelFitting(vin)
%           [p_tot,a1,K0_double,vB,mu_sk,rk]=SmRG_mixtureModelFitting(vin)
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

% check inputs
if nargin<1
    help SmRG_mixtureModelFitting
    return
end

if sum(Vin(:)<0)
    error ('only images with positive values are allowed')
end

% check for nans in data
vin = Vin(:); 
[Vd,i_nan,i_ok]=SmRG_workWithNans(vin);

% pre-allocate p_tot
p_tot=zeros(1,length(vin));

V = uint16(Vd(:));

% Gaussian initial condition
K0 =2*(min(V)); % mean
vB = 2*K0;    % variance

% Neg-Bin initial condition
mu_sk = 2*mean(Vd); % mean
v_sk = var(Vd);     % variance

rk = mu_sk^2/(v_sk-mu_sk);

% Mixing prior probability (0.5 0.5)
pb = 0.5; pa = 1-pb; % prior

% EM initial conditions
n_of_pixel    = length(Vd);
log_lik_a(1)  = 0; log_lik_b(1)  = 0;
log_lik_a(2)  = 10;log_lik_b(2)  = 10;

tol   = 1e-6; 
count = 2;

%% mixture, EM

% iterates until convergence
while (abs(log_lik_a(count)-log_lik_a(count-1))>tol...
        && abs(log_lik_b(count)-log_lik_b(count-1))>tol )
    
    % initialize log likelihood
    psi_a = zeros(n_of_pixel,1);
    psi_b = psi_a;
    
    % initialize log posterior
    b = psi_a;
    a = psi_a;
    
    % E-step
    for ii=1:n_of_pixel
        
        % cast
        K0_double = double(K0); vB_double = double(vB);
        Vtmp = (Vd(ii));
        
        %  Gaussian log likelihood 
        psi_b_tmp = (-normlike([K0_double sqrt(vB_double)],Vtmp));
            if isnan(psi_b_tmp);keyboard;end
            %if (Vtmp-K0_double)<0;K0_double=double(min(V));end
        psi_b(ii,1) = psi_b_tmp;
        
        %  Neg-Bin log likelihood 
        psi_a_tmp = (-nbinlike_mu([mu_sk 1/rk],(Vtmp)));
            if isnan(psi_a_tmp);keyboard;end
        psi_a(ii,1) = psi_a_tmp;
        
        % log posterior
        b(ii,1) = ((psi_b_tmp)+log(pb))-log(exp(psi_b_tmp)*pb+exp(psi_a_tmp)*pa);
        a(ii,1) = ((psi_a_tmp)+log(pa))-log(exp(psi_b_tmp)*pb+exp(psi_a_tmp)*pa);
        
    end
    b1 = exp(b);
    a1 = exp(a);
    if sum(isnan(b1))~=0;keyboard;end
    if sum(isnan(a1))~=0;keyboard;end
    
    % M-step
    
    % update Gaussian Parameters
    K0_double = sum(b1.*Vd)/sum(b1);
    vB = sum(b1.*((Vd-K0_double).^2))/sum(b1);
    K0 = int16(K0_double);
    
    % update Neg-Bin Parameters
    mu_sk = sum(a1.*(Vd))/sum(a1);
    v_sk = sum(a1.*(Vd-mu_sk).^2)/sum(a1);%   pk = mu_sk/v_sk;
    
    rk = mu_sk^2/(v_sk-mu_sk);
    
    % today's posterior becomes tomorrow's prior :)
    pb = mean(b1);
    pa = mean(a1);
    
    % update for tollerance check
    count= count+1;
    log_lik_a(count)=mean(psi_a);
    log_lik_b(count)=mean(psi_b);
end

p_tot(i_ok) =a1;
p_tot(i_nan)=1;   