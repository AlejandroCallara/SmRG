function y = nbinpdf_mu(x,mu,alpha)
%|====================================================================================
%|NBINPDF_MU Negative binomial probability density function.
%|
%|   Y = NBINPDF_MU(X,MU,ALPHA) returns the negative binomial probability density 
%|                              function with parameters R and P at the values in X.
%|                              Note that the density function is zero unless X is 
%|                              an integer.
%|
%|   NOTE:
%|   The size of Y is the common size of the input arguments. A scalar input  
%|   functions as a constant matrix of the same size as the other inputs.    
%|
%|  Last revision:
%|  22 May 2018
%|  Michele Scipioni, University of Pisa
%|
%|====================================================================================


if nargin < 3
    error(message('stats:nbinpdf:TooFewInputs'));
end

[errorcode, x, mu, alpha] = distchck(3,x,mu,alpha);

if errorcode > 0
    error(message('stats:nbinpdf:InputSizeMismatch'));
else
    r = 1./alpha;
    p = r./(r+mu);    
end

% Initialize Y to zero.
if isa(x,'single') || isa(r,'single') || isa(p,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end
if ~isfloat(x)
    x = double(x);
end

% Out of range or missing parameters and missing data return NaN.
% Infinite values for R correspond to a Poisson, but its mean cannot
% be determined from the (R,P) parametrization.
nans = ~(0 < r & isfinite(r) & 0 < p & p <= 1) | isnan(x);
y(nans) = NaN;

% Negative binomial distribution is defined on the non-negative
% integers.  Data outside this support return 0.
k = find(0 <= x & isfinite(x) & x == round(x)  &  ~nans);
if ~isempty(k)
    lognk = gammaln(r(k) + x(k)) - gammaln(x(k) + 1) - gammaln(r(k));
    y(k) = exp(lognk + r(k).*log(p(k)) + x(k).*log1p(-p(k)));
end

% Fix up the degenerate case.
k = find(p == 1  &  ~nans);
if ~isempty(k)
    y(k) = (x(k) == 0);
end
