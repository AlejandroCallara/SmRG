function boot_dip = SmRG_nullOfUnifDist(nboot)
% SmRG_nullOfUnifDist: 
%           creates the null distribution for the hartigans'dip test 
%           sampling nboot times. 
%
% Syntax:
%           boot_dip=SmRG_nullOfUnifDist(nboot)     {nboot = 10000}
%
% Input:
%           nboot: Number of samples to sample from uniform distribution
%
% Output:
%           boot_dip: bootstrapped sample of uniform distribution

% check inputs
if nargin<1
    nboot = 10000;
end

boot_dip =[];

for iboot=1:nboot
    unifpdfboot=sort(unifrnd(0,1,1,255));
    [unif_dip]=HartigansDipTest(unifpdfboot);
    boot_dip=[boot_dip; unif_dip];
end

boot_dip=sort(boot_dip);
