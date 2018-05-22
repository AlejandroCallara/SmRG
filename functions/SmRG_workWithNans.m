function [Vd,i_nan,i_ok]=SmRG_workWithNans(vin)
% SmRG_workWithNans: 
%           returns the indices of Nan entries in vin
%
% Syntax:
%           [Vd,i_nan,i_ok]=SmRG_workWithNans(vin)
%
% Inputs:
%           vin: 1D array of double
%
% Outputs:
%           Vd:     1D array (size(Vd)=size(vin)). Vd is vin without nan entries. 
%           i_nan:  indices of nan entries in vin.
%           i_ok:   indices of non-nan entries in vin

% check inputs 
if nargin <1
    help SmRG_workWithNans
    return
end

vin = vin(:);
i_ok=1:length(vin);
i_nan=find(isnan(vin));
i_ok(i_nan)=[];
Vd=vin(i_ok);
