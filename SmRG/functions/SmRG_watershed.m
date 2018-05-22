function [Jwat] = SmRG_watershed(J)
% SmRG_watershed: 
%           applies watershed on 3D binary image.
% 
% Syntax:
%           [outs] = SmRG_watershed(J)
%
% Input:
%           J: 3D binary image
% Output:
%           Jwat: watershed 

% check input
if nargin<1
    help SmRG_watershed
    return
end

% get dims
[nx, ny, nz] = size(J);

% applies watershed for each slice
for d = 1:nz
    wat_tmp  = J(:,:,d);
    Jdist = bwdist(~wat_tmp);
    tymp = -Jdist;
    mask1= imextendedmin(tymp,1);
    mask2= imimposemin(tymp,(mask1));
    waters = watershed(mask2);
    wat_tmp(waters==0)=0;
    Jwat(:,:,d)=wat_tmp;
end