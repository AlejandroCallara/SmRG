function [background, Vb] = SmRG_getBackground(Vin,method)
% SmRG_getBackground: 
%           extracts background pixels from image Vin with specified method. 
%
% Syntax:
%           [background]     = SmRG_getBackground(Vin)
%           [background, Vb] = SmRG_getBackground(Vin)
%           [background, Vb] = SmRG_getBackground(Vin,method)
%
% Inputs:
%           Vin: 2D image 
%           method: one of {'otsu','mm'}
% Outputs:
%           background:  mean value of background pixels
%           Vb: image with only background pixels

% check inputs 
if nargin <1
    help SmRG_getBackground
    return
end


if nargin <2
    method = 'otsu';
end


% extract background 
if strcmp(method,'otsu')
    
    % Otsu's method
    vv = (Vin-min(Vin(:)))/(max(Vin(:))-min(Vin(:)));
    t = graythresh(vv);
    Vbin = im2bw(vv,t);
    
    % isolate background
    Vb = Vin.*(~Vbin);
    [Vnotnan,i_nan,i_ok]=SmRG_workWithNans(Vb);
    background = round(sum(Vb(:),'omitnan')/length(find(Vnotnan)));
    
elseif strcmp(method,'mm')
     
    % mixture model fitting
    [p_tot,a1] = SmRG_mixtureModelFitting(Vin);
    
    
    % isolate background
    i_background = ((1-p_tot)<=0.95);
    Vbin = false(size(Vin));
    [ii, jj] = ind2sub(size(Vin),find(i_background));
    for kk = 1:length(ii)
        Vbin(ii(kk),jj(kk))=true;
    end 
    
    % outputs
    Vb = Vin.*(~Vbin);
    background = round(sum(Vb(:))/sum(~Vbin(:)));
end