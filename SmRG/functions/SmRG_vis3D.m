function SmRG_vis3D(J)
% SmRG_vis3D: 
%           3D plot of binary 3D dataset. (ideally the output of
%           SmRG_regionGrowing
%
% Syntax:
%           SmRG_getBackground(J)
%
% Inputs:
%           Vin: 2D image 
%           method: one of {'otsu','mm'}
% Outputs:
%           background:  mean value of background pixels
%           Vb: image with only background pixels

% check inputs 
if nargin <1
    help SmRG_vis3D
    return
end

[p,v] = isosurface(J, 0.5);
figure;
patch('Faces',p, 'Vertices',v,'FaceColor','green',...
      'EdgeColor','none','FaceAlpha',1);
axis equal;

camlight; 
camlight(0, -100); 
lighting phong;
disp('plot finished')