function SmRG_vis3D(J,tree,alfa)
% SmRG_vis3D: 
%           3D plot of binary 3D dataset (ideally the output of SmRG_regionGrowing)
%           with optional swc superimposed. 
%           .
%
% Syntax:
%           SmRG_vis3D(J)
%           SmRG_vis3D(J,tree)
%
% Inputs:
%           J: 3D binary dataset 
%           tree = swc Nx7 matrix
% Outputs:
%           background:  mean value of background pixels
%           Vb: image with only background pixels

% check inputs 
if nargin <1
    help SmRG_vis3D
    return
end
figure;
if nargin ==2
    color = [1 0 0];
    hold on
    for i = 1 : size(tree, 1)
        % Draw a line between the current node and its parent
        cnode = tree(i, 3:5);
        parent = tree(i, 7);
        [pidx] = find(tree(:, 1) == parent);
        if numel(pidx) == 1
            pnode = tree(pidx, 3:5);
            h = plot3([cnode(1);pnode(1)], [cnode(2);pnode(2)], [cnode(3); pnode(3)],...
                'Color', color);
            set(h(1),'linewidth',3);
        end
    end
    %The following line of code can adjust axis ratio
    daspect([1 1 1])
    alpha 0
end
[p,v] = isosurface(J, 0.5);

patch('Faces',p, 'Vertices',v,'FaceColor','green',...
      'EdgeColor','none','FaceAlpha',0.5);
axis equal;

camlight; 
camlight(0, -100); 
lighting phong;
disp('plot finished')
