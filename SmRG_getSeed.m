function [seed]=SmRG_getSeed(V)
% SmRG_getSeed: 
%           allows to select a seed by manually clicking on the image.
% 
%           1. This function is called when the user wants to select a
%           slicE in the stack where to search for potential seeds.  
%           2. otherwise it can be used just for visualization of a 3D
%           stack.
%
% Syntax:
%           Z = SmRG_scrollDataset(V)
%
% Input:
%           V: 3D grayscale/logic image
% Output:
%           Z: position of the slider

% get image size
[nx,ny,nz]=size(V);

% first and last slices indexes
minSli = 1;
maxSli = nz;

% set some figure properties
figure(1);
f = gcf;
set(f, 'ToolBar', 'none');
title('Select seed and press Enter')
% set some axis properties
a = axes; 
axis([1 nx 1 ny]);

% display dataset
imagesc(a,V(:,:,1));

% set slider behaviour
Sli = uicontrol('Style','slider','Min',minSli,'Max',maxSli,...
                'SliderStep',[1 1]./(maxSli-minSli),'Value',1,...
                'Position',[20 20 200 20]);
set(Sli,'Callback',@(hObject,eventdata) imagesc(V(:,:,round(get(hObject,'Value'))),'Parent',a) )
% initialize seed with empty vector
seed = [];

% iterate until a seed is selected
while isempty(seed)
    try
        [x,y] = getpts(gca,1);
        z = round(Sli.Value);
        seed  = [round(x),round(y),z];
        close(f)
        disp('Seed with coordinates: ')
        disp(num2str(seed))
        disp('selected as starting point')
    catch ME
        if (strcmp(ME.identifier,'images:getpts:interruptedMouseSelection'))&&(~isvalid(f))
            disp('Seed not selected by user')
            break 
        end
    end
end

