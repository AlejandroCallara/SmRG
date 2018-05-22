function SmRG_visSeed(V,S)
% SmRG_visSeed : 
%           void function. overlaps two 3D datasets for visualization. 
%
% Syntax:
%           SmRG_visSeed(V,S)
%
% Inputs:
%           V: 3D grayscale/logic image
%           S: 3D grayscale/logic image

% check inputs 
if nargin <2
    help SmRG_visSeed
    return
end

% get image size
[nx,ny,nz]=size(V);

% first and last slices indexes
minSli = 1;
maxSli = nz;

% set some figure properties
f = figure;
gcf;
set(f, 'ToolBar', 'none');

% set some axis properties
a = axes; 
axis([1 nx 1 ny]);

% create slider
Sli = uicontrol('Parent',f,'Style','slider',...
    'Min',minSli,'Max',maxSli,...
    'SliderStep',[1 1]./(maxSli-minSli),...
    'Value',1,'Position',[20 5 520 20],...
    'Tag','Sli',...
    'Callback',@Sli_callback);
drawnow;

set( findall( f, '-property', 'Units' ), 'Units', 'Normalized' )

% display dataset
imshowpair(V(:,:,1),S(:,:,1),'Parent',a);

% slider callback
    function Sli_callback(hObject,eventdata,handles)
        z = round(get(hObject,'Value'));
        imshowpair(V(:,:,z),S(:,:,z),'Parent',a);
        drawnow;
    end

while isvalid(f)
    uiwait(f)
end
end