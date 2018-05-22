function Z = SmRG_scrollDataset(V)
% SmRG_scrollDataset: 
%           scrolls dataset until button press.
% 
%           1. This function is called when the user wants to select a
%           slice in the stack where to search for potential seeds.  
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

% check inputs 
if nargin <1
    help SmRG_scrollDataset
    return
end

% get image size
[nx,ny,nz] = size(V);

% first and last slices indexes
minSli = 1;
maxSli = nz;

% set some figure properties
figure(1);
f = gcf;
set(f, 'ToolBar', 'none');

% set some axis properties
a = axes;
axis([1 nx 1 ny]);

% display dataset
imagesc(a,V(:,:,1));
drawnow;

data = struct('val',1);
% set slider behaviour
Sli = uicontrol('Parent',f,'Style','slider',...
    'Min',minSli,'Max',maxSli,...
    'SliderStep',[1 1]./(maxSli-minSli),...
    'Value',1,'Position',[20 5 520 20],...
    'Tag','Sli',...
    'Callback',@Sli_callback);
drawnow;

% set button behaviour
button = uicontrol('Parent', f,'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.91 0.2 0.05 0.1],...
    'String','Z',...
    'Callback',@button_callback);
drawnow;

set( findall( f, '-property', 'Units' ), 'Units', 'Normalized' )

% slider callback
    function Sli_callback(hObject,eventdata,handles)
        data.val = round(get(hObject,'Value'));
        Z=([data.val]);
        imagesc(V(:,:,Z),'Parent',a);
        drawnow;
    end

% button callback
    function button_callback(hObject,eventdata,handles)
        Z=([data.val]);
        close(f)        
    end

% iterate until button press (select Z for seed search)
uiwait(f)
end