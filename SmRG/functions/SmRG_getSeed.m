function [seed]=SmRG_getSeed(V)
% SmRG_getSeed: 
%           allows to select a seed by manually clicking on the image.
%
% Syntax:
%           seed = SmRG_getSeed(V)
%
% Input:
%           V: 3D grayscale/logic image
% Output:
%           seed: seed clicked by the user

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
axis([1 nx 1 ny]);

% display dataset
imagesc(V(:,:,1));

% set slider behaviour
Sli = uicontrol('Parent',f,'Style','slider',...
    'Min',minSli,'Max',maxSli,...
    'SliderStep',[1 1]./(maxSli-minSli),...
    'Value',1,'Position',[20 5 520 20],...
    'Tag','Sli',...
    'Callback',@Sli_callback);
drawnow;

% initialize seed with empty vector
seed = [];


% slider callback
    function Sli_callback(hObject,eventdata,handles)
        z = round(get(hObject,'Value'));
        imagesc(V(:,:,z),'Parent',gca);
        drawnow;
    end

set( findall( f, '-property', 'Units' ), 'Units', 'Normalized' )

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
end
