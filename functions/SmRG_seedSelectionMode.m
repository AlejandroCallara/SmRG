function str_mode = SmRG_seedSelectionMode
% SmRG_seedSelectionMode: 
%           interactive figure to chose seed selection mode.
% 
% Syntax:
%           str_out = SmRG_seedSelectionMode
%
% Output:
%           str_out: can be one of the following:
%                                                'manual'
%                                                'semi-automatic'
%                                                'automatic'
str_mode = 'none';
string1 = 'Manual Seed';
string2 = 'Semi-automatic Seed';
string3 = 'Automatic Seed ';

% create figure and set some properties
f= figure;
set (f, 'ToolBar', 'none');
set (f, 'Units', 'normalized')
set (f, 'Position', [0.4 0.2 0.2 0.5]);


% set button behaviour
button1 = uicontrol('Parent', f,'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.2 0.7 0.6 0.2],...
    'String',string1,...
    'FontSize',15,...
    'Callback',@button1_callback);
drawnow;

% set button behaviour
button2 = uicontrol('Parent', f,'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.2 0.4 0.6 0.2],...
    'String',string2,...
    'FontSize',15,...
    'Callback',@button2_callback);
drawnow;

% set button behaviour
button3 = uicontrol('Parent', f,'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.2 0.1 0.6 0.2],...
    'String',string3,...
    'FontSize',15,...
    'Callback',@button3_callback);
drawnow;

% button1 callback
    function button1_callback(hObject,eventdata,handles)
        str_mode='manual';
        close(f)        
    end

% button2 callback
    function button2_callback(hObject,eventdata,handles)
        str_mode='semi-automatic';
        close(f)       
    end

% button3 callback
    function button3_callback(hObject,eventdata,handles)
        str_mode='automatic';
        close(f)        
    end

% wait for button press
uiwait(f)

end