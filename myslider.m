function [  ] = myslider( command_str,callback,position, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



if strcmp(command_str, 'initialize')
    temp = get(findobj('tag', 'slider'), 'parent');
    if ~isempty(temp)
        delete(temp(:));
    end
    
    if nargin < 3
        position = [800 400 100 40];
    end
    
    if nargin < 2
        callback = '';
    end
    
    h = uicontrol('style', 'slider', 'tag', 'slider',...
        'unit', 'pix', 'position', [position(1), position(2) + position(4)*.5, position(3), position(4)*.5], ...
         varargin{:}, ...
        'callback', ['myslider(''edit_slider'');', callback]);
    slider_text(get(h,'value'), callback,position);
    
elseif strcmp(command_str, 'edit_slider')
    set(findobj('tag', 'edit'), 'string', num2str(get(findobj('tag','slider'),'value')));
    
elseif strcmp(command_str, 'edit_text')
    new_data = str2double(get(findobj('tag', 'edit'), 'string'));
    if isnumeric(new_data) && (~isnan(new_data))
        if new_data > get(findobj('tag', 'slider'), 'max')
            new_data = get(findobj('tag','slider'),'max');
        elseif new_data < get(findobj('tag', 'slider'), 'min')
            new_data = get(findobj('tag','slider'),'min');
        end
        
        new_data = nearest(new_data);
        set(findobj('tag', 'slider'), 'value', new_data);
        set(findobj('tag', 'edit'), 'string', num2str(new_data));
        
    else
        set(findobj('tag', 'edit'), 'string', num2str(get(findobj('tag', 'slider'), 'value')));
    end
    
end


end


function [ ] = slider_text( value,callback,position)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

uicontrol('style', 'edit', 'string', num2str(value),...
    'tag', 'edit', ...
    'unit', 'pix', 'position', [position(1), position(2), position(3), position(4)*.5], ...
    'callback', ['myslider(''edit_text'');', callback]);

end