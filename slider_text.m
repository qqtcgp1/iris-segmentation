function [  ] = slider_text( command_str )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    command_str = 'initialize';
end

if strcmp(command_str, 'initialize')
    uicontrol('style', 'edit', 'value', 22600, 'tag', 'edit', ...
        'unit', 'pix', 'position', [800 380 100 20]);
elseif strcmp(command_str, 'edit')
    set(findobj('tag', 'slider'), 'value', str2num(get(findobj('tag', 'edit'), 'string')));
end

end