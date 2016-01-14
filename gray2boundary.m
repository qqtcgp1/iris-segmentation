function [ boundary, graythresh ] = gray2boundary( gray, graythresh, disksize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin <= 2
    disksize = 2;
end

if isempty(graythresh)
    graythresh = 30;
    boundary = gray >= graythresh;
    figure('tag', '123haha', 'userdata', {gray; boundary});
    imshow(boundary);
    
    myslider('initialize', 'callback_1(''123haha'')', [200 5 100 40], 'units', 'pix', ...
        'min', 1, 'max', 255, 'value', graythresh, 'sliderstep', [1/254 20/254]);
    uicontrol('style', 'push', 'string', 'graythresh ok', 'unit', 'pix', 'posi', [300 10 100 20], 'callb', 'uiresume(gcbf)');
    uiwait(gcf);
    graythresh = get(findobj('tag', 'slider'), 'value');
    boundary = get(findobj('tag', '123haha'), 'userdata');
    boundary = boundary{2};
else
    boundary = (gray >= graythresh);
end


se = strel('disk', disksize);
boundary = imclose(boundary, se);
boundary = ~imclose(~boundary, se);


% remove white and black holes
boundary = imfill(boundary,'holes');
boundary = larger_comp(boundary, 0.3);

boundary = bwmorph(boundary, 'remove');
boundary = bwmorph(boundary, 'fill');
boundary = bwmorph(boundary, 'close');
boundary = bwmorph(boundary, 'thin', Inf);

end

