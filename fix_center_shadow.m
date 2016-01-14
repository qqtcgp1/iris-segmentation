function [ frame, rect_position1, rect_position2 ] = fix_center_shadow( rect_position1, rect_position2, original, graythresh, disksize )
%fix_center_shadow is a substitute of gray2boundary. added functionality of
%removing the center shadow..
%if the position of the shadow rect_position is given (as [left, bottom,
%width, height] format, the shadow will be dealt with silently. If this
%argument is given as an empty array, the function asks the user to crop
%the rectangle of the position in an image. The cropped position will be
%output by the function, fercilitating the processing of other similar
%images (similar as the position of the shadow).

if nargin < 4
    disksize = 2;
end
frame = original>=graythresh;


% if rect_position is not provided, ask user to manually specify the rectangle.
if isempty(rect_position1)
    fig = figure('tag', 'xvyei', 'userd', {original; frame});
    imshow(frame);
    myslider('initialize', 'callback_1(''xvyei'')', [200 5 100 40], 'units', 'pix', ...
        'min', 1, 'max', 255, 'value', graythresh, 'sliderstep', [1/254 20/254]);
    uicontrol('style', 'push', 'string', 'graythresh ok', 'unit', 'pix', 'posi', [300 10 100 20], 'callb', 'uiresume(gcbf)');
    uiwait(gcf);
    h = uicontrol('style', 'text', 'unit', 'pix', 'position', [5 5 100 40], 'string', 'Crop the center shadow');
    temp = imrect;
    rect_position1 = wait(temp);
    rect_position1 = arrayfun(@nearest, rect_position1);
    set(h, 'string', 'Crop (part of) the iris');
    temp = imrect;
    rect_position2 = wait(temp);
    rect_position2 = arrayfun(@nearest, rect_position2);
    frame = get(fig, 'userd');
    frame = frame{2};
end


% erase the entire rectangle to black.
frame(:, rect_position1(1):rect_position1(1) + rect_position1(3)) = 0;

frame = larger_comp(frame, 0.2);
frame = imfill(frame, 'holes');

% make a new copy
frame_copy = frame;
frame_copy(rect_position2(2):rect_position2(2) + rect_position2(4), rect_position2(1):rect_position2(1) + rect_position2(3)) = 0;

% get rid of the white and black islands.
% find the erased part of the cornea (middle part of cornea).
corners = corner(frame_copy, 'MinimumEigenvalue');

[~,I] = min(abs(corners(:,1) - rect_position1(1)));
left_corners(1,:) = corners(I,[2 1]);
corners(I,:) = [];
[~,I] = min(abs(corners(:,1) - rect_position1(1)));
left_corners(2,:) = corners(I,[2 1]);
corners(I,:) = [];
if left_corners(1,1) > left_corners(2,1)
    left_corners = left_corners([2,1],:);
end


[~,I] = min(abs(corners(:,1) - rect_position1(1) - rect_position1(3)));
right_corners(1,:) = corners(I,[2 1]);
corners(I,:) = [];
[~,I] = min(abs(corners(:,1) - rect_position1(1) - rect_position1(3)));
right_corners(2,:) = corners(I,[2 1]);
corners(I,:) = [];

if right_corners(1,1) > right_corners(2,1)
    right_corners = right_corners([2,1],:);
end



% fixed the broken cornea.
frame(min(left_corners(1,1), right_corners(1,1)): max(left_corners(2,1), right_corners(2,1)), min(left_corners(:,2)): max(right_corners(:,2))) = 1;

% do the morphology image processing, if argument disksize is provided.
if nargin == 4
    se = strel('disk', disksize);
    frame = imclose(frame, se);
    frame = ~imclose(~frame, se);
end


% then canny, to get boundary
frame = bwmorph(frame, 'remove');
frame = bwmorph(frame, 'fill');
frame = bwmorph(frame, 'close');
frame = bwmorph(frame, 'thin', Inf);

% branch = bwmorph(frame, 'branch');
% frame = bwmorph(frame, 'skel', Inf);


frame = remove_endpoints(frame);

% flag = 1; N = 20;
% while flag
%     endpoints = bwmorph(frame, 'endpoints');
%     frame(endpoints) = 0;
%     flag = N<= 20 && (sum(sum(endpoints))~= 0);
% end

end



