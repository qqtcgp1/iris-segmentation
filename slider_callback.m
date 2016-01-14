function [ boundary, area ] = slider_callback( normalized, original )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[boundary, ~] = gray2boundary(normalized, ...
    get(findobj('tag', 'slider'),'value'));

set(findobj('type', 'axes'), 'next', 'remove');
imshow(uint8(original)); set(findobj('type', 'axes'), 'next', 'add');
p = bw2points(boundary); scatter(p(:,2), p(:,1), 'w.');
% [iris_pts, area, volume] = get_iris(boundary);
iris_pts = p;


solid = (normalized < get(findobj('tag', 'slider'), 'value'));

% remove white and black holes
solid = imfill(solid,'holes');
solid = ~imfill(~solid,'holes');

se = strel('disk', 2);
solid = imclose(solid, se);
solid = ~imclose(~solid, se);



area = sum(sum((~solid)));
volume = iris_volume( ~solid, 339.5 );
scatter(iris_pts(:,2), iris_pts(:,1), '.r');
set(findobj('tag', 'text_area'), 'string', char('Area',num2str(area)));
set(findobj('tag', 'text_volume'), 'string', char('Volume', num2str(volume)));

end