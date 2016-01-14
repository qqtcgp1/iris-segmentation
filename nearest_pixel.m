function [ output_index ] = nearest_pixel( bw, point )
%nearest_pixel(bw, point) finds the closest white point to point in the
%black and white image bw. point can be in subscript and index. Output is
%in linear index.

if length(point) == 2
    point = sub2ind( size(bw), point(1), point(2) );
end

if bw(point) == 1
    output_index = point; return;
end

white_pixels = find( bw==1 );

[x,y] = ind2sub(size(bw), point);
dist = zeros(length(white_pixels),1);

for i = 1:length(dist)
    [x1,y1] = ind2sub(size(bw),white_pixels(i));
    dist(i) = (x1 - x).^2 + (y1-y).^2;
end

[~, temp] = min(dist);
output_index = white_pixels(temp);

end

