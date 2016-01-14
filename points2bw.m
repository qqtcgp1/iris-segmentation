function [ bw ] = points2bw( points, im_size )
%points2bw converts an nx2 array specifying the white pixels in an image to
%a black and white image.

if size(points, 2) == 2
    points_ = zeros(size(points,1),1);
    for i = 1:length(points_)
        points_(i) = sub2ind(im_size, points(i,1), points(i,2));
    end
end

bw = zeros(im_size(1), im_size(2));
bw(points_) = 1;


end