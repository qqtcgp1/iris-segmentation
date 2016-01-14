function [ out ] = fliplr_points( points, imsize )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


imwidth = imsize(2);
out = points;
out(:,2) =  (imwidth + 1)*ones( size(points,1),1) -points(:,2) ;

end