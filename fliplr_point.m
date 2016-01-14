function [ point_out ] = fliplr_point( point, size2 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
point_out = [point(1), size2+1 - point(2)];

end

