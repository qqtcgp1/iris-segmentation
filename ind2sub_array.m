function [ array_2column ] = ind2sub_array( im_size, array_1column )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

array_2column = zeros( length(array_1column),2);

for i = 1: length(array_1column)
    [array_2column(i,1), array_2column(i,2)] = ind2sub(im_size, array_1column(i));
end

end

