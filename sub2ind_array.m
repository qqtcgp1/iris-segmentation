function [ array_1column ] = sub2ind_array( im_size, array_2column )
%sub2ind_array is an extended version to sub2ind. sub2ind only works for 1
%set of subscript, while this version works for an nx2 array, each row
%giving a set of subscript.
%the output is a column vector, each element is the index converted from
%the same row of array_2column.

assert(size(array_2column,2) == 2);

array_1column = zeros(size(array_2column,1),1);
for i = 1:length(array_1column)
    array_1column(i) = sub2ind(im_size, array_2column(i,1), array_2column(i,2));
end
end