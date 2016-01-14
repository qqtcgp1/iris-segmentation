function [ points ] = bw2points( bw )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

index = find(bw);
points = zeros(length(index), 2);
for i = 1:length(index)
    [points(i,1), points(i,2)] = ind2sub(size(bw), index(i));
end

end

