function [ processed_bw ] = larger_comp( bw,fraction )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

processed_bw = bw; processed_bw(:,:) = 0;

cc = bwconncomp(bw, 4);

M = max(cellfun(@length, cc.PixelIdxList));

index = find(cellfun(@length, cc.PixelIdxList) >= M * fraction);

points = [];

for i = 1:length(index)
     points = [points; cc.PixelIdxList{index(i)}];
end

processed_bw(points) = 1;    



end

