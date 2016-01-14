function [ frame ] = remove_markings_Dewarp( frame, width )
%REMOVE_MARKINGS_DEWARP Summary of this function goes here
%   Detailed explanation goes here


index = (frame == 255);
ref(1:size(frame,1), 1:size(frame,2)) = false;
ref(index) = true;

cc = bwconncomp(ref);
[m, I] = max(cellfun(@length, cc.PixelIdxList));
marking_ = cc.PixelIdxList{I};
marking = zeros(m,2);


for i = 1:m
    [marking(i,1), marking(i,2)] = ind2sub(size(frame), marking_(i));
end



position(1) = min(marking(:,1)) - width;

position(4) = max(marking(:,2)) + width;
position(2) = min(marking(:,2)) - width;
position(3) = max(marking(:,1)) + width;

frame(position(1):position(3), position(2):position(4)) = 0;




end

