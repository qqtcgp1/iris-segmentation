function [ result ] = polygon_overlap_fix( polygon1, polygon2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[intersectx, intersecty] = polybool('intersection', ...
    polygon1(:,1), polygon1(:,2), polygon2(:,1), polygon2(:,2));

if isempty(intersectx)
    result = 0;
    return;
end

result = 2.0* polyarea(intersectx, intersecty) ./ (polyarea( polygon1(:,1), polygon1(:,2) ) + ...
    polyarea( polygon2(:,1), polygon2(:,2)));


end