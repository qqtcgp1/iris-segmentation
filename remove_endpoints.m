function [ bw ] = remove_endpoints( bw )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

points = bw2points(bw);
mark = [0 0];

for i = 1:size(points,1)
    if length(white_neighbors(bw, points(i,:), 8)) <= 1
        mark = [mark; points(i,:)];
    end
end

if size(mark,1) >= 2
    for i = 2:size(mark,1)
        bw( mark(i,1), mark(i,2)) = 0;
    end
    bw = remove_endpoints(bw);

else
    return;
end
    

end

