function [ frame ] = remove_markings( frame )

% 5.remove the mark on top-right
frame(50:72,603:638) = 0;


for i = 1:39
    frame = outlyer_average(frame, [78,305 + (i-1)*2], 8);
end

for i = 1:37
    frame = outlyer_average(frame, [50,307 + (i-1)*2], 8);
end

for i = 1:14
    frame = outlyer_average(frame, [22 + (i-1)*2, 343], 8);
end

for i = 50:77
    frame = outlyer_average(frame, [i, 305]);
    frame = outlyer_average(frame, [i, 381]);
end


end




function [ bw ] = outlyer_average( bw, point, connectivity )
%outlyer_average addresses an outlyer point of an bw image by averaging it
%out with neighboring points. disregards the neighboring points that has the same intensity value as the
%interested point.

if nargin == 2
    connectivity = 8;
end

if length(point) == 1
    [point(1), point(2)] = ind2sub(size(bw), point);
end



% find all neighbors according to connectivity
neighbors = find_neighbors(size(bw), point, connectivity);

% exlude points that has the same intensity as the center point
neighbors( bw(neighbors) == bw(point(1), point(2))  ) = [];

bw(point(1), point(2)) = mean(bw(neighbors));


end