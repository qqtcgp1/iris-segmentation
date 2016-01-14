function [ path ] = find_path( bw, first, last, exclude )
%path = find_path( bw, first, last) finds all points that form a connected
%path from first to last in the black and white image bw. path == [] if no
%path is found. first and last can be input as linear index or 2D subscripts.
%The output is given as an array of linear indices. See help sub2ind
%and ind2sub.

% check for validity. first and last must be white points of bw.
if nargin == 3
    exclude = false;
end

if length(first)==2
    first = sub2ind(size(bw), first(1), first(2));
end

if length(last)==2
    last = sub2ind(size(bw), last(1), last(2));
end

assert( bw(first) && bw(last) );

neighbors_first = white_neighbors(bw, first,8);
neighbors_first = [neighbors_first, first];

neighbors_last = white_neighbors(bw,last,8);
neighbors_last = [neighbors_last, last];




% Create a copy of bw, and in the new copy, flip the points first, last to black.
bw_ = bw; bw_(neighbors_first) = 0; bw_(neighbors_last) = 0;

CC = bwconncomp(bw_,8);
[~, I] = sort( cellfun( @length, CC.PixelIdxList ), 'ascend');
CC.PixelIdxList = CC.PixelIdxList (I);

% find the component that contains one point from each of neighbor list.
% Almost there....

add_back = [first; last];

if exclude
    add_back = [];
end

for i = 1:CC.NumObjects
    pixel_list = CC.PixelIdxList{i};
    flag1 = 0; flag2 = 0;
    
    for j = 1:length(pixel_list)
        [sub1, sub2] = ind2sub(size(bw), pixel_list(j));
        
        for k = 1:length(neighbors_first)
            [Sub1,Sub2] = ind2sub(size(bw),neighbors_first(k));
            if  (sub1 == Sub1 || sub1+1 == Sub1 || sub1 == Sub1+1) && ...
                (sub2 == Sub2 || sub2+1 == Sub2 || sub2 == Sub2+1)
            add_back = [add_back; neighbors_first(k)];
            flag1 = 1; break;
            end
        end
        
        % if flag1; continue; end;
        
        for k = 1:length(neighbors_last)
            [Sub1,Sub2] = ind2sub(size(bw),neighbors_last(k));
            if  (sub1 == Sub1 || sub1+1 == Sub1 || sub1 == Sub1+1) && ...
                (sub2 == Sub2 || sub2+1 == Sub2 || sub2 == Sub2+1)
            add_back = [add_back; neighbors_last(k)];
            flag2 = 1; break;
            end
        end
        
    end
    
    if flag1 && flag2; 
        object_label = i; break; 
    end;
end

if flag1 && flag2
    path = [CC.PixelIdxList{object_label}; add_back];
else path = [];
end

end