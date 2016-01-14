function [ neighbor_set ] = find_neighbors( size, point, connectivity )
%find_neighbors collects all points that are neighbors to the input point.
%connectivity is defaulted to be 8. size represents the dimension of the
%image. point can be input as linear index or 2D subscripts.
%The output is given as an array of linear indices. See help sub2ind
%and ind2sub.

if nargin == 3
    assert(connectivity == 8 || connectivity == 4);
elseif nargin == 2
    connectivity = 8;
elseif nargin < 2
        error('Too few arguments.\n');
end

neighbor_set = [];
size1 = size(1); size2 = size(2);

if length(point)==1
    [I,J] = ind2sub(size, point);
    point = [I J];
end

if point(1)~=1
        neighbor_set = [neighbor_set, sub2ind(size, point(1)-1,point(2)) ];
    
    if connectivity == 8
        if point(2)~=1
            neighbor_set = [neighbor_set, sub2ind(size,point(1)-1,point(2)-1)];
        end
        
        if point(2)~=size2
            neighbor_set = [neighbor_set, sub2ind(size,point(1)-1, point(2)+1)];
        end
    end
end

if (point(2)~=1)
    neighbor_set = [neighbor_set, sub2ind(size,point(1), point(2)-1)];
end

if (point(2)~=size2)
    neighbor_set = [neighbor_set, sub2ind(size,point(1), point(2)+1)];
end


if point(1)~=size1
        neighbor_set = [neighbor_set, sub2ind(size,point(1)+1,point(2))];
    
    if connectivity == 8
        if point(2)~=1
            neighbor_set = [neighbor_set, sub2ind(size,point(1)+1,point(2)-1)];
        end
        
        if point(2)~=size2
            neighbor_set = [neighbor_set, sub2ind(size,point(1)+1, point(2)+1)];
        end
    end
end

end