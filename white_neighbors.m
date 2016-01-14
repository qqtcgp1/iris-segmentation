function [ neighbor_set ] = white_neighbors( bw, point, connectivity )
%neightbor_set = white_neighbors(bw, point, connectivity) finds all bright
%points neighbor to point in the black and white image bw. Uses the
%function find_neighbors. Connectivity is defaulted to be 8.
%point can be input as linear index or 2D subscripts.
%The output is given as an array of linear indices. See help sub2ind
%and ind2sub.

if nargin == 3
    assert(connectivity == 8 || connectivity == 4);
elseif nargin == 2
    connectivity = 8;
elseif nargin < 2
        error('Too few arguments.\n');
end

neighbor_set_ = find_neighbors(size(bw), point, connectivity);

neighbor_set = [];

for i = 1:length(neighbor_set_)
    if bw(neighbor_set_(i))
        neighbor_set = [neighbor_set, neighbor_set_(i)];
    end
end


end