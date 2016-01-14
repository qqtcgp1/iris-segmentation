function [ integer ] = nearest_( x, upper_bound )
% a modified version of nearest(x). Assumes a lower bound of 0, and gives
% an upper bound. Outputs 1 if input <= 0, outputs upper_bound - 1 if input
% >= upper_bound.

if numel(x) > 1
    integer = arrayfun( @ nearest_, x, upper_bound*ones(size(x)) ); return;
end

integer = nearest(x);
if (integer <= 0)
    integer = 1;
elseif (integer >= upper_bound)
    integer = upper_bound - 1;
end

end
