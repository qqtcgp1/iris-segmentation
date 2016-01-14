function [ volume ] = iris_volume( iris_solid, centre )
%UNTITLED Summary of this function goes here
%   centre is the horizontal position of the mid-point of the iris tips.


S = size(iris_solid);
v = zeros(S(2),1);

for i = 1:S(2)
    v(i) = pi*length(find(iris_solid(:,i)))*abs(centre-i);
end

volume = sum(v);


end