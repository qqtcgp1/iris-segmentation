function [ iris_left_bottom, iris_right_bottom, boundary ] = iris_bottom_fit( boundary, critical_points_A )
%iris_bottom_fit finds the bestfit iris bottom boundary as output in
%iris_left_bottom and iris_right_bottom, in linear indices form in a column.
%critical_points_A = [a1; a2; a3; a4; a5; a6] found as in critical_points.


a3 = critical_points_A(3,:);
a6 = critical_points_A(6,:);
a2 = critical_points_A(2,:);
a5 = critical_points_A(5,:);


path_ = find_path(boundary, a6, a3);
path = zeros(length(path_), 2);
 
for i = 1:length(path_)
     [path(i,1), path(i,2)] = ind2sub(size(boundary), path_(i));
end

path_adjusted = path;
for i = 1:length(path_)
    path_adjusted(i,:) = path(i,:) - a3;
end


linear_coeff = path_adjusted(:,2)\path_adjusted(:,1);

iris_right_bottom_ = zeros(size(boundary,2) - a3(2), 2);

for i = 1:size(boundary,2)-a3(2)
    iris_right_bottom_(i, :) = [nearest_(i*linear_coeff, size(boundary,1)) + a3(1), nearest_(i + a3(2), size(boundary,2))];
end


if nargout == 3
    boundary(path_) = 0;
    for i = 1:size(iris_right_bottom_,1)
        boundary(iris_right_bottom_(i,1), ...
            iris_right_bottom_(i,2)) = 1;
    end
end



path_ = find_path(boundary, a5, a2);
path = zeros(length(path_), 2);
 
for i = 1:length(path_)
     [path(i,1), path(i,2)] = ind2sub(size(boundary), path_(i));
end


% do the best-fit
path_adjusted = path;
for i = 1:length(path_)
    path_adjusted(i,:) = path(i,:) - a2;
end

linear_coeff = path_adjusted(:,2)\path_adjusted(:,1);

iris_left_bottom_ = zeros(a2(2), 2);

for i = 1:a2(2)
    iris_left_bottom_(i, :) = [nearest_((-i)*linear_coeff + a2(1), size(boundary,1)), nearest_(-i + a2(2), size(boundary,2))];
end


if nargout == 3
    boundary(path_) = 0;
    for i = 1:size(iris_left_bottom_,1)
        boundary(iris_left_bottom_(i,1), ...
            iris_left_bottom_(i,2)) = 1;
    end
    
    boundary(a2) = 1; boundary(a3) = 1;
end

iris_left_bottom = zeros(size(iris_left_bottom_,1), 1);
for i = 1:size(iris_left_bottom_,1)
    iris_left_bottom(i) = sub2ind(size(boundary), iris_left_bottom_(i,1), iris_left_bottom_(i,2));
end
iris_right_bottom = zeros(size(iris_right_bottom_,1), 1);
for i = 1:size(iris_right_bottom_,1)
    iris_right_bottom(i) = sub2ind(size(boundary), iris_right_bottom_(i,1), iris_right_bottom_(i,2));
end
    

 

end