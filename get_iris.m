function [ iris_pts, area, volume, A ] = get_iris( boundary )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

try
    [ mesh_left, tip_left, tip_right, mesh_right, mesh_left_bottom, ...
        mesh_right_bottom ] = critical_points( boundary );
catch
    throwAsCaller(MException('yyy:ddd', 'cannot do critical_points'));
end

    
A = [mesh_left; tip_left; tip_right; mesh_right; mesh_left_bottom; mesh_right_bottom];

left_top = find_path(boundary, mesh_left, tip_left);
left_top = ind2sub_array(size(boundary), left_top);

left_left = connect_line( mesh_left, mesh_left_bottom);
right_right = connect_line(mesh_right, mesh_right_bottom);
right_top = find_path(boundary, mesh_right, tip_right);
right_top = ind2sub_array(size(boundary), right_top);

[left_bottom, right_bottom] = iris_bottom_fit( boundary, A);
left_bottom = ind2sub_array(size(boundary), left_bottom);
right_bottom = ind2sub_array(size(boundary), right_bottom);

iris_pts = [left_top; left_left; right_right; right_top; left_bottom; right_bottom];

if nargout >= 2
    solid = imfill( points2bw(iris_pts, size(boundary)), 'holes' );
    area = sum(sum(solid));
end

if nargout >= 3
    volume = iris_volume(solid, 0.5*(A(2,2)+A(3,2)));
end

end

