function [ mesh_obj_3D, nx, ny, mesh_obj_2D,boundary_node_set, revolve_obj ] = mesh_iris_3D( top_iris, bottom_iris, origin, num_replicates, plane, axis, ...
    mesh_opt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin <= 5
    plane = 3; axis = 2;
end

origin = origin * ones(1,3);
origin(plane) = 0; origin(axis) = 300;

[mesh_obj_2D, nx, ny,~,boundary_node_set] = mesh_iris_2D(top_iris, bottom_iris, mesh_opt);

revolve_obj = revolve_mesh( mesh_obj_2D, plane, axis, origin, ...
    num_replicates);

mesh_obj_3D = generate3D(revolve_obj);

end

