FILE_Matlab_output = '/Users/yuquan/GitHub/de_terminal_copy/Matlab_output.dat';

scale_factor = 0.027;   % computed from pixel to mm conversion of oct image

num_replicates = 40;

true_data_frame = [5,7,9,11,13,15,17,19,21,25];


origin = 0.5*(meshwork_l{1} + meshwork_r{1});   % mid-point of the image
mesh_opt = struct('mesh_size', [5,5], 'imsiz', size(original{1}), 'dilator_muscle',1, ...
    'dilator_width', 1, 'do_refine',false, 'sphincter_muscle', 1, 'dilator_retreat', 1, 'sphincter_wholewidth', 0);

% mesh object
[mesh_obj, nx, ny, mesh_right, boundary_node_set] = mesh_iris_3D( top_iris_r{1}, bottom_iris_r_smooth{1}, origin(2), num_replicates,3,2,mesh_opt);

mesh_opt.left_right = 'left';
mesh_opt.mesh_size = [nx, ny];
mesh_left_original = mesh_iris_2D( top_iris_l{1},bottom_iris_l_smooth{1}, mesh_opt);
mesh_left_right = mesh_left_original;

mesh_left_right.node_list(:,1) = - mesh_left_original.node_list(:,1) + 2* origin(2);

num_node_2D = length( mesh_left_original.node_list );
right_node = mesh_obj.node_list( 1: num_node_2D,:);


left_right_diff = mesh_left_right.node_list - mesh_right.node_list;

for i = 1:num_replicates/2
    angle = i*(2*pi/num_replicates);
    % modify the ith section, which is mesh_obj.node_list( 1 + (i)* num_node_2D : (i+1)*num_node_2D,:)
    
    modification_2D = (angle / pi) * left_right_diff;
    modification_projected = [modification_2D(:,1) * cos(angle), modification_2D(:,2), modification_2D(:,1)* sin(angle)];    
    mesh_obj.node_list( 1 + (i)* num_node_2D : (i+1)*num_node_2D,:) = ...
        mesh_obj.node_list( 1 + (i)* num_node_2D : (i+1)*num_node_2D,:) + modification_projected;
    
   if i ~= num_replicates/2
        modification_projected = [modification_2D(:,1) * cos(angle), modification_2D(:,2), - modification_2D(:,1)* sin(angle)];
        mesh_obj.node_list( num_node_2D * num_replicates - i*num_node_2D + 1: num_node_2D * num_replicates - (i-1)*num_node_2D ,:) = ...
            mesh_obj.node_list( num_node_2D * num_replicates - i*num_node_2D + 1: num_node_2D * num_replicates - (i-1)*num_node_2D ,:) +...
            modification_projected;
   end
end

mesh_obj = flip_orientation(mesh_obj);


% one-liner transformation, combined in a single equation:
zmin = min(-mesh_obj.node_list(:,2)); l = length(mesh_obj.node_list);
mesh_obj.node_list = scale_factor * ( mesh_obj.node_list * [1,0,0; 0,0,-1; 0,1,0] - [ones(l,1)*origin(2), zeros(l,1), zmin* ones(l,1)]);

% say the transformation is Y = X * M + C, then M and C is as this:
M = scale_factor * [1,0,0; 0,0,-1; 0,1,0];
C = -scale_factor * [ones(l,1)*origin(2), zeros(l,1), zmin* ones(l,1)];


% % Transformations, in steps.
% 
% % modify the coordinates to fit the cylindrical axis
% % reverse y. exchange y and z. so that the three axis are still in positive
% % orientation. (meaning ex x ey = ez, not -ez).
% 
% 
% mesh_obj.node_list(:,2) = - mesh_obj.node_list(:,2);
% mesh_obj.node_list(:,[2 3]) = mesh_obj.node_list(:, [3 2]);
% 
% % translate in the x-y plane, such that the centre coincide with the
% % origin. only need to subtract first column by origin(2)
% 
% mesh_obj.node_list(:,1) = mesh_obj.node_list(:,1) - origin(2);
% 
% % translate in the z direction, so that the minimum z value is 0.
% zmin = min( mesh_obj.node_list(:,3) );
% mesh_obj.node_list(:,3) = mesh_obj.node_list(:,3) - zmin;
% mesh_obj.node_list = mesh_obj.node_list * scale_factor;   %pixel to mm conversion.

fileID = fopen(FILE_Matlab_output, 'w');
fprintf(fileID, 'Number of boundary nodes: %d\n', length(boundary_node_set));
fprintf(fileID, 'Number of nodes: %d\n', mesh_obj.num_nodes);
fprintf(fileID, 'Number of elements: %d\n<boundary node set>\n', mesh_obj.num_elements);
fprintf(fileID, '%6d  ',boundary_node_set);
fprintf(fileID, '\n</boundary node set>\n\n');
% output to txt file (feb format)
print_feb(mesh_obj, fileID);
fclose(fileID);


for i = 1:length(iris_contour_r)
    iris_contour_r{i} = iris_contour_r{i} * scale_factor;
    iris_contour_l{i} = iris_contour_l{i} * scale_factor;
    
    iris_contour_r{i} = iris_contour_r{i}(:,[2 1]);
    iris_contour_l{i} = iris_contour_l{i}(:,[2 1]);
    
    iris_contour_r{i}(:,2) = -iris_contour_r{i}(:,2);
    iris_contour_l{i}(:,2) = -iris_contour_l{i}(:,2);
    
    iris_contour_r{i}(:,2) = iris_contour_r{i}(:,2) - min(iris_contour_r{i}(:,2));
    iris_contour_l{i}(:,2) = iris_contour_l{i}(:,2) - min(iris_contour_l{i}(:,2));
end

FILE_Matlab_output = '/Users/yuquan/GitHub/de_terminal_copy/true_data.dat';
fileID = fopen(FILE_Matlab_output, 'w');
for i = 1:length(true_data_frame)
    % left iris
    fprintf(fileID, 'step %d left, %d points\n', i, length(iris_contour_l{true_data_frame(i)}));
    for j = 1:length(iris_contour_l{true_data_frame(i)})
        fprintf(fileID, '%10g %10g\n', iris_contour_l{true_data_frame(i)}(j,:));
    end
    % right iris
    fprintf(fileID, 'step %d right, %d points\n', i, length(iris_contour_r{true_data_frame(i)}));
    for j = 1:length(iris_contour_r{true_data_frame(i)})  %%%% fixed
        fprintf(fileID, '%10g %10g\n', iris_contour_r{true_data_frame(i)}(j,:));
    end
end