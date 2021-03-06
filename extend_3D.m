% <User Specified Variables>

FILE_Matlab_output = '/Users/yuquan/GitHub/de_terminal_copy/Matlab_output.dat';

scale_factor = 0.027;   % computed from pixel to mm conversion of oct image

num_replicates = 40;

true_data_frame = [5,7,9,11,13,15,17,19,21,25];

Lens_Node_Density = 60; %number of nodes on the half lens curve

Order_Lens_Poly = 3;

% </User Specified Variables>


origin = 0.5*(meshwork_l{1} + meshwork_r{1});   % mid-point of the image


% <UI to mark the lens contour manually>
% Look at the text instruction of the UI.
figure, imshow( 255 - uint8(original{1}) ); hold on
plot( origin(2), origin(1), 'or', 'linew', 7 );
h_text = uicontrol('style', 'text', 'string', 'mark the center', ...
     'unit', 'normalized', 'position', [0.1 0.1 0.2 0.03]);
 
 [center(1),center(2)] = ginput(1);
 plot( center(1), center(2), 'or');
 
 index = 1;
 scatter(iris_contour_l{index}(:,2), iris_contour_l{index}(:,1), '.g');
 scatter(iris_contour_r{index}(:,2), iris_contour_r{index}(:,1), '.r');
 scatter(meshwork_l{index}(2), meshwork_l{index}(1), 'ro', 'linew', 8);
 scatter(meshwork_r{index}(2), meshwork_r{index}(1), 'ro', 'linew', 8);
 scatter(intercept_l{index}(2), intercept_l{index}(1), 'or', 'linew', 8);
 scatter(intercept_r{index}(2), intercept_r{index}(1), 'or', 'linew', 8);
 scatter(extended_cornea_dense{index}(:,2), extended_cornea_dense{index}(:,1), '.m');
 scatter(tip_l{index}(:,2), tip_l{index}(:,1), 'ob', 'linew', 8);
 scatter(tip_r{index}(:,2), tip_r{index}(:,1), 'ob', 'linew', 8);
 
 set(h_text, 'string', 'mark points on the left lens bounary, then press Enter');
 [lens_leftx_marked, lens_lefty_marked] = ginput;
 scatter( lens_leftx_marked, lens_lefty_marked, 'bx');
 
  set(h_text, 'string', 'mark points on the right lens bounary, then press Enter');
 [lens_rightx_marked, lens_righty_marked] = ginput;
 scatter( lens_rightx_marked, lens_righty_marked, 'bx');
 % </UI>
 
 
 % <polynomial fit of the lens contour>
 lens_leftx_marked = [center(1); lens_leftx_marked]; lens_lefty_marked = [center(2); lens_lefty_marked];
 lens_rightx_marked = [center(1); lens_rightx_marked]; lens_righty_marked = [center(2); lens_righty_marked];
 
 lens_pleft = polyfit(lens_leftx_marked, lens_lefty_marked, Order_Lens_Poly);    % Order 3 polynomial
 lens_pright = polyfit(lens_rightx_marked, lens_righty_marked, Order_Lens_Poly); % Order 3 polynomial, to fit the half lens.
 % </polynomial fit of the lens contour>
 
 % <Plot the lens>
 leftxplot = min(lens_leftx_marked) - 100:0.1:max(lens_leftx_marked);
 leftyplot = polyval(lens_pleft, leftxplot);
 rightxplot = min(lens_rightx_marked):0.1:max(lens_rightx_marked)+100;
 rightyplot = polyval(lens_pright, rightxplot);
 plot(leftxplot, leftyplot, 'r');
 plot(rightxplot, rightyplot, 'r');
 % </Plot the lens>
 
 
 % <Left and right ends of the lens curve>
 lens_node_minx = min(lens_leftx_marked) - 30;
 lens_node_maxx = max(lens_rightx_marked) + 30;
 % </Left and right ends of the lens curve>

 
 % the nodes are ordered from the center to the sides.
 lens_nodes_rightx = linspace(center(1)+5, lens_node_maxx, Lens_Node_Density);
 lens_nodes_leftx = flip(linspace(lens_node_minx, center(1)-5, Lens_Node_Density));  % these are row vectors
 
 lens_nodes_righty = polyval(lens_pright, lens_nodes_rightx);  %also row vectors
 lens_nodes_lefty = polyval(lens_pleft, lens_nodes_leftx);
 
 % node coordinates
 lens_nodes_left = [lens_nodes_leftx; lens_nodes_lefty]';
 lens_nodes_right = [lens_nodes_rightx; lens_nodes_righty]';
 
 % <Element node connectivity>
 lens_node_index = [1:Lens_Node_Density-1; 2:Lens_Node_Density]';
 % </Element node connectivity>
 
 % <create 1D lens elements (line segments)>
 for index = 1:length(lens_node_index)
     lens_elements(index) = element( element_type.line2, lens_node_index(index,:), 10 );
 end
  % </create 1D lens elements (line segments)>
  
 
 % <creat 1D lens mesh>
 lens_mesh_obj_left = mesh_class(lens_nodes_left, lens_elements);
 lens_mesh_obj_right = mesh_class(lens_nodes_right, lens_elements);
 % </create 1D lens mesh>
 
 % <3D iris mesh, axial symmetric>
mesh_opt = struct('mesh_size', [5,5], 'imsiz', size(original{1}), 'dilator_muscle',1, ...
    'dilator_width', 1, 'do_refine',false, 'sphincter_muscle', 1, 'dilator_retreat', 1, 'sphincter_wholewidth', 0);

% mesh object: iris object
[mesh_obj, nx, ny, mesh_right, boundary_node_set, revolve_obj] = mesh_iris_3D( top_iris_r{1}, bottom_iris_r_smooth{1}, origin(2), num_replicates,3,2,mesh_opt);
% </3D iris mesh, axial symmetric>

% <mesh the left iris contour, 2D>
mesh_opt.left_right = 'left';
mesh_opt.mesh_size = [nx, ny];
mesh_left_original = mesh_iris_2D( top_iris_l{1},bottom_iris_l_smooth{1}, mesh_opt);
% </mesh the left iris contour, 2D>

% <The non-axial-symmetric 3D mesh>
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
% </Then non-axial-symmetric 3D mesh>

% <The non-axial-symmetric 3D lens mesh>
% revolve the lens mesh, into 3D surface
lens_revolve_obj = revolve_mesh(lens_mesh_obj_right, revolve_obj.plane, revolve_obj.axis, revolve_obj.origin, revolve_obj.num_replicates);
lens_mesh_obj_right_3D = generate3D(lens_revolve_obj);

mesh_left_right_lens = lens_mesh_obj_left;
mesh_left_right_lens.node_list(:,1) = - mesh_left_right_lens.node_list(:,1) + 2* origin(2);


num_node_2D_lens = length( mesh_left_right_lens.node_list );
right_node_lens = lens_mesh_obj_right_3D.node_list( 1: num_node_2D_lens,:);

left_right_diff_lens = mesh_left_right_lens.node_list - lens_mesh_obj_right.node_list;

for i = 1:num_replicates/2
    angle = i*(2*pi/num_replicates);
    % modify the ith section, which is mesh_obj.node_list( 1 + (i)* num_node_2D : (i+1)*num_node_2D,:)
    
    modification_2D = (angle / pi) * left_right_diff_lens;
    modification_projected = [modification_2D(:,1) * cos(angle), modification_2D(:,2), modification_2D(:,1)* sin(angle)];    
    lens_mesh_obj_right_3D.node_list( 1 + (i)* num_node_2D_lens : (i+1)*num_node_2D_lens,:) = ...
        lens_mesh_obj_right_3D.node_list( 1 + (i)* num_node_2D_lens : (i+1)*num_node_2D_lens,:) + modification_projected;
    
   if i ~= num_replicates/2
        modification_projected = [modification_2D(:,1) * cos(angle), modification_2D(:,2), - modification_2D(:,1)* sin(angle)];
        lens_mesh_obj_right_3D.node_list( num_node_2D_lens * num_replicates - i*num_node_2D_lens + 1: num_node_2D_lens * num_replicates - (i-1)*num_node_2D_lens ,:) = ...
            lens_mesh_obj_right_3D.node_list( num_node_2D_lens * num_replicates - i*num_node_2D_lens + 1: num_node_2D_lens * num_replicates - (i-1)*num_node_2D_lens ,:) +...
            modification_projected;
   end
end

mesh_obj = mesh_obj + lens_mesh_obj_right_3D;
% </The non-axial-symmetric 3D lens mesh>


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