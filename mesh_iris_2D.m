function [ mesh_obj_2D, n_int_x, n_int_y, bottom_nodes, boundary_node_set ] = mesh_iris_2D( top_iris, bottom_iris, mesh_opt)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% figure, hold on; axis equal; set(gca, 'ydir', 'reverse');axis([0 600 0 300]);
% for i = 1:length(top_iris)
%     scatter(top_iris(i,2), top_iris(i,1),'.');
%     pause(0.1);
% end
% for i = 1:length(bottom_iris)
%     scatter(bottom_iris(length(bottom_iris)+1-i,2), bottom_iris(length(bottom_iris)+1-i,1), '.');
%     pause(0.1);
% end
% polyarea([top_iris(:,1);flipud(bottom_iris(:,1))], [top_iris(:,2);flipud(bottom_iris(:,2))])
% above is for debugging

boundary_node_top = []; boundary_node_bot = [];
mesh_size = mesh_opt.mesh_size;
imsiz = mesh_opt.imsiz;
if isfield(mesh_opt, 'left_right')
    left_right = mesh_opt.left_right;
else left_right = 'right';
end

if isfield(mesh_opt, 'dilator_muscle')
    dilator_muscle = mesh_opt.dilator_muscle;
    dilator_width = mesh_opt.dilator_width;
    if isfield(mesh_opt, 'do_refine');
        do_refine = mesh_opt.do_refine;
    else do_refine = 0;
    end
else dilator_muscle = 0;
end


if strcmp( left_right, 'left')
    [mesh_obj_2D, n_int_x, n_int_y, bottom_nodes, boundary_node_set] = mesh_iris_2D( fliplr_points( top_iris, imsiz), ...
        fliplr_points( bottom_iris, imsiz), setfield(mesh_opt, 'left_right', 'right'));
    mesh_obj_2D.node_list(:,1) = - mesh_obj_2D.node_list(:,1) + imsiz(2) + 1;
    return
end

if isfield(mesh_opt, 'sphincter_muscle')
    sphincter_muscle = mesh_opt.sphincter_muscle;
else
    sphincter_muscle = 0;
end

if sphincter_muscle && isfield(mesh_opt, 'sphincter_wholewidth')
    sphincter_wholewidth = mesh_opt.sphincter_wholewidth;
else
    sphincter_wholewidth = 0;
end

if isfield(mesh_opt, 'dilator_retreat')
    dilator_retreat = mesh_opt.dilator_retreat;
else
    dilator_retreat = false;
end


    function result = foo(node_point)
        [~, result] = ismembertol(node_point, mesh_obj_2D.node_list, 0.0000001, 'ByRows', 1);
    end

if dilator_muscle
    % sort pts in horizontal order
    [~,ind] = sort( bottom_iris(:,2));
    bottom_iris = bottom_iris(ind,:);
    [~,ind] = sort( top_iris(:,2));
    top_iris = top_iris(ind,:);
    
    if abs(min(bottom_iris(:,2)) - min(top_iris(:,2)))  > 0.1
        bottom_iris = [top_iris(1,:); bottom_iris];
    end
    
    % bottom_iris_modified is the inner boundary
    bottom_iris_modified = bottom_iris - [[0;dilator_width*ones(length(bottom_iris)-1,1)],zeros(length(bottom_iris),1)];
    
    [mesh_obj_2D_, n_int_x, n_int_y, bottom_nodes, boundary_node_set] = mesh_iris_2D( top_iris, bottom_iris_modified, setfield(mesh_opt, 'dilator_muscle', 0));
    dilator_top = [bottom_iris(1,[2 1]); fliplr(bottom_nodes(2:end,:))];
    dilator_bot = fliplr(bottom_nodes(2:end,:)) + [zeros(length(bottom_nodes)-1,1), dilator_width*ones(length(bottom_nodes)-1,1)];
    
    nodes = [dilator_top; dilator_bot];
    dilator_bot = [bottom_iris(1,[2,1]); dilator_bot];
    
    e(1) = element(element_type.tri3, [1, 2, 2+length(bottom_nodes)-1],2);
    for i = 1:length(bottom_nodes)-1 - 1
        e(i+1) = element(element_type.quad4, [1+i, 2+i, 2+i+length(bottom_nodes)-1, 1+i+length(bottom_nodes)-1],2);
    end
    
    if dilator_retreat
        if sphincter_muscle
        e(1).materialID = 4;
        for w = 2:5
            e(w).materialID = 3;
        end
        else
            for w = 1:5
                e(w).materialID = 1;
            end
        end           
    end
        
    mesh_obj_2D = mesh_obj_2D_ + mesh_class(nodes, e); % this concludes the rough mesh of dilator.
    
    
    if do_refine
        n_refine = 3; % this number should be close to sqrt(mesh_size/dilator_width), in order to harvest the best aspect ratio of elements
        
        fine_node_top = cell(length(dilator_bot)-1,1);
        fine_node_bot = fine_node_top;
        for i = 1:length(dilator_bot)-1
            % computes the nodes of the finer mesh of the dilator layer
            fine_node_top{i} = zeros(n_refine-1,2);
            fine_node_bot{i} = zeros(n_refine-1,2);
            for j = 1:n_refine-1
                fine_node_top{i}(j,:) = dilator_top(i,:) + (j/n_refine)*(dilator_top(i+1,:) - dilator_top(i,:));
                fine_node_bot{i}(j,:) = dilator_bot(i,:) + (j/n_refine)*(dilator_bot(i+1,:) - dilator_bot(i,:));
            end
        end
        
        % add these nodes to the node list
        old_n_nodes = length(mesh_obj_2D.node_list);
        mesh_obj_2D.node_list = [mesh_obj_2D.node_list; cell2mat(fine_node_top)];
        mesh_obj_2D.node_list = [mesh_obj_2D.node_list; cell2mat(fine_node_bot)]; clear temp
        
        % find the upper nodes of the 2nd bottom most layer of the main
        % iris. The nodes should be iris_bot2(:,:), and node number is
        % ind(:)
        iris_bot2 = fliplr(bottom_nodes(2:end,:));
        ind2 = zeros(1, length(iris_bot2));
        for j = 1:length(iris_bot2)
            temp = abs(mesh_obj_2D_.node_list(:,1)-iris_bot2(j,1)) < 1e-6;
            [mx,indx] = max(mesh_obj_2D_.node_list(temp,2));
            temp = find(temp, indx+1);
            indx = temp(indx);
            mesh_obj_2D_.node_list(indx,2) = -Inf;
            m = max(mesh_obj_2D_.node_list(temp,2));
            mesh_obj_2D_.node_list(indx,2) = mx;
            iris_bot2(j,2) = m;
        end
        
        % eliminate those divided elements
        mesh_obj_2D.elements = setdiff(mesh_obj_2D.elements, element(element_type.tri3, [foo(dilator_top(1,:)), foo(dilator_top(2,:)), foo(dilator_bot(2,:))]));
        mesh_obj_2D.elements = setdiff(mesh_obj_2D.elements, element(element_type.tri3, [foo(dilator_top(1,:)), foo(dilator_top(2,:)), foo(iris_bot2(1,:))]));
        
        for i = 2:length(dilator_bot)-1
            mesh_obj_2D.elements = setdiff(mesh_obj_2D.elements, element(element_type.quad4, ...
                [foo(dilator_top(i,:)), foo(dilator_top(i+1,:)), ...
                foo(dilator_bot(i+1,:)), foo(dilator_bot(i,:))]));
            
            mesh_obj_2D.elements = setdiff(mesh_obj_2D.elements, element(element_type.quad4, ...
                [foo(dilator_top(i,:)), foo(dilator_top(i+1,:)), foo(iris_bot2(i,:)), foo(iris_bot2(i-1,:))]));
        end
        
        % build the new refined elements!
        % first, the muscle strip. should be easy.
        % use the function foo
        
        e = element(element_type.tri3, [0 0 0]);
        e(1) = [];
        
        for i = 1:length(dilator_bot)-1
            if i == 1
                e = [e element(element_type.tri3, [foo(dilator_top(1,:)), ...
                    foo(fine_node_top{1}(1,:)), foo(fine_node_bot{1}(1,:))],2)];
            else
                e = [e element(element_type.quad4, [foo(dilator_top(i,:)), ...
                    foo(fine_node_top{i}(1,:)), foo(fine_node_bot{i}(1,:)), ...
                    foo(dilator_bot(i,:))], 2)];
            end
            
            for j = 1:n_refine-2
                e = [e element(element_type.quad4, [foo(fine_node_top{i}(j,:)), ...
                    foo(fine_node_top{i}(j+1,:)), foo(fine_node_bot{i}(j+1,:)), ...
                    foo(fine_node_bot{i}(j,:))],2)];
            end
            
            e = [e element(element_type.quad4, [foo(fine_node_top{i}(n_refine-1,:)), ...
                foo(dilator_top(i+1,:)), foo(dilator_bot(i+1,:)), foo(fine_node_bot{i}(n_refine-1,:))], 2)];
        end
        
        e = [e element(element_type.tri3, [1, foo(iris_bot2(1,:)), ...
            foo(fine_node_top{1}(1,:))], 1)];
        
        for i = 1:n_refine-2
            e = [e element(element_type.tri3, [foo(fine_node_top{1}(i,:)),...
                foo(iris_bot2(1,:)), foo(fine_node_top{1}(i+1,:))],1)];
        end
        
        e = [e element(element_type.tri3, [foo(fine_node_top{1}(n_refine-1,:)), ...
            foo(iris_bot2(1,:)), foo(dilator_top(2,:))],1)];
        
        for i = 2:length(dilator_bot)-1
            e = [e element(element_type.tri3, [foo(dilator_top(i,:)), ...
                foo(iris_bot2(i-1,:)), foo(fine_node_top{i}(1,:))],1)];
            e = [e element(element_type.tri3, [foo(iris_bot2(i,:)), ...
                foo(dilator_top(i+1,:)),foo(fine_node_top{i}(end,:))],1)];
            n_left = floor((n_refine - 2)/2.0);
            n_right = n_refine - 2 - n_left;
            for j = 1:n_left
                e = [e element(element_type.tri3, [foo(iris_bot2(i-1,:)),...
                    foo(fine_node_top{i}(j+1,:)), foo(fine_node_top{i}(j,:))],1)];
            end
            e = [e element(element_type.tri3, [foo(iris_bot2(i-1,:)), ...
                foo(iris_bot2(i,:)), foo(fine_node_top{i}(n_left+1,:))],1)];
            for j = 1:n_right
                e = [e element(element_type.tri3, [foo(fine_node_top{i}(j+n_left,:)),...
                    foo(iris_bot2(i,:)), foo(fine_node_top{i}(j+n_left+1,:))],1)];
            end
        end
    end
    
    %mesh_obj_2D.elements = [mesh_obj_2D.elements e];
    return;
end




[~, tip_index_top] = min(top_iris(:,2));
iris_tip = top_iris(tip_index_top,:);
[~, tip_index_bot] = min(top_iris(:,2));


nominal_int_x = 0; nominal_int_y = 0;
n_int_x = 0; n_int_y = 0;

if numel(mesh_size) == 1
    mesh_size = [mesh_size mesh_size];
end

distance_x = max(top_iris(:,2)) - min(top_iris(:,2));
int_input = (numel(mesh_size) == 2);

if int_input
    nominal_int_x = mesh_size(1);
    nominal_int_y = mesh_size(2);
    n_int_x = nearest(distance_x / nominal_int_x);
    n_int_y = zeros(1, n_int_x + 1);
else
    n_int_x = mesh_size(1);
    n_int_y = mesh_size(2:end);
    nominal_int_x = nearest(distance_x / n_int_x);
end


n_node_x = n_int_x + 1;

% initialize top nodes, from x positions.
h = 1;

top_x = zeros(1, n_node_x);
mmma = max( top_iris(:,2));
mmmi = nearest(min( top_iris(:,2)));
if abs(min(bottom_iris(:,2)) - mmmi)  > 0.1
    bottom_iris = [top_iris(1,:); bottom_iris];
end


for i =0 : n_int_x
    top_x(h) = mmmi + nearest(i*distance_x / n_int_x);
    h = h+1;
end


for h = 1: 1+n_int_x
    if top_x(h) > mmma
        top_x(h) = mmma;
    end
    temp = find( abs(top_iris (:,2) - top_x(h))<0.5,1, 'first');
    
    % if can't find the coordinate, the problem is usually due to the last
    % node of top_iris. Thus a run without this node solves the problem.
    if isempty(temp)
        [ mesh_obj_2D, n_int_x, n_int_y, bottom_nodes ] = mesh_iris_2D( top_iris(1:end-1,:), bottom_iris, mesh_opt);
        return;
    end
    top_y(h) = top_iris ( temp, 1);
end

% calculates y-witdh at different x, and hence n_int_y at different x
distance_y = zeros(length(top_x),1);
bottom_nodes = zeros( length(top_x),2);

for i = 1: length(distance_y)
    j =  find(abs(bottom_iris (:,2) - top_x(i)) < 0.5, 1, 'first') ;
    distance_y(i) = bottom_iris(j,1)- top_y(i);
    bottom_nodes(i,:) = bottom_iris(j,:);
end



if int_input
    for i = 1:length(distance_y)
        if distance_y(i) >= nominal_int_y
            n_int_y(i) = nearest(distance_y(i)/nominal_int_y);
        end
    end
else
    n_int_y = mesh_size(2:end);
end

% tip
gcoord_x(1) = iris_tip(2);
gcoord_y(1) = iris_tip(1);
node_counter = 2;

% all coordinates
for i = 2: length(n_int_y)
    j =  find(abs(bottom_iris(:,2) - top_x(i)) < 0.5, 1, 'first');
    
    if int_input
        temp = nominal_int_y;
    else
        temp = nearest( distance_y(i) / n_int_y(i));
    end
    
    
    % temp is the interval
    if (bottom_iris(j,1) - top_y(i)) >= temp || (~int_input)
        for k = 0: n_int_y(i)
            gcoord_x(node_counter) = top_x(i);
            gcoord_y(node_counter) = top_y(i)+ nearest( k * distance_y(i) / n_int_y(i) );
            % unless this node is v close to the boundary!
            if k==0
                boundary_node_top = [boundary_node_top,node_counter];
            elseif k==n_int_y(i)
                gcoord_y(node_counter) = bottom_iris(j,1); %#ok<*AGROW>
                boundary_node_bot = [boundary_node_bot,node_counter];
            end
            
            node_counter = node_counter+1;
        end
    end
end

gcoord = [gcoord_x;gcoord_y];
gcoord = gcoord';

n_nodes_y = n_int_y + 1.0;
sphincterTriSet = [];
sphincterQuadSet = [];


ntri = 1;
nhex = 1;
node_index = 1;
% triangle for the first element at the tip
node_tri (ntri, 1) = node_index;
node_tri (ntri, 2) = node_index+1;
node_tri (ntri, 3) = node_index+2;
ntri = ntri+1;

if n_int_y(2) == 2
    node_tri(ntri,:) = [1,3,4];
    sphincterTriSet = [sphincterTriSet ntri];
    ntri = ntri+1;
else
    sphincterTriSet = [sphincterTriSet ntri-1];
end

node_index = node_index+1;
tri_cuttoff = 0; quad_cuttoff = 0;

for i = 2: length (n_int_y) - 1
    if sphincter_wholewidth && i==6
        tri_cuttoff = ntri - 1;
        quad_cuttoff = nhex - 1;
    end
    if n_nodes_y (i+1) - n_nodes_y(i) == 2
        node_tri (ntri, 1) = node_index;
        node_tri (ntri,2) = node_index + n_nodes_y(i);
        node_tri (ntri,3) = node_tri(ntri,2) + 1;
        ntri = ntri+1;
        for j = 1 : n_nodes_y (i)-1
            node_hex (nhex, 1) = node_index;
            node_hex (nhex,2) = node_index + n_nodes_y(i)+1;
            node_hex(nhex,3) = node_hex(nhex,2)+1;
            node_hex(nhex,4) = node_hex(nhex,1) +1;
            nhex = nhex+1;
            node_index = node_index+1;
        end
        node_tri(ntri,1) = node_index;
        node_tri(ntri,2) = node_index + n_nodes_y (i+1) -1;
        node_tri(ntri,3) = node_tri(ntri,2) +1;
        if (i <= 5)
            sphincterTriSet = [sphincterTriSet, ntri];
        end
        ntri =ntri+1;
        node_index = node_index +1;
    end
    
    if n_nodes_y (i) - n_nodes_y(i+1) == 2
        node_tri (ntri, 1) = node_index;
        node_tri (ntri,2) = node_index + n_nodes_y(i);
        node_tri (ntri,3) = node_tri(ntri,1) + 1;
        ntri = ntri+1;
        node_index = node_index +1;
        for j = 1 : n_nodes_y (i+1)-1
            node_hex (nhex, 1) = node_index;
            node_hex (nhex,2) = node_index + n_nodes_y(i) -1;
            node_hex(nhex,3) = node_hex(nhex,2)+1;
            node_hex(nhex,4) = node_hex(nhex,1) +1;
            nhex = nhex+1;
            node_index = node_index+1;
        end
        node_tri(ntri,1) = node_index;
        node_tri(ntri,2) = node_index + n_nodes_y (i+1) + 1;
        node_tri(ntri,3) = node_tri(ntri,1) + 1;
        if (i <= 5 )
            sphincterTriSet = [sphincterTriSet, ntri];
        end
        
        ntri =ntri+1;
        node_index = node_index + 2;
    end
    
    if n_nodes_y(i+1) - n_nodes_y (i) ==1
        node_tri (ntri,1) = node_index;
        node_tri(ntri,2) = node_index + n_nodes_y(i);
        node_tri(ntri,3) = node_tri(ntri,2) +1;
        ntri= ntri+1;
        for j = 1 : n_nodes_y(i)-1
            node_hex (nhex, 1) = node_index;
            node_hex (nhex,2) = node_index + n_nodes_y(i)+1;
            node_hex(nhex,3) = node_hex(nhex,2)+1;
            node_hex(nhex,4) = node_hex(nhex,1) +1;
            if (i <= 5) && (j == n_nodes_y(i)-1)
                sphincterQuadSet = [sphincterQuadSet, nhex];
            end
            nhex = nhex+1;
            node_index = node_index+ 1;
        end
        node_index = node_index+ 1;
    end
    
    if n_nodes_y(i) - n_nodes_y (i+1) == 1
        node_tri (ntri,1) = node_index;
        node_tri(ntri,2) = node_index + n_nodes_y(i);
        node_tri(ntri,3) = node_tri(ntri,1) +1;
        ntri= ntri+1;
        node_index = node_index+1;
        
        for j = 1 : n_nodes_y (i+1)-1
            node_hex (nhex, 1) = node_index;
            node_hex (nhex,2) = node_index + n_nodes_y(i)-1;
            node_hex(nhex,3) = node_hex(nhex,2)+1;
            node_hex(nhex,4) = node_hex(nhex,1) +1;
            if (i <= 5) && (j == n_nodes_y (i+1)-1)
                sphincterQuadSet = [sphincterQuadSet, nhex];
            end
            nhex = nhex+1;
            node_index = node_index+1;
        end
        node_index = node_index+1;
    end
    
    if n_nodes_y (i) == n_nodes_y (i+1)
        for j = 1 : n_nodes_y (i)-1
            node_hex (nhex, 1) = node_index;
            node_hex (nhex,2) = node_index + n_nodes_y(i);
            node_hex(nhex,3) = node_hex(nhex,2)+1;
            node_hex(nhex,4) = node_hex(nhex,1) +1;
            if i <= 5 && j == n_nodes_y (i)-1
                sphincterQuadSet = [sphincterQuadSet, nhex];
            end
            nhex = nhex+1;
            node_index = node_index+1;
        end
        node_index = node_index+1;
    end
end

% create object r of class revolve_mesh. gcoord is 2d points (size Nx2), 3
% means 2d points are on xy plane, 2 means to be rotated wrt y-axis, [300,
% origin, 0] is the origin to rotate about, node_hex (Nx4) and node_tri
% (Nx3) are 2d node listing specifying elements.
out_of_bound = find( node_hex > length( gcoord ));
if ~isempty(out_of_bound)
    for i = 1:length(out_of_bound)
        [out_of_bound_(i,1), out_of_bound_(i,2)] = ind2sub( size(node_hex), out_of_bound(i));
    end
    to_remove = unique( out_of_bound_(:,1));
    node_hex( to_remove, :) = [];
end

e_tri = element(element_type.tri3, [1 2 3], 1);
e_tri(:) = [];
e_quad = e_tri;
for i = 1:length(node_tri)
    e_tri(i) = element(element_type.tri3, node_tri(i,:),1);
end
for i = 1:length(node_hex)
    e_quad(i) = element(element_type.quad4, node_hex(i,:),1);
end

if sphincter_wholewidth
    sphincterTriSet = 1:tri_cuttoff;
    sphincterQuadSet = 1:quad_cuttoff;
end

if sphincter_muscle
    for i = sphincterTriSet
        e_tri(i).materialID = 4;
    end
    for i = sphincterQuadSet
        e_quad(i).materialID = 3;
    end
end


mesh_obj_2D = mesh_class(gcoord, [e_tri e_quad]);

boundary_node_set = [boundary_node_top, fliplr(boundary_node_bot),1];

% below is for debugging
% figure, hold on; axis([0 600 0 300]); axis equal; set(gca,'ydir', 'reverse');
% for i = boundary_node_set
%     scatter(mesh_obj_2D.node_list(i,1), mesh_obj_2D.node_list(i,2), '.');
%     pause(0.3);
% end

end