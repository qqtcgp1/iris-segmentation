
% 3D mesh node assign
angle = 5/180 *pi; 
nnode = length (gcoord_x);
node_index3d = nnode +1;
last_node = nnode*2;

for i = 1: nnode
    gcoord_z (i) = 0;
end

h = 1;

for i = node_index3d:last_node
    gcoord_x(i) = gcoord_x(h);
    gcoord_y(i) = gcoord_y(h);
    
    distance_adj = sqrt(((center(1)- gcoord_x(h))^2)+(center(2)-gcoord_y(h))^2);
    
    distance_oppo = distance_adj *angle; 
    
    gcoord_z(i) = distance_oppo; 
    
    h = h+1;
end

gcoord = [gcoord_x;gcoord_y; gcoord_z];
gcoord = gcoord';

% frame = zeros(300,600); figure, imshow(frame); hold on;
% scatter(iris_contour_r{1}(:,2), iris_contour_r{1}(:,1), '.r');
plot3(gcoord_x,gcoord_y,gcoord_z,'b*');
axis ([0 600 0 300 0 200]);
hold on;

% 
%  for i = 1: length (gcoord_x)
%      text(gcoord_x(i),gcoord_y(i),gcoord_z(i), int2str(i),'fontsize',8,'color','r');
%  end

% 3D elements
nnel = size(node_rec,2);  %no. of node per element
nel = length(node_rec); %no. of elements
node_hex = zeros(nel,8);
for i = 1:nel
    for j= 1: nnel
    node_hex(i,j) = node_rec(i,j);
    node_hex(i,j+4) = node_rec(i,j)+nnode;
    end
end

nnel = size (node_tri,2);
nel = length(node_tri);

for i = 1:length(node_tri)
    for j = 1: nnel
   node_tetra (i,j)=node_tri(i,j);
   node_tetra(i,j+3)=node_tri(i,j) + nnode;
    end
end

%plot hexahedral element 3D

nel = length(node_hex) ;                  % number of elements
nnode = length(gcoord) ;          % total number of nodes in system
nnel = size(node_hex,2);                % number of nodes per element
% 
% Initialization of the required matrices
X_3d = zeros(nnel,nel) ;
Y_3d = zeros(nnel,nel) ;
Z_3d = zeros(nnel,nel) ;

    
for iel=1:nel   
     for i=1:nnel
     nd(i)=node_hex(iel,i);         % extract connected node for (iel)-th element
     X_3d(i,iel)=gcoord(nd(i),1);    % extract x value of the node
     Y_3d(i,iel)=gcoord(nd(i),2);    % extract y value of the node
     Z_3d(i,iel)=gcoord(nd(i),3) ;   % extract z value of the node
     end  
end
    
%plot hexahedral element    
    
for i = 1: length(X_3d)
    plot3 ([X_3d(1,i), X_3d(2,i)], [Y_3d(1,i),Y_3d(2,i)],[Z_3d(1,i),Z_3d(2,i)]);
    plot3 ([X_3d(2,i),X_3d(3,i)],[Y_3d(2,i),Y_3d(3,i)],[Z_3d(2,i),Z_3d(3,i)]);
    plot3([X_3d(3,i),X_3d(4,i)],[Y_3d(3,i),Y_3d(4,i)],[Z_3d(3,i),Z_3d(4,i)]);
    plot3([X_3d(4,i),X_3d(1,i)],[Y_3d(4,i),Y_3d(1,i)],[Z_3d(4,i),Z_3d(1,i)]);
    plot3 ([X_3d(1+4,i), X_3d(2+4,i)], [Y_3d(1+4,i),Y_3d(2+4,i)],[Z_3d(1+4,i),Z_3d(2+4,i)]);
    plot3 ([X_3d(2+4,i),X_3d(3+4,i)],[Y_3d(2+4,i),Y_3d(3+4,i)],[Z_3d(2+4,i),Z_3d(3+4,i)]);
    plot3([X_3d(3+4,i),X_3d(4+4,i)],[Y_3d(3+4,i),Y_3d(4+4,i)],[Z_3d(3+4,i),Z_3d(4+4,i)]);
    plot3([X_3d(4+4,i),X_3d(1+4,i)],[Y_3d(4+4,i),Y_3d(1+4,i)],[Z_3d(4+4,i),Z_3d(1+4,i)]);
    
    plot3 ([X_3d(1,i),X_3d(1+4,i)], [Y_3d(1,i),Y_3d(1+4,i)],[Z_3d(1,i),Z_3d(1+4,i)]);
    plot3 ([X_3d(2,i),X_3d(2+4,i)], [Y_3d(2,i),Y_3d(2+4,i)],[Z_3d(2,i),Z_3d(2+4,i)]);
    plot3 ([X_3d(3,i),X_3d(3+4,i)], [Y_3d(3,i),Y_3d(3+4,i)],[Z_3d(3,i),Z_3d(3+4,i)]);
    plot3 ([X_3d(4,i),X_3d(4+4,i)], [Y_3d(4,i),Y_3d(4+4,i)],[Z_3d(4,i),Z_3d(4+4,i)]);
end



    
    
nel = length(node_tetra) ;                  % number of elements
nnode = length(gcoord) ;          % total number of nodes in system
nnel = size(node_tetra,2);                % number of nodes per element
%
% Initialization of the required matrices
x_3d = zeros(nnel,nel) ;
y_3d = zeros(nnel,nel) ;
z_3d = zeros(nnel,nel) ;


for iel=1:nel
    for i=1:nnel
        nd(i)=node_tetra(iel,i);         % extract connected node for (iel)-th element
        x_3d(i,iel)=gcoord(nd(i),1);    % extract x value of the node
        y_3d(i,iel)=gcoord(nd(i),2);    % extract y value of the node
        z_3d(i,iel)=gcoord(nd(i),3) ;   % extract z value of the node
    end
end
    
    
for i = 1: length(x_3d)
    plot ([x_3d(1,i), x_3d(2,i)], [y_3d(1,i),y(2,i)]);
    plot ([x_3d(2,i),x_3d(3,i)],[y(2,i),y(3,i)]);
    plot([x_3d(3,i),x_3d(1,i)],[y(3,i),y(1,i)]);
end


        
for i = 1: length(x_3d)
    plot3 ([x_3d(1,i), x_3d(2,i)], [y_3d(1,i),y_3d(2,i)],[z_3d(1,i),z_3d(2,i)]);
    plot3 ([x_3d(2,i),x_3d(3,i)],[y_3d(2,i),y_3d(3,i)],[z_3d(2,i),z_3d(3,i)]);
    plot3([x_3d(3,i),x_3d(1,i)],[y_3d(3,i),y_3d(1,i)],[z_3d(3,i),z_3d(1,i)]);
    
    plot3 ([x_3d(1+3,i), x_3d(2+3,i)], [y_3d(1+3,i),y_3d(2+3,i)],[z_3d(1+3,i),z_3d(2+3,i)]);
    plot3 ([x_3d(2+3,i),x_3d(3+3,i)],[y_3d(2+3,i),y_3d(3+3,i)],[z_3d(2+3,i),z_3d(3+3,i)]);
    plot3([x_3d(3+3,i),x_3d(4,i)],[y_3d(3+3,i),y_3d(4,i)],[z_3d(3+3,i),z_3d(4,i)]);
    
    plot3 ([x_3d(1,i),x_3d(1+3,i)], [y_3d(1,i),y_3d(1+3,i)],[z_3d(1,i),z_3d(1+3,i)]);
    plot3 ([x_3d(2,i),x_3d(2+3,i)], [y_3d(2,i),y_3d(2+3,i)],[z_3d(2,i),z_3d(2+3,i)]);
    plot3 ([x_3d(3,i),x_3d(3+3,i)], [y_3d(3,i),y_3d(3+3,i)],[z_3d(3,i),z_3d(3+3,i)]);
    
end

    
    
    
    
    
    
    
    
    
    
    
    
    