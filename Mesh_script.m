N=1;
int = 10;
num_replicates = 20;
origin = 295.25;

% x coor space out by int size

distance_x = max(top_iris_r{N}(:,2)) - min(top_iris_r{N} (:,2));
n_int_x = nearest(distance_x / int);
n_node_x = n_int_x + 1;


h = 1;
for i =0 : n_int_x
    top_x(h) = min(top_iris_r{N}(:,2)) + nearest(i*distance_x / n_int_x);
   
    if top_x(h) > max( top_iris_r{N}(:,2) )
        top_x(h) = max( top_iris_r{N}(:,2));
    end
    
    top_y(h) = top_iris_r{N} ( find( top_iris_r{N} (:,2) == top_x(h), 1), 1);
    h = h+1;
end


distance_y = zeros(length(top_x),1); n_int_y = distance_y;

for i = 1: length(top_x)
    j = find( bottom_iris_r_smooth {N} (:,2) == top_x(i) );
    distance_y(i) = bottom_iris_r_smooth {N}(j,1)- top_y(i);
    
    if distance_y(i) >= int             
        n_int_y(i) = nearest(distance_y(i)/int);
    end
end


h =1;

gcoord_x(h) = tip_r{N}(1,2);
gcoord_y(h) = tip_r{N}(1,1);
h = h+1;

for i = 1: length(top_x)

    j = find( bottom_iris_r_smooth {N} (:,2) == top_x(i) );
    
    if (bottom_iris_r_smooth {N}(j,1) - top_y(i)) >= int
        
        for k = 0: n_int_y(i)
            gcoord_x(h) = top_x(i);
            gcoord_y (h) = top_y(i)+ nearest( k * distance_y(i) / n_int_y(i) );
            h = h+1;
        end
    end
end

gcoord = [gcoord_x;gcoord_y];
gcoord = gcoord';
for i = 1:length(n_int_y)
n_int_y (i) = n_int_y(i) + 1; 
end


ntri = 1;
nhex = 1; 
node_index = 1;
% triangle for the first element at the tip
node_tri (ntri, 1) = node_index;
node_tri (ntri, 2) = node_index+1;
node_tri (ntri, 3) = node_index+2;
ntri = ntri+1;
node_index = node_index+1;

for i = 2: length (n_int_y)-1
 if n_int_y (i+1) - n_int_y(i) == 2 
        node_tri (ntri, 1) = node_index;
        node_tri (ntri,2) = node_index + n_int_y(i); 
        node_tri (ntri,3) = node_tri(ntri,2) + 1; 
        ntri = ntri+1; 
   for j = 1 : n_int_y (i)-1
        node_hex (nhex, 1) = node_index; 
        node_hex (nhex,2) = node_index + n_int_y(i)+1;
        node_hex(nhex,3) = node_hex(nhex,2)+1; 
        node_hex(nhex,4) = node_hex(nhex,1) +1;
        nhex = nhex+1;
        node_index = node_index+1; 
    end
         node_tri(ntri,1) = node_index; 
         node_tri(ntri,2) = node_index + n_int_y (i+1) -1;
         node_tri(ntri,3) = node_tri(ntri,2) +1; 
         ntri =ntri+1;
         node_index = node_index +1; 
      end

 if n_int_y (i) - n_int_y(i+1) == 2 
        node_tri (ntri, 1) = node_index;
        node_tri (ntri,2) = node_index + n_int_y(i); 
        node_tri (ntri,3) = node_tri(ntri,1) + 1; 
        ntri = ntri+1; 
        node_index = node_index +1;
    for j = 1 : n_int_y (i+1)-1
        node_hex (nhex, 1) = node_index; 
        node_hex (nhex,2) = node_index + n_int_y(i) -1;
        node_hex(nhex,3) = node_hex(nhex,2)+1; 
        node_hex(nhex,4) = node_hex(nhex,1) +1;
        nhex = nhex+1;
        node_index = node_index+1; 
    end
         node_tri(ntri,1) = node_index; 
         node_tri(ntri,2) = node_index + n_int_y (i+1) + 1;
         node_tri(ntri,3) = node_tri(ntri,1) + 1; 
         ntri =ntri+1;
         node_index = node_index + 2; 
 end
      
  if n_int_y(i+1) - n_int_y (i) ==1 
          node_tri (ntri,1) = node_index;
          node_tri(ntri,2) = node_index + n_int_y(i);
          node_tri(ntri,3) = node_tri(ntri,2) +1;
          ntri= ntri+1;
          for j = 1 : n_int_y(i)-1
              node_hex (nhex, 1) = node_index;
              node_hex (nhex,2) = node_index + n_int_y(i)+1;
              node_hex(nhex,3) = node_hex(nhex,2)+1;
              node_hex(nhex,4) = node_hex(nhex,1) +1;
              nhex = nhex+1;
              node_index = node_index+ 1;
          end
           node_index = node_index+ 1;
  end
   
if n_int_y(i) - n_int_y (i+1) == 1 
          node_tri (ntri,1) = node_index;
          node_tri(ntri,2) = node_index + n_int_y(i);
          node_tri(ntri,3) = node_tri(ntri,1) +1;
          ntri= ntri+1;
          node_index = node_index+1;
          
          for j = 1 : n_int_y (i+1)-1
              node_hex (nhex, 1) = node_index;
              node_hex (nhex,2) = node_index + n_int_y(i)-1;
              node_hex(nhex,3) = node_hex(nhex,2)+1;
              node_hex(nhex,4) = node_hex(nhex,1) +1;
              nhex = nhex+1;
              node_index = node_index+1;
          end
         node_index = node_index+1;
end
      
if n_int_y (i) == n_int_y (i+1)
    for j = 1 : n_int_y (i)-1
              node_hex (nhex, 1) = node_index;
              node_hex (nhex,2) = node_index + n_int_y(i);
              node_hex(nhex,3) = node_hex(nhex,2)+1;
              node_hex(nhex,4) = node_hex(nhex,1) +1;
              nhex = nhex+1;
              node_index = node_index+1;
    end
    node_index = node_index+1;
end
end


% center = [tip_r{N}(1,2)-((tip_r{N}(1,2)-tip_l{N}(1,2))/2), tip_r{N}(1,1)];
% 
frame = zeros(300,600); figure, imshow(frame); hold on;
scatter(iris_contour_r{1}(:,2), iris_contour_r{1}(:,1), '.r');
hold on;
% scatter(gcoord_x,gcoord_y,'b*');



% plotting rectangle
for i =1:length(node_hex)
    line(gcoord(node_hex(i,[1 2 3 4 1]),1), gcoord(node_hex(i,[1 2 3 4 1]),2), 'color', 'r');
end



% plotting triangle
for i = 1:length(node_tri)
    line(gcoord(node_tri(i,[1 2 3 1]),1), gcoord(node_tri(i,[1 2 3 1]),2), 'color', 'g');
end

% create object r of class revolve_mesh. gcoord is 2d points (size Nx2), 3
% means 2d points are on xy plane, 2 means to be rotated wrt y-axis, [300,
% origin, 0] is the origin to rotate about, node_hex (Nx4) and node_tri
% (Nx3) are 2d node listing specifying elements.
r = revolve_mesh( gcoord, 3, 2, [300 origin 0], num_replicates, node_hex, node_tri);
[gcoord3D, hex_3D, tri_3D] = generate3D(r);

scatter3( gcoord3D(1:3*length(gcoord),1), gcoord3D(1:3*length(gcoord),2), gcoord3D(1:3*length(gcoord),3), 'xr');


for i = 1:length(hex_3D)
    line(gcoord3D(hex_3D(i,[1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4]),1), gcoord3D(hex_3D(i, [1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4]),2),...
        gcoord3D(hex_3D(i,[1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4]),3),'color', 'b' );
end

for i = 1:length(tri_3D)
    line(gcoord3D(tri_3D(i,[1 2 3 1 4 5 6 4 5 2 3 6]),1), gcoord3D(tri_3D(i, [1 2 3 1 4 5 6 4 5 2 3 6]),2),...
        gcoord3D(tri_3D(i,[1 2 3 1 4 5 6 4 5 2 3 6]),3),'color', 'b' );
end





