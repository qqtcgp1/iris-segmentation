N=1;
int = 5;

% x coor space out by int size
distance_x = max(top_iris_r{N}(:,2)) - min(top_iris_r{N} (:,2));
n_int_x = nearest(distance_x / int);
n_node_x = n_int_x + 1;


% detailed x spacing
h = 1;
for i =0 : n_int_x
    top_x(h) = min(top_iris_r{N}(:,2)) + nearest(i*distance_x / n_int_x);
   
    if top_x(h) > max( top_iris_r{N}(:,2) )
        top_x(h) = max( top_iris_r{N}(:,2));
    end
    
    top_y(h) = top_iris_r{N} ( find( top_iris_r{N} (:,2) == top_x(h), 1), 1);
    h = h+1;
end


% detailed y spacing
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
nrec = 1; 
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
        node_rec (nrec, 1) = node_index; 
        node_rec (nrec,2) = node_index + n_int_y(i)+1;
        node_rec(nrec,3) = node_rec(nrec,2)+1; 
        node_rec(nrec,4) = node_rec(nrec,1) +1;
        nrec = nrec+1;
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
        node_rec (nrec, 1) = node_index; 
        node_rec (nrec,2) = node_index + n_int_y(i) -1;
        node_rec(nrec,3) = node_rec(nrec,2)+1; 
        node_rec(nrec,4) = node_rec(nrec,1) +1;
        nrec = nrec+1;
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
              node_rec (nrec, 1) = node_index;
              node_rec (nrec,2) = node_index + n_int_y(i)+1;
              node_rec(nrec,3) = node_rec(nrec,2)+1;
              node_rec(nrec,4) = node_rec(nrec,1) +1;
              nrec = nrec+1;
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
              node_rec (nrec, 1) = node_index;
              node_rec (nrec,2) = node_index + n_int_y(i)-1;
              node_rec(nrec,3) = node_rec(nrec,2)+1;
              node_rec(nrec,4) = node_rec(nrec,1) +1;
              nrec = nrec+1;
              node_index = node_index+1;
          end
         node_index = node_index+1;
end
      
if n_int_y (i) == n_int_y (i+1)
    for j = 1 : n_int_y (i)-1
              node_rec (nrec, 1) = node_index;
              node_rec (nrec,2) = node_index + n_int_y(i);
              node_rec(nrec,3) = node_rec(nrec,2)+1;
              node_rec(nrec,4) = node_rec(nrec,1) +1;
              nrec = nrec+1;
              node_index = node_index+1;
    end
    node_index = node_index+1;
end
end


center = [tip_r{N}(1,2)-((tip_r{N}(1,2)-tip_l{N}(1,2))/2), tip_r{N}(1,1)];




























