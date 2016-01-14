

S = size(frame); size1 = S(1); size2 = S(2); clear S
cont_ = cont;
for i = 1:size2
    index = find(cont(:,i));
    dif = zeros(length(index)-1,1);
    for j = 1:length(index)-1
        dif(j) = index(j+1) - index(j);
    end
    cont_(index(find(dif<=15) + 1),i) = 0;
end

num_edge = zeros(size2,1);

for i = 1:size2
    num_edge(i) = length(find(cont_(:,i)));
end


flag = true;
while (flag)
    mesh_left = find(num_edge>=3,3);
    flag = (mesh_left(3) - mesh_left(2) ~=1 ) || (mesh_left(2) - mesh_left(1) ~=1);
end
mesh_left = mesh_left(1);

flag = true;
while (flag)
    tip_left = find(num_edge(mesh_left:end) <= 2, 3);
    flag = (tip_left(3) - tip_left(2) ~=1 ) || (tip_left(2) - tip_left(1) ~= 1);
end
tip_left = mesh_left + tip_left(1);

flag = true;
while (flag)
    tip_right = find(num_edge(tip_left:end) >= 3, 3);
    flag = (tip_right(3) - tip_right(2) ~=1) || (tip_right(2) - tip_right(1) ~=1);
end
tip_right = tip_left + tip_right(1);

flag = true;
while (flag)
    mesh_right = find(num_edge(tip_right:end) <= 2, 3);
    flag = (mesh_right(3) - mesh_right(2) ~=1) || (mesh_right(2) - mesh_right(1) ~=1);
end
mesh_right = mesh_right(1) + tip_right;


tip_left = tip_left - 3;
mesh_right = mesh_right - 3;

temp = find(cont_(:,mesh_left)); mesh_left = [temp(2), mesh_left];
temp = find(cont_(:,tip_left)); tip_left = [temp(length(temp)), tip_left];
temp = find(cont_(:,tip_right)); tip_right = [temp(length(temp)), tip_right];
temp = find(cont_(:,mesh_right)); mesh_right = [temp(2), mesh_right];


figure, imshow(cont);
hold on;
 scatter([mesh_left(2) tip_left(2) tip_right(2) mesh_right(2)],...
    [mesh_left(1) tip_left(1) tip_right(1) mesh_right(1)], 'r', 'o', 'linew', 5);