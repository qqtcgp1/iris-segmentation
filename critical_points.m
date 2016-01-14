 function [ meshwork, boundary] = critical_points( boundary, with_cornea, side , interested_fraction, threshold, eps)
%[a1, a2, a3, a4, a5, a6] = critical_points(iris_bw_image, threshold) gives a1 to a6
%to be the six corner and tip points of the iris image. threshold is
%defaulted to be 15.

if nargin < 6
    eps = 3;
end

if nargin < 5
    threshold = 10;
end

if nargin < 4
    interested_fraction = 0.5;
end

if nargin < 3
    side = 'left';
end

S = @(n) size(boundary,n);

if strcmp(side, 'right')
    [meshwork, boundary] = critical_points(fliplr(boundary), with_cornea, 'left', interested_fraction, threshold, eps);
    meshwork = fliplr_point(meshwork, S(2));
    if nargout == 2
        boundary = fliplr(boundary);
    end
    
    return;
end

errormsg = MException('critical_points:fail', 'cannot find critical points.');
interested = nearest(interested_fraction * S(2));



for i = 1: interested
    index = find(boundary(:,i));
    dif = zeros(length(index)-1,1);
    for j = 1:length(index)-1
        dif(j) = index(j+1) - index(j);
    end
    boundary(index(dif<=threshold),i) = 0;
end

num_edge = zeros(interested,1);

for i = 1:interested
    num_edge(i) = length(find(boundary(:,i)));
end


flag = true;
temp = 1;
meshwork(2) = 1;

while (flag)
    temp = temp + meshwork(2) - 1;
    meshwork = find(num_edge(temp:end)>=3, threshold);
    if numel(meshwork)~= threshold; throw(errormsg); end;
    flag = ~prod(diff(meshwork) == 1);
end
meshwork = meshwork(1) + temp - 1;

temp = find(boundary(:,meshwork));
meshwork = [temp(2), meshwork(1)];




 end