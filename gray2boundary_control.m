function [ boundary ] = gray2boundary_control( control, gray, graythresh, disksize, preferred_control, control_tol, recursion_indicator )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 7
    recursion_indicator = 'none';
end

if nargin <= 3
    disksize = 2;
end

boundary = (gray < graythresh);

% remove white and black holes
solid = imfill(boundary,'holes');
solid = ~imfill(~solid,'holes');

se = strel('disk', disksize);
solid = imclose(solid, se);
solid = ~imclose(~solid, se);

boundary = edge(solid,'canny',[],1);

min_control = preferred_control * (1-control_tol);
max_control = preferred_control * (1+control_tol);


[~, area, volume, ~] = get_iris(boundary);

if strcmp(control, 'area')
    control_val = area;
elseif strcmp(control, 'volume')
    control_val = volume;
end

while (control_val > max_control)
    if strcmp(recursion_indicator, 'control_too_small')
        warning('control_tol too small, cannot converge to preferred_area. control_tol has been increased by 0.01 for convergence.');
        boundary = gray2boundary_control(control, gray, graythresh, disksize, preferred_control, control_tol + 0.01);
        return;
    end
    boundary = gray2boundary_control(control, gray, graythresh + 1, disksize, preferred_control, control_tol, 'area_too_large');
    return;
end

while (control_val < min_control)
    if strcmp(recursion_indicator, 'area_too_large')
        warning('control_tol too small, cannot converge to preferred_control. control_tol has been increased by 0.01 for convergence.');
        boundary = gray2boundary_control(control, gray, graythresh, disksize, preferred_control, control_tol + 0.01);
        return;
    end
    boundary = gray2boundary_control(control, gray, graythresh - 1, disksize, preferred_control, control_tol, 'area_too_small');
    return;
end

boundary = edge(boundary,'canny',[],1);

end