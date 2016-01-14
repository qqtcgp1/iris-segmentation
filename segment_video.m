function [ segmented ] = segment_video( control, original, perferred_control, threshold, control_tol, verbose )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 6
    verbose = true;
end

if nargin < 5
    control_tol = 0.1;
end


for i = 1:length(original)
    try
        segmented{i} = gray2boundary_control(control, original{i}, threshold ,2, perferred_control,control_tol);
    catch
        segmented{i} = gray2boundary(original{i},30,2);
    end
    
    if verbose
        fprintf([num2str(i), 'th frame completed\n']);
    end
end


end