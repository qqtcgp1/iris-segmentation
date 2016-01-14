function [ result, optimum ] = polygon_overlap( polygon1, polygon2, max_offset, algorithm )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if numel(max_offset) == 1
    max_offset = [max_offset, max_offset];
end
obj_fun = @(x) (-polygon_overlap_fix( polygon1 + repmat(x,[size(polygon1,1),1]), polygon2));

if nargin <= 3
    algorithm = 'patternsearch';
end

if polyarea( polygon1(:,1),polygon1(:,2)) >  polyarea( polygon2(:,1),polygon2(:,2))
    result = polygon_overlap( polygon2, polygon1, max_offset);
    return;
end


if strcmp(algorithm, 'patternsearch')
    optimstruct = struct('objective', obj_fun, 'x0', [0 0], 'Aineq', [], ...
        'bineq', [], 'Aeq', [], 'beq' ,[], 'lb', -max_offset, 'ub', max_offset, ...
        'nonlcon', [], 'rngstate', [], 'solver', 'patternsearch', 'options', []);
elseif strcmp(algorithm, 'ga')
    optimstruct = struct('fitnessfcn', obj_fun, 'nvars', 2, 'Aineq', [], ...
        'bineq', [], 'Aeq', [], 'beq', [], 'lb', -max_offset, 'ub', max_offset, ...
        'nonlcon', [], 'intcon', [], 'rngstate', [], 'solver', 'ga', 'options', []);
elseif strcmp(algorithm, 'de')
    optimstruct = struct('objective', obj_fun, 'I_D', 2, 'I_NP', 30, 'F_weight', 0.8, ...
        'FVr_minbound', [-10 -10], 'FVr_maxbound', [10 10], 'F_VTR', -0.99999, 'I_itermax', 5e9, ...
        'I_bnd_constr', 0, 'F_CR', 0.8, 'I_strategy', 1, 'solver', 'de');
end

polygon1(:,1) = polygon1(:,1) - min( polygon1(:,1));
polygon1(:,2) = polygon1(:,2) - min( polygon1(:,2));
polygon2(:,1) = polygon2(:,1) - min( polygon2(:,1));
polygon2(:,2) = polygon2(:,2) - min( polygon2(:,2));



rng('shuffle');
[optimum, result] = feval( str2func(optimstruct.solver), optimstruct);


end