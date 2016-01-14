classdef LS_Problem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        A;
        b;
    end
    
    properties (Dependent)
        n_equations;
        n_unknowns;
    end
    
    methods
        function value = get.n_equations(obj)
            value = size( obj.A, 1);
        end
        
        function value = get.n_unknowns(obj)
            value = size(obj.A, 2);
        end
        
        function [] = check_dim(obj)
            assert(isvector(obj.b));
            assert( length(obj.b) == obj.n_equations );
            assert( obj.n_equations >= obj.n_unknowns );
        end
        
        function obj = LS_Problem( A, b )
            if nargin == 0
                A = []; b = [];
            end
            obj.A = A; obj.b = b;
        end
        
        function x = solve(obj)
            check_dim(obj);
            x = obj.A \ obj.b;
        end
        
        function constrained = add_constraint(obj, v, y)
            constrained = constrained_LS( obj, v, y );
        end
        
    end
    
end

