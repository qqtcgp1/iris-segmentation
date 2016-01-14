classdef constrained_LS < LS_Problem
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        constrained_coeff;
        constrained_rhs;
    end
    
    methods
        function obj = constrained_LS( LS, constrained_coeff, constrained_rhs)
            if nargin == 0
                new_A = []; new_b = []; constrained_coeff = []; constrained_rhs = [];
            else
            
            A1 = LS.A(:,1);
            A_prime = LS.A(:,2:end);
            new_b = LS.b - A1*constrained_rhs / constrained_coeff(1);
            new_A = A_prime - A1/constrained_coeff(1) * row_vector(constrained_coeff(2:end));
            end
            
            obj@LS_Problem(new_A, new_b);
            obj.constrained_coeff = constrained_coeff;
            obj.constrained_rhs = constrained_rhs;
        end
        
        function result = redundant_coeff(obj, solution)
            result = (obj.constrained_rhs - dot( obj.constrained_coeff(2:end), solution)) / obj.constrained_coeff(1);
        end
        
        function x = solve(obj)
            x_ = solve@LS_Problem (obj);
            %x = (obj.constrained_rhs - dot( obj.constrained_coeff(2:end), x_ )) / obj.constrained_coeff(1);
            x = [redundant_coeff(obj,x_); x_];
        end
        
    end
    
end


function v = row_vector( vector )

assert( isvector(vector) );

if size(vector,1)~=1
    v = vector';
else v = vector;
end
end

