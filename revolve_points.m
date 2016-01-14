classdef revolve_points
    %This class is used to revolve a set of 2D points about an axis in 360
    %degrees.
    
    properties (Access = protected)
        base3D
        plane  %plane of the 2D base
        axis   % axis of rotation
        origin
        num_replicates
    end
    
    properties
        sectional_angle_degree
        sectional_angle_radian
        num_basepoints
        num_totalpoints
    end
    
    methods
        function value = get.sectional_angle_degree(obj)
            value = 360.0/obj.num_replicates;
        end
        
        function value = get.sectional_angle_radian(obj)
            value = 2.0*pi / obj.num_replicates;
        end
        
        function value = get.num_basepoints(obj)
            value = size(obj.base3D, 1);
        end
        
        function value = get.num_totalpoints(obj)
            value = obj.num_basepoints * obj.num_replicates;
        end
    end
    
    methods
        function obj = revolve_points( base2D, plane, axis, origin, num_replicates)
            obj.plane = plane;
            obj.axis = axis;
            obj.origin = origin;
            obj.num_replicates = num_replicates;
            
            obj.base3D = zeros( size(base2D,1), 3);
            if plane == 1
                obj.base3D(:, 2:3) = base2D;
            elseif plane == 2
                obj.base3D(:, [3 1]) = base2D;
            elseif plane == 3
                obj.base3D(:, 1:2) = base2D;
            else error('the argument plane can only be 1, 2 or 3.');
            end
        end
        
        function result3D = generate3D(obj)
            result3D = zeros(obj.num_totalpoints,3);
            reference = zeros( size(obj.base3D));
            for i = 1:obj.num_basepoints
                reference(i,:) = obj.base3D(i,:) - obj.origin;
            end
            result3D(1:obj.num_basepoints,:) = reference( 1:obj.num_basepoints,:);
            for i = 2:obj.num_replicates
                angle = (i-1)* obj.sectional_angle_radian;
                result3D(1+(i-1)*obj.num_basepoints: i*obj.num_basepoints, obj.axis) = reference( :, obj.axis);
                third_axis = setdiff( [1,2,3], [obj.axis, obj.plane]);
                result3D(1+(i-1)*obj.num_basepoints: i*obj.num_basepoints, third_axis) = cos(angle)*reference( :, third_axis);
                result3D(1+(i-1)*obj.num_basepoints: i*obj.num_basepoints, obj.plane) = sin(angle)*reference( :, third_axis);
            end
            for i = 1:obj.num_totalpoints
                result3D(i,:) = result3D(i,:) + obj.origin;
            end
        end

    end
    
end

