classdef element_type
    properties
        n_nodes
        int_representation
        plot_order
        is_2D
    end
    
    properties (Dependent)
        is_3D
    end
    
    methods (Access = private)
        function obj = element_type(n_nodes, int_representation, plot_order, is_2D)
            obj.n_nodes = n_nodes;
            obj.int_representation = int_representation;
            obj.plot_order = plot_order;
            obj.is_2D = is_2D;
        end
    end
    
    methods
        function TF = get.is_3D(obj)
            TF = ~ obj.is_2D;
        end
        
        function sorted_array = sort(obj_array)
            [~, ix] = sort([obj_array.int_representation]);
            sorted_array = obj_array(ix);
        end
        
        function result = version3D(obj)
            switch obj
                case element_type.tri3
                    result = element_type.penta6;
                case element_type.quad4
                    result = element_type.hex8;
                otherwise
                    error('Error\n');
            end
        end
    end
    
    enumeration
        tri3(3,1, [1 2 3 1], 1)
        quad4(4,2,[1 2 3 4 1], 1)
        tet4(4,3, [1 2 3 1 4 2 3 4], 0)
        penta6(6,4, [1 2 3 1 4 5 6 4 5 2 3 6], 0)
        hex8(8,5, [1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4], 0)
    end
end

