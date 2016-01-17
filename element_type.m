classdef element_type
    properties
        n_nodes
        int_representation
        plot_order
        dimensionality
    end
    
    properties (Dependent)
        is_1D
        is_2D
        is_3D
    end
    
    methods (Access = private)
        function obj = element_type(n_nodes, int_representation, plot_order, dimensionality)
            obj.n_nodes = n_nodes;
            obj.int_representation = int_representation;
            obj.plot_order = plot_order;
            obj.dimensionality = dimensionality;
        end
    end
    
    methods
        function TF = get.is_3D(obj)
            TF = (obj.dimensionality == 3);
        end
        
        function TF = get.is_2D(obj)
            TF = (obj.dimensionality == 2);
        end
        
        function TF = get.is_1D(obj)
            TF = (obj.dimensionality == 1);
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
                    
                case element_type.line2
                    result = element_type.quad4;
                    
                otherwise
                    error('Error\n');
            end
        end
    end
    
    enumeration
        tri3(3,1, [1 2 3 1], 2)
        quad4(4,2,[1 2 3 4 1], 2)
        tet4(4,3, [1 2 3 1 4 2 3 4], 3)
        penta6(6,4, [1 2 3 1 4 5 6 4 5 2 3 6], 3)
        hex8(8,5, [1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4], 3)
        
        line2(2,6, [1 2], 1);
    end
end

