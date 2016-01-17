classdef revolve_mesh < revolve_points
    %This is used to revolve a set of points, along with the meshed elements,
    %about an axis. The original 2D elements are modified to 3D elements.
    
    properties (Access = private)
        elements_2D
    end
    
    properties (Dependent)
        num_elem_2D
        num_elem_3D
    end
    
    methods
        function value = get.num_elem_2D(obj)
            value = numel(obj.elements_2D);
        end
        
        function value = get.num_elem_3D(obj)
            value = obj.num_elem_2D * obj.num_replicates;
        end
    end               
    
    methods
        function obj = revolve_mesh(mesh_obj, plane, axis, origin, num_replicates)
            base2D = mesh_obj.node_list;
            elements_ = mesh_obj.elements;
            obj = obj@revolve_points(base2D, plane, axis, origin, num_replicates);
            obj.elements_2D = elements_;
        end

        
        function mesh_3D = generate3D(obj)
            nodes = generate3D@revolve_points(obj);
            elements_3D = obj.elements_2D(1);
            for i = 1:obj.num_elem_3D
                elements_3D(i) = elements_3D(1);
            end
           
            function e_out = update_elements( e_in, index )
                if numel(e_in) == 1
                    index2 = index + 1;
                    if index2 == obj.num_replicates
                        index2 = 0;
                    end
                    e_out = element( version3D(e_in.type), ...
                        [e_in.node_index + index*obj.num_basepoints, ...
                        e_in.node_index + index2* obj.num_basepoints], ...
                        e_in.materialID);
                else
                    e_out = e_in;
                    for j = 1:numel(e_in)
                        e_out(j) = update_elements( e_in(j), index);
                    end
                end
            end
            
            for i = 1:obj.num_replicates - 1
                elements_3D(1+(i-1)*obj.num_elem_2D: i*obj.num_elem_2D) = ...
                    update_elements( obj.elements_2D, i-1);
            end
            
            elements_3D( 1+(obj.num_replicates-1)*obj.num_elem_2D: obj.num_replicates*obj.num_elem_2D) = ...
                update_elements( obj.elements_2D, obj.num_replicates - 1);
            mesh_3D = mesh_class( nodes, elements_3D);
        end
    end
end

