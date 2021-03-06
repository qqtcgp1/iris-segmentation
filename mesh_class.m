classdef mesh_class
    %mesh_class is a type containing meshing information needed for a
    % geometry. A list of nodes are stored in node_list. Different 2D/3D
    % element types can be stored. Elements are stored by referencing to
    % a few nodes in node_list. An element can be referenced by an index
    % specifying the overall position among all the element types, using
    % the function ElementAccess. Plot and flip_orientation functions are
    % available. We can also output the mesh information into text with the
    % function print_feb, in feb file format.
    
    properties (Access = public)
        node_list
        elements
    end
    
    properties (Dependent, SetAccess = 'private')
        dimensionality
        is3D
        is2D
        is1D
        num_nodes
        num_elements
    end
    
    methods    
        function D = get.dimensionality(obj)
            temp = [obj.elements.type];
            D = max( [temp.dimensionality] );
        end
            
        function obj = mesh_class(node_list, elements)          
            obj.node_list = node_list;
            obj.elements = elements;
        end
        
            
        
        function TF = get.is2D(obj)
            TF = (obj.dimensionality == 2);
        end
        
        function TF = get.is3D(obj)
            TF = (obj.dimensionality == 3);
        end
        
        function TF = get.is1D(obj)
            TF = (obj.dimensionality == 1);
        end
        
        function result = get.num_nodes(obj)
            result = length(obj.node_list);
        end
        
        function result = get.num_elements(obj)
            result = length(obj.elements);
        end
        
        
        function e = ElementAccess(obj,a)
            e = obj.elements(a);
            e = e.node_index;
        end
        
        
        function [] = plot(obj, markertype, lineproperties)
            if nargin < 3
                lineproperties = { 'color', 'b' };
            end
            
            if nargin < 2
                markertype = '.r';
            end
            
            if obj.is3D
                scatter3( obj.node_list(:,1), obj.node_list(:,2), obj.node_list(:,3), markertype);
            else
                scatter( obj.node_list(:,1), obj.node_list(:,2), markertype);
            end
            
            hold on

             if size(obj.node_list,2) == 2
                 red = [];
                for i = 1:length(obj.elements)
                    e = obj.elements(i);
                    if e.materialID == 1
                        line( obj.node_list(e.node_index(e.plot_order),1), ...
                            obj.node_list( e.node_index(e.plot_order), 2), lineproperties{:});
                    else red = [red, i];
%                     elseif e.materialID == 2
%                         line( obj.node_list(e.node_index(e.plot_order),1), ...
%                             obj.node_list( e.node_index(e.plot_order), 2), 'color', 'y');
%                     else 
%                         line( obj.node_list(e.node_index(e.plot_order),1), ...
%                             obj.node_list( e.node_index(e.plot_order), 2), 'color', 'r');
                    end
                end
                
                if ~isempty(red)
                    for e = obj.elements(red)
                        line( obj.node_list(e.node_index(e.plot_order),1), ...
                            obj.node_list( e.node_index(e.plot_order), 2), 'color', 'r');
                    end
                end
             
            elseif size(obj.node_list, 2) == 3
                red = [];
                for i = 1:length(obj.elements)
                    e = obj.elements(i);
                    if e.materialID == 1
                        line( obj.node_list( e.node_index(e.plot_order), 1), ...
                            obj.node_list( e.node_index(e.plot_order), 2), ...
                            obj.node_list( e.node_index(e.plot_order), 3), lineproperties{:});
%                     elseif e.materialID == 2
%                         line( obj.node_list(e.node_index(e.plot_order),1), ...
%                             obj.node_list( e.node_index(e.plot_order), 2), ...
%                             obj.node_list( e.node_index(e.plot_order), 3), 'color', 'y');
%                     else 
%                         line( obj.node_list(e.node_index(e.plot_order),1), ...
%                             obj.node_list( e.node_index(e.plot_order), 2), ...
%                             obj.node_list( e.node_index(e.plot_order), 3), 'color', 'r');
                    else red = [red, i];
                    end
                end
                if ~isempty(red)
                    for e = obj.elements(red)
                        line( obj.node_list(e.node_index(e.plot_order),1), ...
                            obj.node_list( e.node_index(e.plot_order), 2), ...
                            obj.node_list( e.node_index(e.plot_order), 3),'color', 'r');
                    end
                end
            end
        end
        
        function obj = flip_orientation(obj)
            for i = 1: numel(obj.elements)
                obj.elements(i).node_index = flip(obj.elements(i).node_index);
            end
        end
        
        
        
        function obj_out = revolve(obj, num_replicates, plane, axis, origin)
            revolve_obj = revolve_mesh( obj.node_list, plane, axis, origin, ...
                num_replicates, one_type(obj.elements, element_type.tri3), ...
                one_type(obj.elements, element_type.quad4));
            
            [node_3D, hex3D, tri3D] = generate3D(revolve_obj);
            
            for i = 1:length(hex3D)
                elements_(i) = element( element_type.hex8, hex3D(i,:));
            end
            
            for i = length(hex3D) + 1: length(hex3D) + length(tri3D)
                elements_(i) = element( element_type.penta6, tri3D( i-length(hex3D),:));
            end
            
            obj_out = mesh_class(node_3D, elements_);
        end
        
        function result = one_type(obj, e_type)
            index = [];
            for j = 1:length(obj.elements)
                if obj.elements(j).type == e_type
                    index = [index j];
                end
            end
            result_ = obj.elements(index);
            result = zeros(length(result_), e_type.n_nodes);
            
            for j = 1:length(result_)
                result(j,:) = result_(j).node_index;
            end
        end
        
        
        
        function obj = plus(obj1, obj2)
            [ia,~] = ismembertol( obj2.node_list, obj1.node_list, 0.0001, 'ByRows', true);
            node_sum = [obj1.node_list; obj2.node_list( ia==0, :)];
            
            function e_out = update_element(e_in)
                if numel( e_in) > 1
                    e_out = e_in;
                    for i = 1:numel(e_in)
                        e_out(i) = update_element(e_in(i));
                    end
                else
                    e_out = element(e_in.type, arrayfun( @update_nodeid, e_in.node_index), e_in.materialID);
                end
            end
            
            function y = update_nodeid(x)
                if isempty(x)
                    y = []; return;
                end
                [~, y] = ismembertol( obj2.node_list(x,:), node_sum, 0.0001, 'ByRows',true);
            end
            
            element_sum = unique([obj1.elements, update_element(obj2.elements )]);

            obj= mesh_class(node_sum, element_sum);
        end
            
        
        function [] = print_feb(obj, fileID)
            fprintf(fileID, '\t<Geometry>\n\t\t<Nodes>\n');
            for i = 1:obj.num_nodes
                fprintf(fileID, '\t\t\t<node id="%i"> ', i);
                print_vector_feb( obj.node_list(i,:), fileID);
                fprintf(fileID, '</node>\n');
            end
            
            fprintf(fileID, '\t\t</Nodes>\n\t\t<Elements>\n');
            
            for i = 1:length(obj.elements)
                fprintf(fileID,['\t\t\t<', char(obj.elements(i).type), ...
                    ' id="%i" mat="%i"> '], i, obj.elements(i).materialID);
                print_vector_feb( obj.elements(i).node_index, fileID);
                fprintf(fileID, ['</' char(obj.elements(i).type) '>\n']);
            end            
            fprintf(fileID, '\t\t</Elements>\n\t</Geometry>\n');
        end
    end
end


function [ ] = print_vector_feb( vect , fileID)

for x = vect(1:end-1)
    fprintf(fileID, '%g, ', x);
end
fprintf( fileID, '%g',vect(end));
end
