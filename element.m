classdef element
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties
        type
        node_index
        materialID
    end
   
    properties (Dependent, SetAccess  = private, Hidden)
        is_2D
        is_3D
        plot_order
    end
        
    
    methods
        function obj = element( type, node_index, materialID )
            if nargin == 2
                materialID = 1;
            end
            obj.type = type;
            assert( isvector(node_index));
            if size(node_index,2) == 1
                node_index = node_index';
            end
            obj.node_index = node_index;
            obj.materialID = materialID;
        end
        
        function result = get.plot_order( obj )
            result = obj.type.plot_order;
        end
        
        function TF = get.is_2D(obj)
            TF = obj.type.is_2D;
        end
        
        function TF = get.is_3D(obj)
            TF = ~ obj.is_2D;
        end
        
        function TF = eq(obj, obj2)
            if numel(obj2) > numel(obj)
                TF = eq(obj2, obj);
                return;
            end
            
            if numel(obj) > 1
                if ( all(size(obj) == size(obj2)))
                    TF = arrayfun( @eq, obj, obj2);
                else assert( numel(obj2) == 1);
                    siz = size(obj);
                    obj3( 1:siz(1),1:siz(2)) = obj2;
                    TF = arrayfun( @eq, obj, obj3);
                end
                return
                    
            else
                TF = (obj.type == obj2.type) && all(sort(obj.node_index) == sort(obj2.node_index));
            end
        end
        
        function TF = ne(obj, obj2)
            TF = ~eq(obj,obj2);
        end
        
        function [sorted_array] = sort(obj_array, is_recursion)
            if nargin < 2
                is_recursion = 0;
            end
            
            if ~is_recursion
                % sort according to element type
                temp = [obj_array.type];
                [~,ia] = sort([temp.int_representation]);
                sorted_array = obj_array(ia);
                
                % sort each group of objects of the same type separately
                temp = [sorted_array.type];
                [~,ia,~] = unique( [temp.int_representation] );
                for i = 1:length(ia)-1
                    sorted_array(ia(i):ia(i+1)-1) = sort( sorted_array(ia(i):ia(i+1)-1), true );
                end
                sorted_array(ia(end):end) = sort( sorted_array(ia(end):end), true);
            else
                % is same element type. sort according to the node_index array
                % use of sortrows()
                matrix = zeros( numel( obj_array), obj_array(1).type.n_nodes);
                for i = 1:numel(obj_array)
                    matrix(i,:) = obj_array(i).node_index;
                end
                
                [~, index] = sortrows(matrix);
                sorted_array = obj_array(index);
            end
        end
            
            
        
        function obj = flip_orientation(obj)
            obj.node_index = flip( obj.node_index);
        end
    end
    
end

