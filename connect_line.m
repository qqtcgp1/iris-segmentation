function [ line, bw_w_line ] = connect_line( first, last, bw )
%line = connect_line(first, last) find the line (as columns of coordinates)
%connecting first and last. input bw and output bw_w_line are optional. (used when
%wishing to add the line to the bw image.) Coordinates in 'line' will only
%be turned into integers if input bw is provided.

if first(1) > last(1)
    temp = first; first = last; last = temp;
end

slope = (last(2) - first(2) )/(last(1) - first(1) );

if abs(slope) > 1.0
    if nargin == 3
        line_ = connect_line( [first(2), first(1)], [last(2), last(1)], bw');
    end
    if nargin == 2
        line_ = connect_line( [first(2), first(1)], [last(2), last(1)]);
    end
    
    line = [line_(:,2), line_(:,1)];
    if nargout == 2
        assert(nargin == 3);
        line_index = sub2ind_array(size(bw), line);
        bw_w_line = bw;
        bw_w_line(line_index) = 1;
    end
    return;    
end
    
line = zeros( last(1) - first(1) + 1, 2);


for i = first(1):last(1)
    line(i - first(1) + 1,:) = [i, first(2) + slope* ( i-first(1)) ];
end

if nargin == 3
    line(:,1) = nearest_(line(:,1), size(bw,1));
    line(:,2) = nearest_(line(:,2), size(bw,2));
end
    
if nargout == 2
    assert( nargin == 3);
    line_index = sub2ind_array( size(bw), line);
    bw_w_line = bw;
    bw_w_line(line_index) = 1;
end

end