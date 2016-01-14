function intercept = iris_cornea_intercept( contour, extended_cornea, meshwork, side )

if nargin <= 2
    side = 'left';
end

if ~iscell(contour)
    siz = size(contour);
    contour = {bw2points(contour), siz };
end

S = contour{2};
contour = contour{1};


if strcmp(side, 'right')
    fliplr_point_modified = @( point2 ) ( S(2) + 1 - point2 );
    
    contour(:,2) = arrayfun(fliplr_point_modified, contour(:,2));
    extended_cornea(:,2) = arrayfun(fliplr_point_modified, extended_cornea(:,2));
    
    intercept = fliplr_point(iris_cornea_intercept({contour, S}, flip(extended_cornea), fliplr_point(meshwork, S(2)), 'left'), S(2));
    return;
end



for j = meshwork(2):-1:1
    column = find(contour(:,2) == j);
    flag = false;
    
    if isempty(column); continue;
        
    else
        if extended_cornea(j) > contour(column(end), 1) && ...
            numel(find_path(points2bw(contour, S), meshwork, [contour(column(end),1), j])) > 100
            flag = true;
        end
    end
    
    if flag
        intercept = [contour(column(end), 1), j]; break;
    end;
end


end