clear
UPPER_CORNEA_FIT_ORDER = 6;
LOWER_CORNEA_FIT_ORDER = 6;
LOWER_IRIS_FIT_ORDER = 5;
N = 25;



crop_posi1 = []; crop_posi2 = [];
meshwork_appro_l = zeros(N,2); meshwork_appro_r = meshwork_appro_l;
h_fig = meshwork_appro_l(:,1); h_axes = h_fig;

incorrect_path = [];

for i = 1:N
    % 2.Read the frame
    frame = imread(['/Users/yuquan/Documents/research works/Iris/Dewarp/BEIYAN WANG_OS01/00', num2str(41+i), '_dewarp.tif']);
    
    % 4.turn to gray scale
    frame = rgb2gray(frame);
    
    % 5.remove the mark on top-right
    frame = remove_markings_Dewarp(frame, 3);
    
    
    % normalize the intensity (by the intensity of the first frame.
    frame = double(frame);
    mean_intensity(i) = mean(mean(frame));
    
    original{i} = frame;
    normalized{i} = frame *( mean_intensity(1) / mean_intensity(i) );
    
    
    % Through an GUI, get rid of the interference of the cneter shadow.
    % Also get cornea boundary via threshold method.
    [boundary_cornea{i}, crop_posi1, crop_posi2] = fix_center_shadow(crop_posi1, crop_posi2, normalized{i}, 5, 2);
    
    h_fig(i) = figure; imshow(boundary_cornea{i});
    h_axes(i) = gca;
    
    
    % approximated position of the meshwork. More accurate positions to be
    % determined later when the iris contour is found.
    meshwork_appro_l(i,:) = critical_points(boundary_cornea{i}, true, 'left', 0.5, 15, 3);
    meshwork_appro_r(i,:) = critical_points(boundary_cornea{i}, true, 'right', 0.5, 15, 3);
    
    %     ml = critical_points(boundary{i}, false, 'left', 0.6, 10, 3);
    %     mr = critical_points(boundary{i}, false, 'right', 0.6, 10, 3);
    
    hold on;
    A = [meshwork_appro_l(i,:); meshwork_appro_r(i,:)];
    %scatter(A(:,2), A(:,1), 'ro', 'linew', 8);
    
    
    % so lower cornea boundary is defined by the path between appro meshwork
    % points
    lower_cornea_ind{i} = find_path(boundary_cornea{i}, meshwork_appro_l(i,:), meshwork_appro_r(i,:));
    lower_cornea_sub{i} = ind2sub_array(size(frame), lower_cornea_ind{i});
    
    % When find_path fails, store the index i into incorrect_path
    if length(lower_cornea_ind{i}) > 650
        incorrect_path = [incorrect_path, i];
    end
    
    temp = lower_cornea_sub{i};
    scatter(temp(:,2), temp(:,1), '.b');
end


% set of indices of correct_path
correct_path = setdiff(1:N, incorrect_path, 'stable');
graythresh = [];


% loops again
for i = 1:N
    
    % using GUI to determine gray threshold value
    [boundary_iris{i}, graythresh] = gray2boundary(normalized{i}, graythresh);
    
    % labelling connected componets
    label = bwlabel(boundary_iris{i});
    points = bw2points(boundary_iris{i});
    
    
    % if more than 3 connected components, only keep the largest two (which
    % will be the left and right part of iris).
    if max(max(label)) >= 3
        boundary_iris{i}(:) = false;
        [~,ind] = min(points(:,2));
        ind = points(ind,:);
        boundary_iris{i} (label == label(ind(1), ind(2))) = true;
        [~,ind] = max(points(:,2));
        ind = points( ind, :);
        boundary_iris{i} (label == label(ind(1), ind(2))) = true;
    end
    
    
    if ismember(i, incorrect_path)
        % if incorrect_path (cornea contour not found), then use the last
        % correct path to approximate
        previous = find( correct_path < i, 1,'last');
        fit_coeff{i} = fit_coeff{previous};
        extended_cornea_dense{i} = extended_cornea_dense{previous};
        
    else
        % fit the cornea contour into a polynomial of order
        % LOWER_CONREA_FIT_ORDER
        raw_cornea = lower_cornea_sub{i};
        x = raw_cornea (:,2);
        y_raw = raw_cornea(:,1);
        
        [curve, goodness] = fit(x, y_raw, ['poly' num2str(LOWER_CORNEA_FIT_ORDER)]);
        fit_coeff{i} = coeffvalues(curve);
    end
    
    
    % extend the cornea polynomial to the whole horizontal range
    extended_cornea{i} = polyval(fit_coeff{i}, 1:size(frame,2));
    extended_cornea{i} = nearest(extended_cornea{i});
    
    % find meshwork more accurately
    [meshwork_l{i}] = critical_points(boundary_iris{i}, false, 'left', 0.55, 15, 3);
    [meshwork_r{i}] = critical_points(boundary_iris{i}, false, 'right', 0.55, 15, 3);
    
    boundary_temp = bw2points(boundary_iris{i});
    
    
    % lower cutoff of the iris = intercept of extended cornea and
    % iris_contour
    intercept_l{i} = iris_cornea_intercept( {boundary_temp, size(boundary_iris{i}) }, extended_cornea{i}, meshwork_l{i}, 'left');
    intercept_r{i} = iris_cornea_intercept( {boundary_temp, size(boundary_iris{i}) }, extended_cornea{i}, meshwork_r{i}, 'right');
    
    % erase lines that are well above extended_cornea
    for k = 1:size(boundary_iris{i}, 2)
        boundary_iris{i}( 1: min(floor(extended_cornea{i}(k))  - 15, size(boundary_iris{i},1)), k) = 0;
    end
    boundary_temp = bw2points(boundary_iris{i});
    
    
    iris_contour_l{i} = find_path(boundary_iris{i}, intercept_l{i}, meshwork_l{i});
    iris_contour_l{i} = ind2sub_array(size(boundary_iris{i}), iris_contour_l{i});
    
    iris_contour_r{i} = find_path(boundary_iris{i}, intercept_r{i}, meshwork_r{i});
    iris_contour_r{i} = ind2sub_array(size(boundary_iris{i}), iris_contour_r{i});
    
    [~,ind] = max(iris_contour_l{i}(:,2));
    tip_l{i} = iris_contour_l{i}(ind,:);
    [~,ind] = min(iris_contour_r{i}(:,2));
    tip_r{i} = iris_contour_r{i}(ind,:);
    
    bottom_iris_l{i} = find_path( points2bw(iris_contour_l{i}, size(boundary_iris{i})), tip_l{i}, intercept_l{i}, true );
    bottom_iris_l{i} = ind2sub_array( size(boundary_iris{i}), bottom_iris_l{i} );
    bottom_iris_l{i} = [bottom_iris_l{i}; intercept_l{i}];
    
    bottom_iris_r{i} = find_path( points2bw(iris_contour_r{i}, size(boundary_iris{i})), tip_r{i}, intercept_r{i}, true );
    bottom_iris_r{i} = ind2sub_array( size(boundary_iris{i}), bottom_iris_r{i} );
    bottom_iris_r{i} = [bottom_iris_r{i}; intercept_r{i}];
    top_iris_r{i} = find_path(boundary_iris{i}, tip_r{i}, meshwork_r{i});
    top_iris_r{i} = ind2sub_array(size(boundary_iris{i}), top_iris_r{i});
    temp = top_iris_r{i}(:,2);
    [temp,ind] = sort(temp);
    top_iris_r{i} = top_iris_r{i}(ind,:);
    
    to_delete = (temp(1:end-1) == temp(2:end));
    to_delete = to_delete + [0; to_delete(1:end-1)];
    top_iris_r{i}(logical(to_delete),:) = [];
    
    
    top_iris_l{i} = find_path(boundary_iris{i}, tip_l{i}, meshwork_l{i});
    top_iris_l{i} = ind2sub_array(size(boundary_iris{i}), top_iris_l{i});
    temp = top_iris_l{i}(:,2);
    [temp,ind] = sort(temp);
    top_iris_l{i} = top_iris_l{i}(ind,:);
    
    to_delete = (temp(1:end-1) == temp(2:end));
    to_delete = to_delete + [0; to_delete(1:end-1)];
    top_iris_l{i}(logical(to_delete),:) = [];
    
    % fit the lower iris into a polynomial
    LS_A_l = zeros( size(bottom_iris_l{i},1), LOWER_IRIS_FIT_ORDER+1 );
    for mm = 1:size(bottom_iris_l{i},1)
        for nn = 1:LOWER_IRIS_FIT_ORDER+1
            LS_A_l(mm,nn) = bottom_iris_l{i}(mm,2).^(LOWER_IRIS_FIT_ORDER + 1-nn);
        end
    end
    
    LS_l = LS_Problem(LS_A_l, bottom_iris_l{i}(:,1));
    
    temp = zeros(1,LOWER_IRIS_FIT_ORDER+1);
    for mm = 1:LOWER_IRIS_FIT_ORDER+1
        temp(mm) = tip_l{i}(2)^ (LOWER_IRIS_FIT_ORDER+1 - mm);
    end
    
    LS_l = LS_l.add_constraint( temp, tip_l{i}(1) );
    
    poly_coeff = solve(LS_l);
    temp = (min(bottom_iris_l{i}(:,2)) : max(bottom_iris_l{i}(:,2)) )';
    bottom_iris_l_smooth{i} = [ polyval(poly_coeff, temp),temp ];
    
    
    % same thing for the right lower iris
    LS_A_r = zeros( size(bottom_iris_r{i},1), LOWER_IRIS_FIT_ORDER+1 );
    for mm = 1:size(bottom_iris_r{i},1)
        for nn = 1:LOWER_IRIS_FIT_ORDER+1
            LS_A_r(mm,nn) = bottom_iris_r{i}(mm,2).^(LOWER_IRIS_FIT_ORDER + 1-nn);
        end
    end
    
    LS_r = LS_Problem(LS_A_r, bottom_iris_r{i}(:,1));
    
    temp = zeros(1,LOWER_IRIS_FIT_ORDER+1);
    for mm = 1:LOWER_IRIS_FIT_ORDER+1
        temp(mm) = tip_r{i}(2)^ (LOWER_IRIS_FIT_ORDER+1 - mm);
    end
    
    LS_r = LS_r.add_constraint( temp, tip_r{i}(1) );
    
    poly_coeff = solve(LS_r);
    temp = (min(bottom_iris_r{i}(:,2)) : max(bottom_iris_r{i}(:,2)) )';
    bottom_iris_r_smooth{i} = [ polyval(poly_coeff, temp),temp ];
    
    [~,ind] = sort( bottom_iris_r_smooth{i}(:,2));
    bottom_iris_r_smooth{i} = bottom_iris_r_smooth{i}(ind,:);
    [~,ind] = sort( bottom_iris_l_smooth{i}(:,2));
    bottom_iris_l_smooth{i} = bottom_iris_l_smooth{i}(ind,:);
    
    clear mm nn temp

    iris_contour_r{i} = [top_iris_r{i}; flipud(bottom_iris_r_smooth{i})];
    iris_contour_l{i} = [flipud(top_iris_l{i}); bottom_iris_l_smooth{i}];
    
    mid_point(i) = 0.5*(meshwork_l{i}(2) + meshwork_r{i}(2));
%     left_2D{i} = iris_contour_l{i};
%     left_2D{i}(:,2) = -iris_contour_l{i}(:,2) + mid_point(i);
%     right_2D{i} = iris_contour_r{i};
%     right_2D{i}(:,2) = -iris_contour_r{i}(:,2) + mid_point(i);
%     
%     left_3D{i} = zeros(length(left_2D{i}),3);
%     left_3D{i}(:,1) = iris_contour_l{i}(:,2) - mid_point(i);
%     left_3D{i}(:,3) = -iris_contour_l{i}(:,1);
%     right_3D{i} = zeros(length(right_2D{i}),3);
%     right_3D{i}(:,1) = iris_contour_r{i}(:,2) - mid_point(i);
%     right_3D{i}(:,3) = -iris_contour_r{i}(:,1);
%     
%     A_{i} = [meshwork_l{i}(2) - mid_point(i),0,meshwork_l{i}(1)];
%     B{i} = [meshwork_r{i}(2) - mid_point(i),0,meshwork_r{i}(1)];
%     D{i} = (A_{i}+B{i})/2 + [0, norm(A_{i}-B{i},2)/2,0];
%     D_{i} = (A_{i}+B{i})/2 - [0, norm(A_{i}-B{i},2)/2,0];
    
    
    
%     iris_contour_l{i} = nearest(union(union(setdiff(iris_contour_l{i}, bottom_iris_l{i}, 'rows'), bottom_iris_l_smooth{i}, 'rows'), tip_l{i}, 'rows'));
%     iris_contour_r{i} = nearest(union(union(setdiff(iris_contour_r{i}, bottom_iris_r{i}, 'rows'), bottom_iris_r_smooth{i}, 'rows'), tip_r{i}, 'rows'));
        
    % get a denser path for the polynomial (so that the line is connected),
    if ~ismember(i, incorrect_path)
        extended_cornea_dense{i} = [extended_cornea{i}; 1:size(boundary_iris{i},2)]';
        extended_cornea_dense{i} = round(extended_cornea_dense{i});
        ind = find(extended_cornea_dense{i}(:,1) >= 1 & extended_cornea_dense{i}(:,1) <= size(boundary_iris{i},1));
        extended_cornea_dense{i} = extended_cornea_dense{i}(ind,:);
        extended_cornea_dense{i} = round(extended_cornea_dense{i});
        y_min = ceil(min(extended_cornea_dense{i}(:,1))) + 3; y_max = floor(max(extended_cornea_dense{i}(:,1)));
        for y = y_min:1:y_max
            r = roots([fit_coeff{i}(1:LOWER_CORNEA_FIT_ORDER), fit_coeff{i}(LOWER_CORNEA_FIT_ORDER+1) - y]);
            r = r(abs(imag(r)) < 1e-4);
            r = round(real(r));
            r = r( r>=1 & r<= size(boundary_iris{i}, 2));
            assert(length(r) == 2);
            extended_cornea_dense{i} = [extended_cornea_dense{i}; y, r(1); y, r(2)];
        end
    end
    
    
%     % fix the iris shape by cutting it by extended_cornea.
%     
%     
%     % first get the filled shape
%     bw_l{i} = points2bw( [round(iris_contour_l{i}); extended_cornea_dense{i} ], size(original{i}));
%     bw_l{i} = bwmorph(bwmorph(bw_l{i}, 'dilate'), 'thin' );
%     bw_l{i} = bwmorph(imfill(bw_l{i}, 'holes'), 'open');
%     pts = bw2points(bw_l{i}); ind = [];
%     
%     
%     % then cut by extended_cornea
%     for w = 1:size(pts, 1)
%         if pts(w, 1) < extended_cornea{i}(pts(w,2)) + 2
%             ind = [ind, w];
%         end
%     end
%     pts(ind, :) = [];
%     bw_l{i} = points2bw(pts, size(original{i}));
%     bw_l{i} = bwmorph( bwmorph(bw_l{i}, 'clean'), 'remove' );
%     iris_contour_l{i} = bw2points(bw_l{i});
%     
%     
%     % first get the filled shape
%     bw_r{i} = points2bw( [round(iris_contour_r{i}); extended_cornea_dense{i} ], size(original{i}));
%     bw_r{i} = bwmorph(bwmorph(bw_r{i}, 'dilate'), 'thin' );
%     bw_r{i} = bwmorph(imfill(bw_r{i}, 'holes'), 'open');
%     pts = bw2points(bw_r{i}); ind = [];
%     
%     % then cut by extended_cornea
%     for w = 1:size(pts, 1)
%         if pts(w, 1) < extended_cornea{i}(pts(w,2)) + 2
%             ind = [ind, w];
%         end
%     end
%     pts(ind, :) = [];
%     bw_r{i} = points2bw(pts, size(original{i}));
%     bw_r{i} = bwmorph( bwmorph(bw_r{i}, 'clean'), 'remove' );
%     iris_contour_r{i} = bw2points(bw_r{i});
%     
%     iris_contour{i} = [iris_contour_l{i}; iris_contour_r{i}];
    
    
    %     % all the plots!
%     figure, imshow(uint8(original{i})); hold on;
%     scatter(iris_contour_l{i}(:,2), iris_contour_l{i}(:,1), '.g');
%     scatter(iris_contour_r{i}(:,2), iris_contour_r{i}(:,1), '.r');
%     scatter(meshwork_l{i}(2), meshwork_l{i}(1), 'ro', 'linew', 8);
%     scatter(meshwork_r{i}(2), meshwork_r{i}(1), 'ro', 'linew', 8);
%     scatter(intercept_l{i}(2), intercept_l{i}(1), 'or', 'linew', 8);
%     scatter(intercept_r{i}(2), intercept_r{i}(1), 'or', 'linew', 8);
%     
%     
%     scatter(extended_cornea_dense{i}(:,2), extended_cornea_dense{i}(:,1), '.m');
%     scatter(tip_l{i}(:,2), tip_l{i}(:,1), 'ob', 'linew', 8);
%     scatter(tip_r{i}(:,2), tip_r{i}(:,1), 'ob', 'linew', 8);
    
    bottom_iris_r_smooth{i}(bottom_iris_r_smooth{i}(:,2) > meshwork_r{i}(2),:) = [];
    bottom_iris_l_smooth{i}(bottom_iris_l_smooth{i}(:,2) < meshwork_l{i}(2),:) = [];
    
    iris_contour_r{i} = [top_iris_r{i}; flipud(bottom_iris_r_smooth{i})];
    iris_contour_l{i} = [flipud(top_iris_l{i}); bottom_iris_l_smooth{i}];
end

