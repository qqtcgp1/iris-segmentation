figure, imshow( 255 - uint8(original{1}) ); hold on

origin = 0.5*(meshwork_l{1} + meshwork_r{1});
plot( origin(2), origin(1), 'or', 'linew', 7 );

h_text = uicontrol('style', 'text', 'string', 'mark the center', ...
     'unit', 'normalized', 'position', [0.1 0.1 0.2 0.03]);
 
 [center(1),center(2)] = ginput(1);
 plot( center(1), center(2), 'or');
 
 
 i = 1;
 scatter(iris_contour_l{i}(:,2), iris_contour_l{i}(:,1), '.g');
 scatter(iris_contour_r{i}(:,2), iris_contour_r{i}(:,1), '.r');
 scatter(meshwork_l{i}(2), meshwork_l{i}(1), 'ro', 'linew', 8);
 scatter(meshwork_r{i}(2), meshwork_r{i}(1), 'ro', 'linew', 8);
 scatter(intercept_l{i}(2), intercept_l{i}(1), 'or', 'linew', 8);
 scatter(intercept_r{i}(2), intercept_r{i}(1), 'or', 'linew', 8);
 
 
 scatter(extended_cornea_dense{i}(:,2), extended_cornea_dense{i}(:,1), '.m');
 scatter(tip_l{i}(:,2), tip_l{i}(:,1), 'ob', 'linew', 8);
 scatter(tip_r{i}(:,2), tip_r{i}(:,1), 'ob', 'linew', 8);
 
 set(h_text, 'string', 'mark points on the left lens bounary');
 [leftx, lefty] = ginput;
 scatter( leftx, lefty, 'bx');
 
  set(h_text, 'string', 'mark points on the right lens bounary');
 [rightx, righty] = ginput;
 scatter( rightx, righty, 'bx');
 
 leftx = [center(1); leftx]; lefty = [center(2); lefty];
 rightx = [center(1); rightx]; righty = [center(2); righty];
 
 pleft = polyfit(leftx, lefty, 3);
 pright = polyfit(rightx, righty, 3);
 
 leftxplot = min(leftx) - 100:0.1:max(leftx);
 leftyplot = polyval(pleft, leftxplot);
 rightxplot = min(rightx):0.1:max(rightx)+100;
 rightyplot = polyval(pright, rightxplot);
 
 plot(leftxplot, leftyplot, 'r');
 plot(rightxplot, rightyplot, 'r');