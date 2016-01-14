clear

frame = imread('0.jpg');

frameg = rgb2gray(frame);

frameg = imcrop(frameg,[170,175,695,315]);

frameg = remove_markings (frameg);


% %median filter

% frameg_median =medfilt2(frameg,[3 3]);
% 
% %Gaussian Filter
% h = fspecial('gaussian',[3 3],0.5);
% frameg_gaussian =imfilter(frameg_median,h);

% %gray to bw image
% %F = graythresh(frameg);
% %framebw=im2bw(frameg_gaussian,F);
% 
% %averaging filter
%  H = fspecial('average',[4 4]);
% frameg_ave = imfilter(frameg_gaussian,H);

image = (frameg > 40);

%remove black spots

image =imfill(image,'holes');

%remove white spots
image = ~imfill(~image,'holes');

se = strel('disk',2);


% Use a small black ball to roll over insdie the white edge..
image = imclose(image,se);
 
image = ~imclose(~image, se);
 
clear se;

% figure; imshow(image);

% %second average filtering
% H = fspecial('average',[3 3]);
% frameg_ave2 = imfilter(I,H);
% 
% 
% %second median filtering
% frameg_median2 =medfilt2( frameg_ave2,[3 3]);
% 
% 
% %segmentation
% a =fspecial('sobel');
% A=imfilter(frameg_median2,a);

image_outline = edge( image,'sobel',[],'both');

% 
% figure, imshow(frame);
% hold on;
% edge = [nan nan];
% for i=1:size(image,1)
%     for j=1:size(image,2)
%         if image(i,j)
%             edge = [edge; i,j];
%         end
%     end
% end

 points = bw2points (image_outline); 

figure (1); imshow (frameg);
hold on;
% imshow(image_outline); title('Sobel');
% hold on; 
scatter(points(:,2),points(:,1),'g.');



%
% C = edge(frameg_median2,'canny',[],1);
% 
% figure (3);imshow(C); title('Canny')
% 
% figure (2);imshow(frameg);title('original');


