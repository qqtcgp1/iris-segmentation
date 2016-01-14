

%%
%Load image
clear
vid = VideoReader('IR007A.avi');

frameg = rgb2gray(read(vid, 1));

%%Image must be in grayscale!
%median filter
frameg_median =medfilt2(frameg,[3 3]);

imshow(frameg_median)
%Gaussian Filter
h = fspecial('gaussian',[3 3],0.5);
frameg_gaussian =imfilter(frameg_median,h);

% % gray to bw image
% F = graythresh(frameg);
% framebw=im2bw(frameg_gaussian,F);

%averaging filter
H = fspecial('average',[3 3]);
frameg_ave = imfilter(frameg_gaussian,H);

%%
%Thresholding
I = (frameg_ave < 5);
%%

%second average filtering
H = fspecial('average',[4 4]);
frameg_ave2 = imfilter(I,H);


%second median filtering
frameg_median2 =medfilt2(frameg_ave2,[9 9]);

%%

%segmentation by Sobel and Canny
a =fspecial('sobel');
A = imfilter(frameg_median2,a);

B = edge( frameg_median2,'sobel',[],'both');
%%
%Plot the image
figure (1);imshow(B); title('Sobel');


C = edge(frameg_median2,'canny',[],1);

figure (3);imshow(C); title('Canny')

figure (2);imshow(frameg);title('original');
%%
