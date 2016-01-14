frameg = imread('frame.png');

% frameg=rgb2gray(frame);



%median filter
frameg_median =medfilt2(frameg,[8 8]);

%Gaussian Filter
h = fspecial('gaussian',[8 8],0.5);
frameg_gaussian =imfilter(frameg_median,h);

%gray to bw image
%F = graythresh(frameg);
%framebw=im2bw(frameg_gaussian,F);

%averaging filter
 H = FSPECIAL('average',[7 7]);
frameg_ave = imfilter(frameg_gaussian,H);


I = (frameg_ave < 16);


%second average filtering
H = FSPECIAL('average',[4 4]);
frameg_ave2 = imfilter(BW,H);


%second median filtering
frameg_median2 =medfilt2( frameg_ave2,[9 9]);


%segmentation
a =fspecial('sobel');
A=imfilter(frameg_median2,a);

B = edge( frameg_median2,'sobel',[],'both');
figure (1);imshow(B); title('Sobel');


C = edge(frameg_median2,'canny',[],1);

figure (3);imshow(C); title('Canny')

figure (2);imshow(frameg);title('original');