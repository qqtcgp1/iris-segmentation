clear 

video = VideoReader('IR007A.avi');

for i = 1 : 11*15
% 2.Read the first frame
frame = read(video,i);

%convert frame to grayscale
frame = rgb2gray(frame);

frame = imcrop(frame, [170,175,695,315]);

frame = remove_markings (frame);
%median filter
frame = medfilt2(frame,[8 8]);

%Gaussian Filter
h = fspecial('gaussian',[8 8],0.5);
frame =imfilter(frame,h);

%Averaging filter

H = fspecial('average',[7 7]);

frame = imfilter(frame,H);


%Thresholding

frame = (frame < 10);

%second average filtering
H = fspecial('average',[4 4]);

frame = imfilter(frame,H);

%second median filtering
frame =medfilt2(frame,[9 9]);

%%

%segmentation by Sobel and Canny

seg_line = edge(frame,'sobel',[],'both');

seg_line (110:end,:) = 0;

boundary{i} = seg_line;

fprintf([num2str(i), 'th completed.\n']); 
end

writerObj = VideoWriter('seg.avi');

open(writerObj);

for i = 1:11*15
  
     writeVideo(writerObj,double(boundary{i}));
end

close(writerObj);

