clear;

%
video =  VideoReader ('IR007A.avi');
% seg_video = VideoReader('seg.avi');

n = 10; 

for i = 1 : n
    
       % 2.Read the first frame
    frame = read(video,i);
    
    %convert frame to grayscale
    frame = rgb2gray(frame);
    
    frame = imcrop(frame, [170,175,695,315]);
    
    frame = remove_markings (frame);
    %median filter
    
    original {i} = frame; 
    
    
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
    
    %segmentation by Sobel
    frame = ~frame;
    
    frame (120:end,:) = 0;
    
    seg = edge(frame,'canny');
    
    boundary {i} = seg;
       
    fprintf([num2str(i), 'th completed.\n']);
    
end

fixedPoints = bw2points (boundary{1}); 

for i = 1: n
    
    point = bw2points( boundary{i}); 

% temp1 = boundary {i} ; 
%find point by rows.. 
% 
% point = [nan nan];
%     
%     for k = 1: size(temp1, 2)
%         for j = 1:size (temp1,1)
%             if temp1 (j,k)
%                 point = [point; j,k ];
%             end
%         end
%     end
%     
%     point(1,:) = [];
%     edge {i} = point;
    
    for j = 1 : size (point, 1)
        if point (j, 2) > 150
            point (j, 1) = point (j, 1);
            point (j,2) = point (j, 2);
        end
    end
    






    tform (i) = fitgeotrans (point,fixedPoints, 'Similarity');
           
    data {i} = imwarp(original{i}, tform(i) );
    
    temp = data{i};
    
    [size1, size2] = size(temp);
    
    temp = temp(1+floor(size1/2)-158: floor(size1/2)+158,  1+floor(size2/2)-348:floor(size2/2)+348);
    
    data{i} = temp;
    
    fprintf([num2str(i), 'th completed.\n']);
    
        point(:,:) = 0;
        
end

writerObj = VideoWriter('register_latest.avi');

open(writerObj);

for k = 1:length(original) 
    
     writeVideo(writerObj,data{k});
     
end

close(writerObj);
 
    




% for i = 1:150
%     
%  
%     tform (i) = fitgeotrans (edge{i}, edge{1}, 'Similarity');
%            
%     data {i} = imwarp(original{i}, tform(i) );
%     
%     fprintf([num2str(i), 'th completed.\n']);
% end

writerObj = VideoWriter('register_latest.avi');

open(writerObj);

for k = 1:length(original) 
     writeVideo(writerObj,data{k});
end

close(writerObj);

scatter (point(:,2),point(:,1),'r.')'
set(gca, 'ydir', 'rev');
equal axis
 
axis equal

% video =  VideoReader ('IR007A.avi');
%
% for i = 1:50
%     frame = read(video,i);
%
%     frame = rgb2gray(frame);
%
%     frame = imcrop(frame,[170,175,695,315]);
%
%     original {i} = frame;
% end
%
% tform = (



% for i = 1:50
%     frame = read(video,i);
%
%     frame = rgb2gray(frame);
%
%     frame = imcrop(frame,[170,175,695,315]);
%
%     frame = remove_markings (frame);
%
%     original {i} = frame;
% end
%
% cpselect(original{2},original{1} );
%

