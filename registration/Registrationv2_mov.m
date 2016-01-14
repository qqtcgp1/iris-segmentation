clc, clear;

[optimizer, metric] = imregconfig('multimodal');
% optimizer.GradientMagnitudeTolerance = 1.0e-4;
optimizer.GrowthFactor = 1.05;
optimizer.Epsilon = 1e-7;
% optimizer.MinimumStepLength = 1.0e-5;
% optimizer.MaximumStepLength = 6.25e-2;
optimizer.InitialRadius = 6.25e-4;
optimizer.MaximumIterations = 100;
% optimizer.RelaxationFactor = 5.0e-1;
seg_video = VideoReader('seg.avi');
video =  VideoReader ('IR007A.avi');

for i = 1:50
    frame = read(seg_video, i);

    frame = rgb2gray(frame);

    % bframe = imcrop(frame, [170,175,600,300]);
   
    segment{i} = frame;
    
end

for i = 1:50 
    frame = read(video,i);
    
    frame = rgb2gray(frame);
    
    frame = imcrop(frame, [170,175,695,315]);
    
    frame = remove_markings(frame);
    
    original {i} = frame;
end  
  
    data{1} = original{1};

for i = 2:length(original)
    % Transform vectors (rotation and translation)
    temp = segment{i};
    ref = segment{i-1};
  
    tform(i) = imregtform(temp, ref, 'Rigid', optimizer, metric);
    % re-size frames
    
    fprintf([num2str(i) ' 1st completed\n']);
    
    
    data {i} = imwarp(original{i}, tform(i) );
    
    temp = data{i};
    
    [size1, size2] = size(temp);
    
    temp = temp(1+floor(size1/2)-158: floor(size1/2)+158,  1+floor(size2/2)-348:floor(size2/2)+348);
    
    data{i} = temp;
end
    
writerObj = VideoWriter('register.avi');

open(writerObj);

for k = 1:length(original) 
     writeVideo(writerObj,data{k});
end

close(writerObj);