clc, clear;

[optimizer, metric] = imregconfig('multimodal');
% optimizer.GradientMagnitudeTolerance = 1.0e-4;
optimizer.GrowthFactor = 1.05;
optimizer.Epsilon = 1.5e-6;
% optimizer.MinimumStepLength = 1.0e-5;
% optimizer.MaximumStepLength = 6.25e-2;
optimizer.InitialRadius = 6.25e-4;
optimizer.MaximumIterations = 100;
% optimizer.RelaxationFactor = 5.0e-1;
n = 11*15;
video = VideoReader('IR007A.avi');

for i = 1:n
    frame = read(video, i);
    
    frame = imcrop(frame, [170,175,695,315]);
    frame = rgb2gray(frame);
    
    frame = remove_markings(frame);
    original{i} = frame;
end

data{1} = original{1};

ref = data{1};
ref = ref(1:109,:);


for i = 2:length(original)
    temp = original{i};
    ref = original{i-1};
    ref = ref(1:109,:);
    tform(i) = imregtform(temp(1:109,:),ref, 'rigid', optimizer, metric);
    fprintf([num2str(i) ' 1st completed\n']);
    data{i} = imwarp(original{i}, tform(i) );
    temp = data{i};
    [size1, size2] = size(temp);
    
    temp = temp(1+floor(size1/2)-158: floor(size1/2)+158,  1+floor(size2/2)-348:floor(size2/2)+348);
    data{i} = temp;
end
    
writerObj = VideoWriter('peaks.avi');

open(writerObj);

for k = 1:length(original) 
     writeVideo(writerObj,data{k});
end

close(writerObj);