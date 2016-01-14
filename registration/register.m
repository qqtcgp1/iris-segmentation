clear

% 1.Read Video
video = VideoReader('IR004B.avi');
video_registered = video;
n = get(video, 'NumberOfFrames');

[optimizer, metric] = imregconfig('multimodal');
% optimizer.GradientMagnitudeTolerance = 1.0e-4;
optimizer.GrowthFactor = 1.05;
optimizer.Epsilon = 1.5e-6;
% optimizer.MinimumStepLength = 1.0e-5;
% optimizer.MaximumStepLength = 6.25e-2;
optimizer.InitialRadius = 6.25e-3;
optimizer.MaximumIterations = 100;
% optimizer.RelaxationFactor = 5.0e-1;


for i = 1:n
    frame = read(video, i);
    
    frame = imcrop(frame, [170,175,695,315]);
    frame = rgb2gray(frame);
    
    frame(50:72,603:638) = 0;
    original{i} = frame;
end

original([82,83]) = [];
data{1} = original{1};
figure, imshow(data{1}); M(1) = getframe;

for i = 2:length(original)
    
    data{i} = imregister(original{i}, original{1}, 'rigid', optimizer, metric);
    
    imshow(data{i});
    M(i) = getframe;
    fprintf([num2str(i) ' 1st completed\n']);
end

% data2{1} = data{1}; imshow(data2{1}); N(1) = getframe;
% 
% for i = 2:length(original)
%     data2{i} = imregister(data{i}, data{i-1}, 'translation', optimizer, metric);
%     
%     imshow(data2{i}); N(i) = getframe;
%     fprintf([num2str(i) ' 2nd completed\n']);
% end
    