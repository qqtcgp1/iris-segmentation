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


for i = 1:51
    % 2.Read the frame
    frame = imread(['/Users/yanghua/Documents/FYP/Iris/Dewarp/BEIYAN WANG_OS01/00', num2str(41+i), '_dewarp.tif']);
    % 4.turn to gray scale
    frame = rgb2gray(frame);
    % 5.remove the mark on top-right
    frame = remove_markings_Dewarp(frame,3);
    original{i} = frame;
end


for i = 2:length(original)
    data{i} = imregister(original{i}, original{1}, 'rigid', optimizer, metric);
    fprintf([num2str(i) ' 1st completed\n']);
end


writerObj = VideoWriter('register.avi');

open(writerObj);

for k = 1:length(original)
    writeVideo(writerObj,double(data{k}));
end

close(writerObj);