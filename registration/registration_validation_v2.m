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
    
    frame = imcrop(frame, [50 1 500 120]);
    
    cornea {i} = frame; 
end

data{1} = original{1};


for i = 2:length(original)
    tform(i) = imregtform(cornea {i},cornea{1}, 'rigid', optimizer, metric);
    
    fprintf([num2str(i) ' 1st completed\n']);
    
    data{i} = imwarp(original{i}, tform(i) );
%     temp = data{i};
%     [size1, size2] = size(temp);
%     
%     temp = temp(1+floor(size1/2)-158: floor(size1/2)+158,  1+floor(size2/2)-348:floor(size2/2)+348);
%     data{i} = temp;
end



% for i = 2:length(original)
%     data{i} = imregister(cornea{i}, cornea{1}, 'rigid', optimizer, metric);
%     
%     fprintf([num2str(i) ' 1st completed\n']);
% end

% data {1} = original {1}; 
% 
% writerObj = VideoWriter('register1.avi');
% 
% open(writerObj);
% 
% for k = 1:length(original)
%     writeVideo(writerObj,data{k});
% end

close(writerObj);
calculate the difference btw register and unregister MI

 for i = 1: 34
      temp = data{i};
      temp = temp (1:120, :); 
      
    ref = original{1};
    ref = ref(1:120,:);
register (i) = MI (temp, ref);
unreg = original {i}; 
unreg = unreg(1:120, :); 
unregister (i) = MI(unreg, ref);
difference (i) = register(i) - unregister(i);
 end

 
 
 
 
 