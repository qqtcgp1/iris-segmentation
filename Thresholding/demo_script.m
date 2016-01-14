% Demo to threshold a double or integer image.
% By Image Analyst, November 2010

% Change the current folder to the folder of this m-file.
if(~isdeployed)
	cd(fileparts(which(mfilename)));
end

% Check that user has the Image Processing Toolbox installed.
versionInfo = ver; % Capture their toolboxes in the variable.
hasIPT = false;
for k = 1:length(versionInfo)
	if strcmpi(versionInfo(k).Name, 'Image Processing Toolbox') > 0
		hasIPT = true;
	end
end
if ~hasIPT
	% User does not have the toolbox installed.
	message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
	end
end

% User has the Image Processing Toolbox.
% Continue with the demo.  Do some initialization stuff.
clc;    % Clear command window.
clear;  % Delete all variables.
close all;      % Close all figure windows except those created by imtool.
imtool close all;       % Close all figure windows created by imtool.
workspace;      % Make sure the workspace panel is showing.
fontSize = 20;

% Read in a standard MATLAB gray scale demo image.
folder = 'C:\Program Files\MATLAB\R2010b\toolbox\images\imdemos';
baseFileName = 'cameraman.tif';
fullFileName = fullfile(folder, baseFileName);
if ~exist(fullFileName, 'file')
	% File doesn't exist.  Try it without the folder.
	% It might be able to find it in a folder off the search path.
	fullFileName = baseFileName;
	if ~exist(fullFileName, 'file')
		% Can't find it off the search path either.
		errorMessage = sprintf('Error: cannot find demo image %s', baseFileName);
		uiwait(msgbox(errorMessage));
		return;
	end
end
grayImage = imread(fullFileName);
% Get the dimensions of the image.  numberOfColorBands should be = 1.
[rows columns numberOfColorBands] = size(grayImage);

% Display the original gray scale image.
subplot(2, 3, 1);
imshow(grayImage, []);
axis off;
title('Original Grayscale Image', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
set(gcf, 'name', 'Thresholding Demo by ImageAnalyst', 'numbertitle', 'off') 

message = sprintf('Thresholding demo by ImageAnalyst.\n\nDo you want to use an integer image or a floating point image?');
button = questdlg(message, 'Image Type?', 'Integer', 'Floating Point', 'Cancel', 'Integer');
drawnow;	% Refresh screen to get rid of dialog box remnants.
if strcmpi(button, 'Cancel')
	close(gcf);	% Get rid of window.
	return;
end
if strcmpi(button, 'Floating Point')
	% Convert to double in the range -5000 to + 15000
	% Get input min and max.
	minGL = double(min(grayImage(:)));
	maxGL = double(max(grayImage(:)));
	% Scale the image
	imageToThreshold = 20000 * (double(grayImage) - minGL) / (maxGL - minGL) - 5000;
	% Verify them
	minDblGL = min(imageToThreshold(:));
	maxDblGL = max(imageToThreshold(:));
	fprintf('Before scaling, min gray level = %.1f, max gray level = %.1f\nAfter scaling,  min gray level = %.1f, max gray level = %.1f\n', ...
		minGL, maxGL, minDblGL, maxDblGL);
	startingLowThreshold = -1000;
	startingHighThreshold = -10000;
else
	% Integer image.  Just leave it alone.
	imageToThreshold = grayImage;
	startingLowThreshold = 7;
	startingHighThreshold = 23;
	
	% Let's compute and display the histogram, just for fun.
	[pixelCount grayLevels] = imhist(grayImage);
	subplot(2, 3, 2); 
	bar(pixelCount);
	title('Histogram of Original Image', 'FontSize', fontSize);
	xlim([0 grayLevels(end)]); % Scale x axis manually.
	grid on;
end

%====================== KEY PART RIGHT HERE!!!! ===================================================
% Threshold with starting range startingLowThreshold to startingHighThreshold.
[lowThreshold highThreshold] = threshold(startingLowThreshold, startingHighThreshold, imageToThreshold);
%====================== KEY PART RIGHT HERE!!!! ===================================================


% Binarize the image.
binaryImage = (imageToThreshold > lowThreshold) & (imageToThreshold < highThreshold);
subplot(2, 3, 3);
imshow(binaryImage, []);
axis off;
title('Binarized Image, White Pixels = Mask Pixels', 'FontSize', fontSize);

% Compute max and min of the original image.
minValue = min(imageToThreshold(:));
maxValue = max(imageToThreshold(:));

% Make the image inside the mask have a value of zero.
maskedImage = imageToThreshold;
maskedImage(binaryImage) = 0;
subplot(4, 3, 7);
imshow(maskedImage, []);
axis off;
title('Zero Value Inside the Mask', 'FontSize', fontSize);

% Make the image inside the mask have the min value.
maskedImage = imageToThreshold;
maskedImage(binaryImage) = minValue;
subplot(4, 3, 8);
imshow(maskedImage, []);
axis off;
title('Min Value Inside the Mask', 'FontSize', fontSize);

% Make the image inside the mask have the max value.
maskedImage = imageToThreshold;
maskedImage(binaryImage) = maxValue;
subplot(4, 3, 9);
imshow(maskedImage, []);
axis off;
title('Max Value Inside the Mask', 'FontSize', fontSize);

% Now do the same thing but OUTSIDE the mask.
outsideMask = ~binaryImage;

% Make the image outside the mask have a value of zero.
maskedImage = imageToThreshold;
maskedImage(outsideMask) = 0;
subplot(4, 3, 10);
imshow(maskedImage, []);
axis off;
title('Zero Value Outside the Mask', 'FontSize', fontSize);

% Make the image outside the mask have the min value.
maskedImage = imageToThreshold;
maskedImage(outsideMask) = minValue;
subplot(4, 3, 11);
imshow(maskedImage, []);
axis off;
title('Min Value Outside the Mask', 'FontSize', fontSize);

% Make the image outside the mask have the max value.
maskedImage = imageToThreshold;
maskedImage(outsideMask) = maxValue;
subplot(4, 3, 12);
imshow(maskedImage, []);
axis off;
title('Max Value Outside the Mask', 'FontSize', fontSize);

% Alert user we're done.
uiwait(msgbox('Done with demo.'));

