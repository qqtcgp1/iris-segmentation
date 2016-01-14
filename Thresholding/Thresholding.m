frame = imread(['C:/Users/Administrator/Documents/FYP/BEIYAN WANG_OS01/00', num2str(41+1), '_dewarp.tif']);

% 4.turn to gray scale
frame = rgb2gray(frame);
% 5.remove the mark on top-right
frame = remove_markings_Dewarp(frame, 3);


figure (3); imshow (frame); 
grayImage = frame; 

% Get the dimensions of the image.  numberOfColorBands should be = 1.
[rows columns numberOfColorBands] = size(grayImage);

% Display the original gray scale image.
% subplot(2, 3, 1);
figure(1); imshow(grayImage);
% axis off;
% title('Original Grayscale Image', 'FontSize', fontSize);
% set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
% set(gcf, 'name', 'Thresholding Demo by ImageAnalyst', 'numbertitle', 'off') 

% message = sprintf('Thresholding demo by ImageAnalyst.\n\nDo you want to use an integer image or a floating point image?');
% button = questdlg(message, 'Image Type?', 'Integer', 'Floating Point', 'Cancel', 'Integer');
% drawnow;	% Refresh screen to get rid of dialog box remnants.
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
	
% 	% Let's compute and display the histogram, just for fun.
% 	[pixelCount grayLevels] = imhist(grayImage);
% 	subplot(2, 3, 2); 
% 	bar(pixelCount);
% 	title('Histogram of Original Image', 'FontSize', fontSize);
% 	xlim([0 grayLevels(end)]); % Scale x axis manually.
% 	grid on;
end

%====================== KEY PART RIGHT HERE!!!! ===================================================
% Threshold with starting range startingLowThreshold to startingHighThreshold.
[lowThreshold highThreshold] = threshold(startingLowThreshold, startingHighThreshold, imageToThreshold);
%====================== KEY PART RIGHT HERE!!!! ===================================================


% Binarize the image.
binaryImage = (imageToThreshold > lowThreshold) & (imageToThreshold < highThreshold);
figure (2); imshow(binaryImage);
% axis off;
% title('Binarized Image, White Pixels = Mask Pixels', 'FontSize', fontSize);
% 
% % Compute max and min of the original image.
% minValue = min(imageToThreshold(:));
% maxValue = max(imageToThreshold(:));
% 
% % Make the image inside the mask have a value of zero.
% maskedImage = imageToThreshold;
% maskedImage(binaryImage) = 0;
% subplot(4, 3, 7);
% imshow(maskedImage, []);
% axis off;
% title('Zero Value Inside the Mask', 'FontSize', fontSize);
% 
% % Make the image inside the mask have the min value.
% maskedImage = imageToThreshold;
% maskedImage(binaryImage) = minValue;
% subplot(4, 3, 8);
% imshow(maskedImage, []);
% axis off;
% title('Min Value Inside the Mask', 'FontSize', fontSize);
% 
% % Make the image inside the mask have the max value.
% maskedImage = imageToThreshold;
% maskedImage(binaryImage) = maxValue;
% subplot(4, 3, 9);
% imshow(maskedImage, []);
% axis off;
% title('Max Value Inside the Mask', 'FontSize', fontSize);
% 
% % Now do the same thing but OUTSIDE the mask.
% outsideMask = ~binaryImage;
% 
% % Make the image outside the mask have a value of zero.
% maskedImage = imageToThreshold;
% maskedImage(outsideMask) = 0;
% subplot(4, 3, 10);
% imshow(maskedImage, []);
% axis off;
% title('Zero Value Outside the Mask', 'FontSize', fontSize);
% 
% % Make the image outside the mask have the min value.
% maskedImage = imageToThreshold;
% maskedImage(outsideMask) = minValue;
% subplot(4, 3, 11);
% imshow(maskedImage, []);
% axis off;
% title('Min Value Outside the Mask', 'FontSize', fontSize);
% 
% % Make the image outside the mask have the max value.
% maskedImage = imageToThreshold;
% maskedImage(outsideMask) = maxValue;
% subplot(4, 3, 12);
% imshow(maskedImage, []);
% axis off;
% title('Max Value Outside the Mask', 'FontSize', fontSize);



