clear
%Read Imgae 153 and 0
image =imread('0.jpg');

%crop the images
image = imcrop(image, [170,175,695,315]);
grayimage= rgb2gray(image);

[rows, columns, numberOfColorBands] = size(grayimage)

if numberOfColorBands > 1
	% It's not really gray scale like we expected - it's color.
	% Convert it to gray scale by taking only the green channel.
	image = grayImage(:, :, 2); % Take green channel.
end 

% Take the Fourier Transform.
F = fft2(grayimage);
realF = log(real(F));

% Sum up the values.
magImage = abs(F).^2;
energy = sum(magImage(:));

threshold = energy*5;