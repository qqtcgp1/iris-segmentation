for i = 1 : 250
    currentFrame = imread(int2str(i-1),'jpg');
    currentFrame = rgb2gray(currentFrame);
    currentFrame = imcrop(currentFrame, [170,175,695,315]);
    currentFrame(50:72,603:638) = 0;
    combinedString=strcat(int2str(i-1),'.jpg');
    imwrite(currentFrame,combinedString);
   end