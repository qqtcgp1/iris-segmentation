function [ F ] = superimpose( log_filename, vid_filename, vid_frames, node_numbers, M,C )
%superimpose oct video and FE outcome in the same movie
%   oct video has the name of vid_filename, and frames to be delt with is
%   vid_frames. FE outcome comes from log_filename


filestring = fileread(log_filename);

vid = VideoReader(vid_filename);

%use strsplit to split the text into segments, each segment is a time step,
%except for the first segment.
steps = strsplit(filestring, 'x;y;z'); figure;

imshow( read(vid,1));
F(1) = getframe;

for i = 1:length(vid_frames)
    imshow( read(vid, vid_frames(i)) );
    hold on;
    
    points = zeros( length(node_numbers),3);

    for j = 1:length(node_numbers)
        cursor = strfind(steps{i+1}, [char(10), num2str(node_numbers(j)), ' ']) - 1;
        cursor = cursor + length([char(10), num2str(node_numbers(j)), ' ']);
        A = sscanf(steps{i+1}(cursor:end), '%f',3);

        points(j,1) = A(1) - C(1,1); points(j,2) = A(2) - C(1,2); points(j,3) = A(3) - C(1,3);
    end
    
    X = (points) * inv(M);
    
    for j = 1:length(node_numbers) - 1
        line( X(j:j+1,1), X(j:j+1,2), 'color', 'r');hold on;
    end
    line( X([end 1],1), X([end 1],2), 'color', 'r');
    
    F(i+1) = getframe; hold off;
end

end

