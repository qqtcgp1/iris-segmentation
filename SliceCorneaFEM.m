function SliceCorneaFEM

% Initialization parameters 1 (start)

vid1 = VideoReader('health control.avi'); % Choose video
numx = 120; % Number of x segments 120
numy = 5; % Number of y segments 5
numz = 8000; % Number of y segments 8000
cct = 0.5; % Thickness of central cornea (mm)
density = 1000; % Density of cornea

numberofvar = 2;
materialvar = 1; %Isotropic Elastic
var1 = 2*(10^6); % Young's modulus of cornea
var2 = 0.49; % Poisson's ratio of cornea
var3 = []; % NA 
materialtype = 'isotropic elastic';

% numberofvar = 3;
% materialvar = 2; %Veronda-Westmann
% var1 = 10000;
% var2 = 10000;
% var3 = 10000;
% materialtype = 'Veronda-Westmann';

visco = 1; %1=yes; 0= no;

if materialvar == 1
    viscotype = 'viscoelastic';
else
    viscotype = 'uncoupled viscoelastic';
end 

numberofvisco = 1;

gvar(1) = 0.95; %default=0
gvar(2) = 0;
gvar(3) = 0;
gvar(4) = 0;
gvar(5) = 0;
gvar(6) = 0;
tvar(1) = 0.000231; %default=1
tvar(2) = 1;
tvar(3) = 1;
tvar(4) = 1;
tvar(5) = 1;
tvar(6) = 1;

IOP = 1053; %1053 Intraocular Pressure of cornea
rairjet = 1.53*(10^-3); %1.53*(10^-3) Radius of airjet
airjetp = 8000; %8000 Pressure of airjet

smooth1 = 0.3;

timestepmodel = 138; %138
timestepsize = 0.000231; %0.000231
penalty = 500000;
analysismethod = 'static'; %'dynamic' or 'static'
chosenname = 'trial1'; % Name chosen

loadfunction(1,1) = 0.00230375675675676;
loadfunction(1,2) = 0;

loadfunction(2,1) = 0.00383221;
loadfunction(2,2) = 0;

loadfunction(3,1) = 0.00516940540540541;
loadfunction(3,2) = 0.0230769230769231;

loadfunction(4,1) = 0.00646175675675676;
loadfunction(4,2) = 0.196153846153846;

loadfunction(5,1) = 0.00921502702702703;
loadfunction(5,2) = 0.515384615384615;

loadfunction(6,1) = 0.0119682972972973;
loadfunction(6,2) = 0.796153846153846;

loadfunction(7,1) = 0.0158453513513514;
loadfunction(7,2) = 1;

loadfunction(8,1) = 0.0170815135135135;
loadfunction(8,2) = 0.946153846153846;

loadfunction(9,1) = 0.0188233783783784;
loadfunction(9,2) = 0.7;

loadfunction(10,1) = 0.0215204594594595;
loadfunction(10,2) = 0.230769230769231;

loadfunction(11,1) = 0.0226442432432432;
loadfunction(11,2) = 0.0192307692307692;

loadfunction(12,1) = 0.0232623243243243;
loadfunction(12,2) = -0.0346153846153846;

loadfunction(13,1) = 0.0239365945945946;
loadfunction(13,2) = -0.0615384615384615;

loadfunction(14,1) = 0.0247232432432432;
loadfunction(14,2) = -0.0423076923076923;

loadfunction(15,1) = 0.0251165675675676;
loadfunction(15,2) = -0.00769230769230769;

loadfunction(16,1) = 0.0257908378378378;
loadfunction(16,2) = 0.0269230769230769;

loadfunction(17,1) = 0.026127972972973;
loadfunction(17,2) = 0.0423076923076923;

loadfunction(18,1) = 0.0265212972972973;
loadfunction(18,2) = 0;

% Initialization parameters 1 (end)

currentos = -1; 
if ispc == 1
    currentos = 1;
else
    if isunix ==1
        currentos =0;
    end
end 

numx2 = numx - 1;
numx3 = floor(numx/2)-1;
numy2 = numy - 1; 
numz2 = numz/4; 
numz3 = 1+1;  %1+1
anglez = (pi/(2*numz2))/(numz3-1);
totalnumnodes = floor(numx/2)*numy*numz3+numy;

fname1 = char(strcat(chosenname, '_', 'Slice' )); % File name
fname2 = char(strcat('Original_', fname1, '.txt'));
file1 = fopen(fname2, 'wt');

fnamevid = char(strcat('Video_', fname1, '.txt'));

filegendir1=mfilename('fullpath');
filegendir2=mfilename('path');
filegendir3=size(filegendir1);
filegendir4=size(filegendir2);

for filegencounter1 = 1 : filegendir3(2) - filegendir4(2)
    filegendir5(filegencounter1) = filegendir1(filegencounter1);
end

filelocation1 = char(strcat(filegendir5));

% Read video into matrix 
vid2 = read(vid1); 
vid2 = squeeze(vid2);

if (size(vid2,4)>1)
    [m,n,n2, numIm] = size(vid2);
else
    [m,n,numIm] = size(vid2);
end

timestepvideo = numIm;

fnameinitial = char(strcat('InitialData', '.txt'));
fileinitial = fopen(fnameinitial, 'wt');

fprintf(fileinitial, '%d,\n', numx );
fprintf(fileinitial, '%d,\n', numy );
fprintf(fileinitial, '%d,\n', numz );
fprintf(fileinitial, '%f,\n', cct );
fprintf(fileinitial, '%f,\n', density );
fprintf(fileinitial, '%f,\n', numberofvar );
fprintf(fileinitial, '%f,\n', materialvar );
fprintf(fileinitial, '%s,\n', materialtype );
fprintf(fileinitial, '%f,\n', var1 );
fprintf(fileinitial, '%f,\n', var2 );
fprintf(fileinitial, '%f,\n', var3 );
fprintf(fileinitial, '%f,\n', visco );
fprintf(fileinitial, '%s,\n', viscotype );
fprintf(fileinitial, '%d,\n', numberofvisco );
for viscounter = 1:size(gvar,2)
    fprintf(fileinitial, '%f,\n', gvar(viscounter) );
end
for viscounter = 1:size(tvar,2)
    fprintf(fileinitial, '%f,\n', tvar(viscounter) );
end
fprintf(fileinitial, '%f,\n', IOP );
fprintf(fileinitial, '%f,\n', rairjet );
fprintf(fileinitial, '%f,\n', airjetp );
fprintf(fileinitial, '%f,\n', smooth1 );
fprintf(fileinitial, '%f,\n', timestepmodel );
fprintf(fileinitial, '%f,\n', timestepsize );
fprintf(fileinitial, '%f,\n', timestepvideo );
fprintf(fileinitial, '%s,\n', analysismethod );
fprintf(fileinitial, '%s,\n', chosenname );

fnametimestepchoice = char(strcat('TimeChoice_', fname1, '.txt'));
filetimestepchoice = fopen(fnametimestepchoice, 'wt');

for timestepchoicecount1 = 1 : 138
    fprintf(filetimestepchoice, '%d,\n', timestepchoicecount1*2.31*(10^-4));
end

fnameload = char(strcat('Loadfunction_', fname1, '.txt'));
fileload = fopen(fnameload, 'wt');

for loadcount1 = 1 :size(loadfunction,1)
    fprintf(fileload, '%f, %f\n', loadfunction(loadcount1,1),loadfunction(loadcount1,2));
end

% for timestepchoicecount1 = 1 : 13
%     fprintf(filetimestepchoice, '%d,\n', timestepchoicecount1*2.31*(10^-3));
% end

%for im = 1:1
for im = 1:numIm

    % Read image
    if (size(vid2,4)>1)
        image = (vid2(:,:,1,im)+vid2(:,:,2,im)+vid2(:,:,3,im))/3;
    else
        image = vid2(:,:,im);
    end

    % Crop image to remove wordings on the image
    %image2 = image(40:end-40,15:end-15);
    
    image2 = image;
    
    for delrow1 = 6:16
        for delcol1 = 7:171
            image2(delrow1, delcol1) = 0;
        end
    end
    
    for delrow2 = 6:16
        for delcol2 = 257:321
            image2(delrow2, delcol2) = 0;
        end
    end
    
    for delrow3 = 6:16
        for delcol3 = 385:571
            image2(delrow3, delcol3) = 0;
        end
    end
        
    for delrow4 = 182:200 %195
        for delcol4 = 8:165
            image2(delrow4, delcol4) = 0;
        end
    end
    
    image2(181,22) = 0;
    image2(181,23) = 0;
    image2(180,23) = 0;
    image2(180,24) = 0;
    image2(179,24) = 0;
    image2(178,24) = 0;
    image2(178,25) = 0;
    image2(177,25) = 0;
    image2(177,26) = 0;
    image2(176,26) = 0;
    image2(176,27) = 0;
    image2(175,26) = 0;
    image2(175,27) = 0;
    image2(174,27) = 0;
    image2(173,27) = 0;
    image2(181,28) = 0;
    image2(174,28) = 0;
    image2(173,28) = 0;
    image2(172,28) = 0;
    image2(181,29) = 0;
    image2(175,29) = 0;
    image2(174,29) = 0;
    image2(173,29) = 0;
    image2(172,29) = 0;
    image2(171,29) = 0;
    image2(170,29) = 0;
    image2(181,30) = 0;
    image2(180,30) = 0;
    image2(177,30) = 0;
    image2(176,30) = 0;
    image2(175,30) = 0;
    image2(174,30) = 0;
    image2(173,30) = 0;
    image2(172,30) = 0;
    image2(171,30) = 0;
    image2(170,30) = 0;
    image2(169,30) = 0;
    image2(180,31) = 0;
    image2(179,31) = 0;
    image2(178,31) = 0;
    image2(177,31) = 0;
    image2(169,31) = 0;
    image2(168,31) = 0;
    image2(167,31) = 0;
    image2(179,32) = 0;
    image2(167,32) = 0;
    image2(166,32) = 0;
    image2(180,33) = 0;
    image2(179,33) = 0;
    image2(166,33) = 0;
    image2(165,33) = 0;
    image2(164,33) = 0;
    image2(181,34) = 0;
    image2(180,34) = 0;
    image2(166,34) = 0;
    image2(165,34) = 0;
    image2(164,34) = 0;
    image2(181,35) = 0;
    image2(167,35) = 0;
    image2(166,35) = 0;
    image2(181,36) = 0;
    image2(169,36) = 0;
    image2(168,36) = 0;
    image2(167,36) = 0;
    image2(169,37) = 0;
    image2(170,37) = 0;
    image2(170,38) = 0;
    image2(171,38) = 0;
    image2(172,38) = 0;
    image2(172,39) = 0;
    image2(173,39) = 0;
    image2(173,40) = 0;
    image2(174,40) = 0;
    image2(175,40) = 0;
    image2(175,41) = 0;
    image2(176,41) = 0;
    image2(177,41) = 0;
    image2(177,42) = 0;
    image2(178,42) = 0;
    image2(178,43) = 0;
    image2(179,43) = 0;
    image2(180,43) = 0;
    image2(180,44) = 0;
    image2(181,44) = 0;
    image2(181,45) = 0;
    
    for delrow5 = 187:195
        for delcol5 = 530:555
            image2(delrow5, delcol5) = 0;
        end
    end
    
    for delrow6 = 190:195
        for delcol6 = 560:571
            image2(delrow6,delcol6) = 0;
        end 
    end 

    % Lowpass filter
    H = fspecial('gaussian',3);
    image2 = imfilter(image2,H);

    % Smooth image
    H = fspecial('average',3);
    image2 = imfilter(image2,H);

    % Threshold using Otsu method
    [m,n] = size(image2);
    
    for gg=1:n
        th = graythresh(image2(:,gg));
        bw(:,gg) = im2bw(image2(:,gg),th);
    end

    % Remove disconnected points
    CC = bwconncomp(bw);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    ind = find(numPixels<10);
    bw2=bw;
    
    for ee2 = 1:length(ind)
        bw2(CC.PixelIdxList{ind(ee2)}) = 0;
    end

    % Fill in small gaps with morphological processing
    se = strel('disk',2);
    bw3 = imclose(bw2,se);
    bw3 = bwareaopen(bw3, 10);

    % Obtain boundary curve perimeter
    bw4 = bwperim(bw3);

    % Remove perimeter at the edge
    bw4(1:3,:)=0;
    bw4(end-2:end,:)=0;
    bw4(:,1:3)=0;
    bw4(:,end-2:end)=0;

    % Find the two connected boundaries
    CC = bwconncomp(bw4);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [Y,idx] = sort(numPixels,'descend');
    [r1,c1] = ind2sub(size(bw4),CC.PixelIdxList{idx(1)});
    [r2,c2] = ind2sub(size(bw4),CC.PixelIdxList{idx(2)});
    
    if r1(15)<r2(15)
        upperpoints = [r1,c1];
        lowerpoints= [r2,c2];
    else
        upperpoints = [r2,c2];
        lowerpoints = [r1,c1];
    end
    
%     imshow(image2)
%     hold on;
%     plot(upperpoints(:,2), upperpoints(:,1));
%     plot(lowerpoints(:,2), lowerpoints(:,1));
%     axis equal;
%     hold off;
    
    conpointeru =[];
    conpointerl =[];
    
    concounter = -1;
    
    for con = 1 : (size(upperpoints,1)-1)
        
        if ((upperpoints(con,2) == upperpoints(con+1,2)) && (con == 1))
            concounter = concounter +2;
            conpointeru(concounter) = con;
            conpointeru(concounter+1) = con+1;            
        end
        
        if ((upperpoints(con,2) == upperpoints(con+1,2)) && (con ~= 1) && (upperpoints(con,2) ~= upperpoints(con-1,2)))
            concounter = concounter +2;
            conpointeru(concounter) = con;
            conpointeru(concounter+1) = con+1;            
        end
        
        if ((upperpoints(con,2) == upperpoints(con+1,2)) && (con ~= 1) && (upperpoints(con,2) == upperpoints(con-1,2)))
            conpointeru(concounter+1) = con+1;        
        end
        
    end
    
    concounter = -1;
    
    for con = 1 : (size(lowerpoints,1)-1)
        
        if ((lowerpoints(con,2) == lowerpoints(con+1,2)) && (con == 1))
            concounter = concounter +2;
            conpointerl(concounter) = con;
            conpointerl(concounter+1) = con+1;            
        end
        
        if ((lowerpoints(con,2) == lowerpoints(con+1,2)) && (con ~= 1) && (lowerpoints(con,2) ~= lowerpoints(con-1,2)))
            concounter = concounter +2;
            conpointerl(concounter) = con;
            conpointerl(concounter+1) = con+1;            
        end
        
        if ((lowerpoints(con,2) == lowerpoints(con+1,2)) && (con ~= 1) && (lowerpoints(con,2) == lowerpoints(con-1,2)))
            conpointerl(concounter+1) = con+1;        
        end
        
    end
    
    if ~isempty(conpointeru)
        
        for con3 = 1 : floor(size(conpointeru,2)/2)
            
            con4 = 2*con3-1;
            consum  = 0;
            
            for con5 = conpointeru(con4) : conpointeru(con4+1)
                consum = consum + upperpoints(con5,1);
            end
            
            conave = consum/(upperpoints(conpointeru(con4+1))-upperpoints(conpointeru(con4))+1);
            
            for con5 = conpointeru(con4) : conpointeru(con4+1)
                upperpoints(con5,1)=conave;
            end        
        end
    end
    
    if ~isempty(conpointerl)     
        
        for con3 = 1 : floor(size(conpointerl,2)/2)
            
            con4 = 2*con3-1;
            consum  = 0;
            
            for con5 = conpointerl(con4) : conpointerl(con4+1)
                consum = consum + lowerpoints(con5,1);
            end
            
            conave = consum/(lowerpoints(conpointerl(con4+1))-lowerpoints(conpointerl(con4))+1);
            
            for con5 = conpointerl(con4) : conpointerl(con4+1)
                lowerpoints(con5,1)=conave;
            end
        end
    end 
    
    for con3 = 1 : floor(size(conpointerl,2)/2)
        
        con5 = ceil(size(conpointerl,2)) - (2*con3-1);
        con6 = ((conpointerl(con5)+1) : 1 : (conpointerl(con5+1)));
        
        for con7 = 1:size(con6,2)
            lowerpoints(con6(size(con6,2)+1-con7),:) = [];
        end
        
    end
    
    for con3 = 1 : floor(size(conpointeru,2)/2)
        
        con5 = ceil(size(conpointeru,2)) - (2*con3-1);
        con6 = ((conpointeru(con5)+1) : 1 : (conpointeru(con5+1)));
        
        for con7 = 1:size(con6,2)
            upperpoints(con6(size(con6,2)+1-con7),:) = [];
        end
        
    end

%     upperpointssp(:,1) = spline(upperpoints((1:smooth1:end),2),upperpoints((1:smooth1:end),1), upperpoints(:,2));
%     upperpointssp(:,2) = upperpoints(:,2);
%     upperpoints = upperpointssp; 
%     
%     lowerpointssp(:,1) = spline(lowerpoints((1:smooth1:end),2),lowerpoints((1:smooth1:end),1), lowerpoints(:,2));
%     lowerpointssp(:,2) = lowerpoints(:,2);
%     lowerpoints = lowerpointssp;

    if im ==1
        standardx = upperpoints(:,2);
    end

    upperpointsfit = fit(upperpoints(:,2), upperpoints(:,1), 'smoothingspline', 'SmoothingParam', smooth1);
    lowerpointsfit = fit(lowerpoints(:,2), lowerpoints(:,1), 'smoothingspline', 'SmoothingParam', smooth1);
    
    upperpointssp(:,2) = standardx;
    upperpointssp(:,1) = upperpointsfit(upperpointssp(:,2));
    lowerpointssp(:,2) = standardx;
    lowerpointssp(:,1) = lowerpointsfit(lowerpointssp(:,2));
    
    upperpoints = upperpointssp; 
    lowerpoints = lowerpointssp;    

%     hold on; 
%     figure(1);
%     axis equal
%     plot(upperpointsfit, upperpoints(:,2), upperpoints(:,1));
%     hold off;
    
%     hold on; 
%     figure(1);
%     axis equal
%     plot(upperpointssp(:,2), upperpointssp(:,1));
%     plot(lowerpointssp(:,2), lowerpointssp(:,1));
%     hold off;
    
    if im == 1
        midcounter1(1) = 0;
        midcounter2(1) = 0;
        midcounter1(2) = 0;
        midcounter2(2) = 0;

        for im6 = 1:size(upperpoints,1)
            if max(-upperpoints(:,1)) == -upperpoints(im6,1)
                midcounter1(1) = im6;
            end   
            if max(-lowerpoints(:,1)) == -lowerpoints(im6,1)
                midcounter2(1) = im6;
            end
        end
    
        for im7 = 1:size(upperpoints,1)
            if max(-upperpoints(:,1)) == -upperpoints(size(upperpoints,1)-im7+1,1)
                midcounter1(2) = size(upperpoints,1)-im7+1;
            end   
            if max(-lowerpoints(:,1)) == -lowerpoints(size(upperpoints,1)-im7+1,1)
                midcounter2(2) = size(upperpoints,1)-im7+1;
            end
        end
        
        midcounter1(3) = (midcounter1(1)+midcounter1(2))/2;
        midcounter2(3) = (midcounter2(1)+midcounter2(2))/2;
    
        midcounter3(1)= ceil((midcounter1(3)+midcounter2(3))/2);
    end
    
    %if less than or equals to size(upperpoints,1)/2, then left half 
    
    if im ==1
        editratio = (-upperpoints(midcounter3,1) - -lowerpoints(midcounter3,1))*1000/cct;
    end
    
    if upperpoints(size(upperpoints,1),2)- upperpoints(midcounter3(1),2) >= upperpoints(midcounter3(1),2) - upperpoints(1,2)
        symlength1 = upperpoints(midcounter3(1),2) - upperpoints(1,2);
        xendpoint1 = upperpoints(midcounter3(1),2) + symlength1;
        for im8 = 1:size(upperpoints,1)
            if upperpoints(im8,2)<=xendpoint1
                xendpoint2 = im8;
            end
        end
        x1 = upperpoints(1:xendpoint2,2);
        x2 = lowerpoints(1:xendpoint2,2);
        y1 = -upperpoints(1:xendpoint2,1);
        y2 = -lowerpoints(1:xendpoint2,1);
    else
        symlength1 = upperpoints(end,2)-upperpoints(midcounter3(1),2);
        xendpoint1 = upperpoints(midcounter3(1),2) - symlength1;
        for im8 = 1:size(upperpoints,1)
            if upperpoints(im8,2)<=xendpoint1
                xendpoint2 = im8;
            end
        end
        x1 = upperpoints(xendpoint2:end,2);
        x2 = lowerpoints(xendpoint2:end,2);
        y1 = -upperpoints(xendpoint2:end,1);
        y2 = -lowerpoints(xendpoint2:end,1);
    end
    
    if xendpoint2 > midcounter3(1)
        midcounter4(1) = midcounter3(1);
    else
        midcounter4(1) = midcounter3(1) - (xendpoint2-1);
    end
    
    %hold on;
    %if size(upperpoints,1)/2 >= midcounter3(1)
        %plot(upperpoints(1:xendpoint2,2), -upperpoints(1:xendpoint2,1));
        %plot(lowerpoints(1:xendpoint2,2), -lowerpoints(1:xendpoint2,1));
    %else
        %plot(upperpoints(xendpoint2:end,2), -upperpoints(xendpoint2:end,1));
        %plot(lowerpoints(xendpoint2:end,2), -lowerpoints(xendpoint2:end,1));
    %end
    %axis equal;
    %hold off;
    
    %hold on; 
    %plot(x1,y1);
    %plot(x2,y2);
    %axis equal;
    %hold off;
    
    upperpoints2(:,1,im) = x1./editratio;
    lowerpoints2(:,1,im) = x2./editratio;
    upperpoints2(:,2,im) = y1./editratio;
    lowerpoints2(:,2,im) = y2./editratio;
    
    for vidcounter1 = 1 :midcounter4(1)
        for vidcounter2 = 1 :2
            upperpoints3(vidcounter1,vidcounter2,im) = upperpoints2(vidcounter1,vidcounter2,im);
        end
    end
    
    for vidcounter1 = 1 :midcounter4(1)-1
        decreup3(vidcounter1) = upperpoints3(midcounter4(1)-vidcounter1+1,1,im) - upperpoints3(midcounter4(1)-vidcounter1,1,im);
    end
    
    for vidcounter1 = midcounter4(1) : size(x1)
        for vidcounter2 = 1 :2
            upperpoints4(vidcounter1-midcounter4(1)+1,vidcounter2,im) = upperpoints2(vidcounter1,vidcounter2,im);
        end
    end
    
    upperpoints4sp(1,1,im) = upperpoints4(1,1,im);
    
    for vidcounter1 = 1 :midcounter4(1)-1
        upperpoints4sp(vidcounter1+1,1,im) = upperpoints4(vidcounter1)+decreup3(vidcounter1);
    end
    
    upperpoints4sp(:,2,im) = spline(upperpoints4(:,1,im),upperpoints4(:,2,im), upperpoints4sp(:,1,im));
    
    upperpoints4 = upperpoints4sp;
    
    for vidcounter3 = 1 : size(upperpoints4,1)
        for vidcounter4 = 1 :2
            upperpoints5(vidcounter3,vidcounter4,im) = upperpoints4(size(upperpoints4,1)-vidcounter3+1,vidcounter4,im);
        end
    end
    
    for vidcounter5 = 1 : size(upperpoints3,1)
        for vidcounter6 = 1 : 2
            if vidcounter6 ==1
                upperpoints6(vidcounter5,vidcounter6,im) = upperpoints3(vidcounter5,vidcounter6,im);
            else
                upperpoints6(vidcounter5,vidcounter6,im) = (upperpoints3(vidcounter5,vidcounter6,im)+upperpoints5(vidcounter5,vidcounter6,im))/2;
            end
        end
    end
    
    for vidcounter1 = 1 :midcounter4(1)
        for vidcounter2 = 1 :2
            lowerpoints3(vidcounter1,vidcounter2,im) = lowerpoints2(vidcounter1,vidcounter2,im);
        end
    end
    
    for vidcounter1 = 1 :midcounter4(1)-1
        decrelow3(vidcounter1) = lowerpoints3(midcounter4(1)-vidcounter1+1,1,im) - lowerpoints3(midcounter4(1)-vidcounter1,1,im);
    end
    
    for vidcounter1 = midcounter4(1) : size(x1)
        for vidcounter2 = 1 :2
            lowerpoints4(vidcounter1-midcounter4(1)+1,vidcounter2,im) = lowerpoints2(vidcounter1,vidcounter2,im);
        end
    end
    
    lowerpoints4sp(1,1,im) = lowerpoints4(1,1,im);
    
    for vidcounter1 = 1 :midcounter4(1)-1
        lowerpoints4sp(vidcounter1+1,1,im) = lowerpoints4(vidcounter1)+decrelow3(vidcounter1);
    end
    
    lowerpoints4sp(:,2,im) = spline(lowerpoints4(:,1,im),lowerpoints4(:,2,im), lowerpoints4sp(:,1,im));
    
    lowerpoints4 = lowerpoints4sp;
    
    for vidcounter3 = 1 : size(lowerpoints4,1)
        for vidcounter4 = 1 :2
            lowerpoints5(vidcounter3,vidcounter4,im) = lowerpoints4(size(lowerpoints4,1)-vidcounter3+1,vidcounter4,im);
        end
    end
    
    for vidcounter5 = 1 : size(lowerpoints3,1)
        for vidcounter6 = 1 : 2
            if vidcounter6 ==1
                lowerpoints6(vidcounter5,vidcounter6,im) = lowerpoints3(vidcounter5,vidcounter6,im);
            else
                lowerpoints6(vidcounter5,vidcounter6,im) = (lowerpoints3(vidcounter5,vidcounter6,im)+lowerpoints5(vidcounter5,vidcounter6,im))/2;
            end
        end
    end
    
    for im9 = 1 :midcounter4(1)-1
        x3(im9) = x1(im9);
    end
    
    for im10 = midcounter4(1)+1 : size(x1,1)
        x5(im10-midcounter4) = x1(im10);
    end
    
    x4 = x1(midcounter4);

    x3 = x3./editratio;
    x4 = x4./editratio;
    x5 = x5./editratio;

    x4a = x4.*editratio;
    
    minsplinealpha = x3(1);
    maxsplinealpha = x4(1);
    sizesplinealpha = floor(numx/2);
    totallength = maxsplinealpha - minsplinealpha; 
    
%     alphalength(1) = totallength/(1-2*(1-(1/((1-sin(anglez))^(sizesplinealpha-1)))));
%     
%     alphalength(2) = 2*alphalength(1)*sin(anglez)/(1-sin(anglez));
%       
%     for alphacount = 3:sizesplinealpha
%         alphalength(alphacount) = alphalength(alphacount-1) + alphalength(alphacount-1)*sin(anglez)/(1-sin(anglez));
%     end

    for alphacount = 1:sizesplinealpha
        alphalength(alphacount) = totallength/sizesplinealpha;
    end
    
    for alphacount2 = 1 :size(alphalength,2)
        alphalengthsum(alphacount2) = 0;
    end
    
    for alphacount3 = 1 : size(alphalength,2)
        for alphacount4 = 1:alphacount3
            alphalengthsum(alphacount3) = alphalengthsum(alphacount3)+alphalength(alphacount4);
        end
        alphalength2(sizesplinealpha-alphacount3+1) = maxsplinealpha - alphalengthsum(alphacount3);
    end
    
    x3a = alphalength2.*editratio;
    
    minsplinealpha = x4(1);
    maxsplinealpha = x5(size(x5,2));
    sizesplinealpha = floor(numx/2);
    totallength = maxsplinealpha - minsplinealpha; 
    
%     alphalength(1) = totallength/(1-2*(1-(1/((1-sin(anglez))^(sizesplinealpha-1)))));
%     
%     alphalength(2) = 2*alphalength(1)*sin(anglez)/(1-sin(anglez));
%       
%     for alphacount = 3:sizesplinealpha
%         alphalength(alphacount) = alphalength(alphacount-1) + alphalength(alphacount-1)*sin(anglez)/(1-sin(anglez));
%     end
    
    for alphacount = 1:sizesplinealpha
        alphalength(alphacount) = totallength/sizesplinealpha;
    end
    
    for alphacount2 = 1 :size(alphalength,2)
        alphalengthsum(alphacount2) = 0;
    end
    
    for alphacount3 = 1 : size(alphalength,2)
        for alphacount4 = 1:alphacount3
            alphalengthsum(alphacount3) = alphalengthsum(alphacount3)+alphalength(alphacount4);
        end
        alphalength2(alphacount3) = minsplinealpha + alphalengthsum(alphacount3);
    end
    
    x5a = alphalength2.*editratio;
    
    for lengthcountfinal1 = 1 : size(x3a,2)
        xalpha(lengthcountfinal1) = x3a(lengthcountfinal1);
    end
    
    xalpha(size(x3a,2)+1) = x4a;
    
    for lengthcountfinal2 = 1 : size(x5a,2)
        xalpha(lengthcountfinal2+size(x3a,2)+1) = x5a(lengthcountfinal2);
    end
    
    x1 = xalpha;
    x2 = xalpha;
    
    y1 = -spline(upperpoints(:,2),upperpoints(:,1), x1);
    y2 = -spline(lowerpoints(:,2),lowerpoints(:,1), x2);
    
    %if im ==1
    %hold on;
    %scatter(x1, y1);
    %scatter(x2, y2);
    %plot(upperpoints(:,2), -upperpoints(:,1));
    %plot(lowerpoints(:,2), -lowerpoints(:,1));
    %grid on;
    %axis equal;
    %hold off;
    %end
    
    x1a = xalpha./editratio;
    x2a = xalpha./editratio;
    
    x3 = xalpha(1:floor(numx/2))./editratio;
    x4 = xalpha(floor(numx/2)+1)./editratio;
    if mod(numx,2)==0
        x5 = xalpha(floor(numx/2)+2:numx+1)./editratio;
    else
        x5 = xalpha(floor(numx/2)+2:numx)./editratio;
    end
    
    for im6 = 1 : floor(numx/2)
        y3(im6) = y1(im6)/editratio;
        y4(im6) = y2(im6)/editratio;
    end
        
    ylength = (y4 - y3)/numy2;
    yspecial = zeros(numy, floor(numx/2));
    yspecial(1,:) = y3;
    yspecial(numy,:) = y4;
    for yspecialcounter = 2 : numy2  
        yspecial(yspecialcounter,:) = (yspecialcounter-1) * ylength+y3;
    end

    ymid1 = y1(floor(numx/2)+1)/editratio;
    ymid2 = y2(floor(numx/2)+1)/editratio;
        
    ymidlength = (ymid2-ymid1)/numy2;
    ymidspecial(1,1) = ymid1;
    ymidspecial(numy,1) = ymid2;
    for yspecialmidcounter = 2:numy2
        ymidspecial(yspecialmidcounter,1) = (yspecialmidcounter-1)* ymidlength+ymid1;
    end
    
    if mod(numx,2)==0
        for im7 = floor(numx/2)+2 : numx+1
            y5(im7-(floor(numx/2)+1)) = y1(im7)/editratio;
            y6(im7-(floor(numx/2)+1)) = y2(im7)/editratio;
        end
    else
        for im7 = floor(numx/2)+2 : numx
            y5(im7-(floor(numx/2)+1)) = y1(im7)/editratio;
            y6(im7-(floor(numx/2)+1)) = y2(im7)/editratio;
        end
    end    
    
    ylength2 = (y6 - y5)/numy2;
    yspecial2 = zeros(numy, floor(numx/2));
    yspecial2(1,:) = y5;
    yspecial2(numy,:) = y6;

    for yspecialcounter2 = 2 : numy2  
        yspecial2(yspecialcounter2,:) = (yspecialcounter2-1) * ylength2+y5;
    end
    
%     hold on; 
%     axis equal;
%     for yspecialplotcounter = 1 : numy
%         plot(x3-x4, yspecial(yspecialplotcounter, :), 'r*','MarkerSize',4);
%         plot(x5-x4, yspecial2(yspecialplotcounter, :), 'b*','MarkerSize',4);
%         %plot(-x5+x4, yspecial2(yspecialplotcounter, :), 'b*','MarkerSize',4);
%     end
%     hold off;

    mean = 0;
    sigma = rairjet;
    
    airjetpshape1x = x3-x4;
    airjetpshape1x(size(x3,2)+1) = 0;

    airjetpshape1y = normpdf(airjetpshape1x,mean, sigma)./normpdf(mean,mean, sigma).*8000;

    for airjetpmatcount1 = 1 : (size(airjetpshape1y,2)-1)
        airjetpmat1x(airjetpmatcount1) = (airjetpshape1x(airjetpmatcount1) + airjetpshape1x(airjetpmatcount1+1))/2;
        airjetpmat1y(airjetpmatcount1) = (airjetpshape1y(airjetpmatcount1) + airjetpshape1y(airjetpmatcount1+1))/2;
    end
    
    %scatter(airjetpshape1x,airjetpshape1y);
    %scatter(airjetpmat1x,airjetpmat1y);
    
    for yspecialcounter3 = 1 : numy
        for yspecialcounter4 = 1: floor(numx/2)
            yspecial3(yspecialcounter3, yspecialcounter4) = yspecial2(yspecialcounter3, floor(numx/2)+1-yspecialcounter4);
        end
    end
    
    yspecial4 = (yspecial3+yspecial)/2;
    
    yspecialmin = yspecial4(1,1);
    
    for yspecialrow = 1 : size(yspecial4,1)
        for yspecialcol = 1 : size(yspecial4,2)
            if yspecialmin >= yspecial4(yspecialrow, yspecialcol)
                yspecialmin = yspecial4(yspecialrow, yspecialcol);
            end
        end
    end
    
    yspecial4 = yspecial4 - yspecialmin;
    
    alteredy = yspecial4;
    
%     for frame0 = 1:5
%         frame = frame0*20;
%         if frame == im 
%             figure(im)
%             hold on; 
%             axis equal;
%             for yspecialplotcounter = 1 : numy
%                 plot(x3-x4, yspecial4(yspecialplotcounter, :), 'r*','MarkerSize',4);    
%             end
%             hold off;
%         end
%     end
    
    ymidspecial = ymidspecial - yspecialmin;
    
    alteredymid = ymidspecial;
    
    datafinal = zeros(numy+1, floor(numx/2));
    datafinal(1,:) = x3-x4;

    for yspecialcounter4 = 2 : numy+1
        datafinal(yspecialcounter4,:) = yspecial4(yspecialcounter4-1,:);
    end

    for zspecialcounter1 = 1 : floor(numx/2)
        for zspecialcounter2 = 1 : numz3
            alteredz(zspecialcounter1, zspecialcounter2) = -datafinal(1, zspecialcounter1)*sin(anglez*(zspecialcounter2-1));
            alteredx(zspecialcounter2, zspecialcounter1) = datafinal(1, zspecialcounter1)*cos(anglez*(zspecialcounter2-1));
        end
    end
    
    enlargecount = 0.005;
    
    palteredx(1) = datafinal(1, 1)-enlargecount;
    palteredx(2) = (datafinal(1, 1)-enlargecount)*cos(anglez*(numz3-1));

    palteredx(3) = 0+enlargecount;
    palteredx(4) = (0+enlargecount)*cos(anglez*(numz3-1));

    palteredz(1) = -(datafinal(1, 1)-enlargecount)*sin(anglez*(1-1));
    palteredz(2) = -(datafinal(1, 1)-enlargecount)*sin(anglez*(numz3-1));

    palteredz(3) = -(0+enlargecount)*sin(anglez*(1-1));
    palteredz(4) = -(0+enlargecount)*sin(anglez*(numz3-1));

    palteredy(1) = 0 - enlargecount;
    palteredy(2) = alteredymid(1) + enlargecount; 
    
    for im2 = 1 : size(alteredx,1)
        for im3 = 1 : size(alteredx,2)
            alteredx2(im2,im3,im) = alteredx(im2,im3);
        end
    end
    
    for im2 = 1 : size(alteredy,1)
        for im3 = 1 : size(alteredy,2)
            alteredy2(im2,im3,im) = alteredy(im2,im3);
        end
    end
    
    for im2 = 1 : size(alteredz,1)
        for im3 = 1 : size(alteredz,2)
            alteredz2(im2,im3,im) = alteredz(im2,im3);
        end
    end
    
    for im2 = 1 : size(alteredymid,1)
        for im3 = 1 : size(alteredymid,2)
            alteredymid2(im2,im3,im) = alteredymid(im2,im3);
        end
    end
    
    nodecount = 0;

    for printcount3 = 1 : numz3
        for printcount1 = 1 : floor(numx/2)
            for printcount2 = 1 : numy
                nodecount = nodecount + 1;
                origmatrix(nodecount, 1, im) = alteredx2(printcount3,printcount1,im);
                origmatrix(nodecount, 2, im) = alteredy2(printcount2,printcount1,im);
                origmatrix(nodecount, 3, im) = alteredz2(printcount1, printcount3,im);
                fprintf(file1, '%d,%d,%d\n',origmatrix(nodecount, 1, im), origmatrix(nodecount, 2, im), origmatrix(nodecount, 3, im));
            end
        end
    end

    for printcountmid = 1 : numy
        nodecount = nodecount + 1;
        origmatrix(nodecount, 1, im) = 0;
        origmatrix(nodecount, 2, im) = alteredymid2(printcountmid,1,im);
        origmatrix(nodecount, 3, im) = 0;
        fprintf(file1, '%d,%d,%d\n',origmatrix(nodecount, 1, im), origmatrix(nodecount, 2, im), origmatrix(nodecount, 3, im));
    end

    shiftx = x4;
    shifty = yspecialmin;
    
    upperpoints3(:,1,im)= upperpoints3(:,1,im)-shiftx;
    upperpoints5(:,1,im)= upperpoints5(:,1,im)-shiftx;
    upperpoints6(:,1,im)= upperpoints6(:,1,im)-shiftx;

    upperpoints3(:,2,im)= upperpoints3(:,2,im)-shifty;
    upperpoints5(:,2,im)= upperpoints5(:,2,im)-shifty;
    upperpoints6(:,2,im)= upperpoints6(:,2,im)-shifty;
        
    lowerpoints3(:,1,im)= lowerpoints3(:,1,im)-shiftx;
    lowerpoints5(:,1,im)= lowerpoints5(:,1,im)-shiftx;
    lowerpoints6(:,1,im)= lowerpoints6(:,1,im)-shiftx;

    lowerpoints3(:,2,im)= lowerpoints3(:,2,im)-shifty;
    lowerpoints5(:,2,im)= lowerpoints5(:,2,im)-shifty;
    lowerpoints6(:,2,im)= lowerpoints6(:,2,im)-shifty;
    
%     for frame0 = 1:5
%         frame = frame0*25;
%         if im == frame
%          
%             figure(frame)
%             hold on; 
%             plot(lowerpoints3(:,1,frame),lowerpoints3(:,2,frame), 'r');
%             plot(lowerpoints3(:,1,frame),lowerpoints5(:,2,frame), 'b');
%             plot(lowerpoints6(:,1,frame),lowerpoints6(:,2,frame), 'g');
%             plot(upperpoints3(:,1,frame),upperpoints3(:,2,frame), 'r');
%             plot(upperpoints3(:,1,frame),upperpoints5(:,2,frame), 'b');
%             plot(upperpoints6(:,1,frame),upperpoints6(:,2,frame), 'g');
%             axis equal;
%             hold off;
%         end
%     end
    
end

%for frame = 73:73;
    %figure(1)
    %axis equal;        
    %scatter3(origmatrix(:,1,frame), origmatrix(:,3,frame), origmatrix(:,2,frame), 4);
    %pause(1)
%end
  
filevid = fopen(fnamevid, 'wt');

for im = 1:numIm
    for printcountvid1 = 1 : size(upperpoints6,1)
        fprintf(filevid, '%f,%f,%f,%f\n',upperpoints6(printcountvid1,1,im), upperpoints6(printcountvid1,2,im), lowerpoints6(printcountvid1,1,im), lowerpoints6(printcountvid1,2,im) );
    end
end

fname3 = char(strcat(fname1, '.feb'));
file2 = fopen(fname3, 'wt');

if currentos == 1
    fname4=char(strcat('febio',{' '},fname3));
else
    if currentos == 0
        fname4=char(strcat('febio.lnx64',{' '},fname3));
    end
end

fprintf(file2, '<?xml version="1.0" encoding="ISO-8859-1"?>\n');
fprintf(file2, '<febio_spec version="1.2">\n');
fprintf(file2, '\t<Globals>\n');
fprintf(file2, '\t\t<Constants>\n');
fprintf(file2, '\t\t\t<T>0</T>\n');
fprintf(file2, '\t\t\t<R>0</R>\n');
fprintf(file2, '\t\t\t<Fc>0</Fc>\n');
fprintf(file2, '\t\t</Constants>\n');
fprintf(file2, '\t</Globals>\n');
fprintf(file2, '\t<Material>\n');

if visco == 1
    
    fprintf(file2, '\t\t<material id="1" name="Material1" type="%s">\n', viscotype);    
    for viscounter = 1:numberofvisco
        fprintf(file2, '\t\t\t<g%d>%f</g%d>\n',viscounter, gvar(viscounter), viscounter);
        fprintf(file2, '\t\t\t<t%d>%f</t%d>\n',viscounter, tvar(viscounter), viscounter);
    end
    fprintf(file2, '\t\t\t<elastic type="%s">\n', materialtype);
    fprintf(file2, '\t\t\t\t<density>%e</density>\n', density);
    if materialvar == 1
        fprintf(file2, '\t\t\t\t<E>%e</E>\n', var1);
        fprintf(file2, '\t\t\t\t<v>%e</v>\n', var2);
    end
    if materialvar == 2
        fprintf(file2, '\t\t\t\t<c1>%e</c1>\n', var1);
        fprintf(file2, '\t\t\t\t<c2>%e</c2>\n', var2);
        fprintf(file2, '\t\t\t\t<k>%e</k>\n', var3);
    end
    fprintf(file2, '\t\t\t</elastic>\n');
    
else

    fprintf(file2, '\t\t<material id="1" name="Material1" type="%s">\n', materialtype);
    fprintf(file2, '\t\t\t<density>%e</density>\n', density);
    if materialvar == 1
        fprintf(file2, '\t\t\t<E>%e</E>\n', var1);
        fprintf(file2, '\t\t\t<v>%e</v>\n', var2);
    end
    if materialvar == 2
        fprintf(file2, '\t\t\t<c1>%e</c1>\n', var1);
        fprintf(file2, '\t\t\t<c2>%e</c2>\n', var2);
        fprintf(file2, '\t\t\t<k>%e</k>\n', var3);
    end
    
end

fprintf(file2, '\t\t</material>\n');

fprintf(file2, '\t\t<material id="2" name="Symmetry" type="rigid body">\n');
fprintf(file2, '\t\t\t<density>%e</density>\n', density);
fprintf(file2, '\t\t\t<center_of_mass>0,0,0</center_of_mass>\n');
fprintf(file2, '\t\t</material>\n');

fprintf(file2, '\t</Material>\n');
fprintf(file2, '\t<Geometry>\n');
fprintf(file2, '\t\t<Nodes>\n');
            
printcountspecial = 0;

for printcount3 = 1 : numz3
    for printcount1 = 1 : floor(numx/2)
        for printcount2 = 1 : numy
                printcountspecial = printcountspecial + 1;
                fprintf(file2, '\t\t\t<node id="%d"> %e,%e, %e</node>\n', printcountspecial, alteredx2(printcount3,printcount1,1), alteredy2(printcount2,printcount1,1), alteredz2(printcount1, printcount3,1));
                %testmatrix(printcountspecial, 1) = alteredx2(printcount3,printcount1,1);
                %testmatrix(printcountspecial, 2) = alteredy2(printcount2,printcount1,1);
                %testmatrix(printcountspecial, 3) = alteredz2(printcount1, printcount3,1);
        end
    end
end

midcountindic = printcountspecial+1;

for printcountmid = 1 : numy
    printcountspecial = printcountspecial + 1;
    fprintf(file2, '\t\t\t<node id="%d"> %e,%e, %e</node>\n', printcountspecial, 0, alteredymid2(printcountmid,1,1), 0);
    %testmatrix(printcountspecial, 1) = 0;
    %testmatrix(printcountspecial, 2) = alteredymid2(printcountmid,1,1);
    %testmatrix(printcountspecial, 3) = 0;
end

%hold on;
%scatter3(testmatrix(:,1), testmatrix(:,3), testmatrix(:,2), 4);
%axis equal;
%hold off;

fprintf(file2, '\t\t\t<node id="%d"> %e,%e, %e</node>\n', totalnumnodes+1, palteredx(2), palteredy(2), palteredz(2));
fprintf(file2, '\t\t\t<node id="%d"> %e,%e, %e</node>\n', totalnumnodes+2, palteredx(2), palteredy(1), palteredz(2));
fprintf(file2, '\t\t\t<node id="%d"> %e,%e, %e</node>\n', totalnumnodes+3, palteredx(4), palteredy(2), palteredz(4));
fprintf(file2, '\t\t\t<node id="%d"> %e,%e, %e</node>\n', totalnumnodes+4, palteredx(4), palteredy(1), palteredz(4));


fprintf(file2, '\t\t</Nodes>\n');
fprintf(file2, '\t\t<Elements>\n');

elecounter1 = 0;

for elecountz = 1:(numz3-1)
    for elecountx = 1 : (floor(numx/2)-1)
        for elecounty = 1 : (numy-1)
            
            eleprinter1 = numy*(elecountx)+elecounty+numy*floor(numx/2)*elecountz;
            eleprinter2 = numy*(elecountx)+elecounty+numy*floor(numx/2)*(elecountz-1);
            eleprinter3 = eleprinter2 +1;
            eleprinter4 = eleprinter1 +1;
            eleprinter5 = numy*(elecountx-1)+elecounty+numy*floor(numx/2)*elecountz;
            eleprinter6 = numy*(elecountx-1)+elecounty+numy*floor(numx/2)*(elecountz-1);
            eleprinter7 = eleprinter6 +1;
            eleprinter8 = eleprinter5 +1;
            
            elecounter1 = elecounter1+1;
            
            fprintf(file2, '\t\t\t<hex8 id="%d" mat="1">\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d</hex8>\n', elecounter1, eleprinter1, eleprinter2, eleprinter3, eleprinter4, eleprinter5, eleprinter6, eleprinter7, eleprinter8);
            
        end
    end
end

for elecountz = 1:numz3-1
    for elecountx = floor(numx/2) : floor(numx/2)
        for elecounty = 1 : (numy-1)
            
            eleprinter1 = numy*(elecountx-1)+elecounty+1+ numy*floor(numx/2)*(elecountz);
            eleprinter2 = numy*floor(numx/2)*numz3+elecounty+1;
            eleprinter3 = numy*(elecountx-1)+elecounty+1+ numy*floor(numx/2)*(elecountz-1);
            eleprinter4 = numy*(elecountx-1)+elecounty+ numy*floor(numx/2)*(elecountz);
            eleprinter5 = numy*floor(numx/2)*numz3+elecounty;
            eleprinter6 = numy*(elecountx-1)+elecounty+ numy*floor(numx/2)*(elecountz-1);

            elecounter1 = elecounter1+1;
            
            fprintf(file2, '\t\t\t<penta6 id="%d" mat="1">\t%d,\t%d,\t%d,\t%d,\t%d,\t%d</penta6>\n', elecounter1, eleprinter1, eleprinter2, eleprinter3, eleprinter4, eleprinter5, eleprinter6);
            
        end
    end
end

elecounter1 = elecounter1+1;
fprintf(file2, '\t\t\t<quad4 id="%d" mat="2">\t%d,\t%d,\t%d,\t%d,</quad4>\n', elecounter1, totalnumnodes+2, totalnumnodes+1, totalnumnodes+3, totalnumnodes+4);

fprintf(file2, '\t\t</Elements>\n');

fprintf(file2, '\t\t<ElementData>\n');
fprintf(file2, '\t\t\t<element id = "%d">\n', elecounter1);
fprintf(file2, '\t\t\t\t<thickness>0.01,0.01,0.01,0.01</thickness>\n');
fprintf(file2, '\t\t\t</element>\n');
fprintf(file2, '\t\t</ElementData>\n');


fprintf(file2, '\t</Geometry>\n');
fprintf(file2, '\t<Boundary>\n');
fprintf(file2, '\t\t<fix>\n');

for elecountz = 1:numz3
    for elecountx = 1 : 1
        for elecounty = numy-1 : numy-1
            eleprinterboundary = numy*(elecountx-1)+elecounty+1+numy*floor(numx/2)*(elecountz-1);
            fprintf(file2, '\t\t\t<node id="%d" bc="xyz"/>\n', eleprinterboundary);
        end
    end
end

for elecountz = 1:1
    for elecountx = 1 : floor(numx/2)
        for elecounty = 1 : numy
            eleprinterboundary = numy*(elecountx-1)+elecounty+numy*floor(numx/2)*(elecountz-1);
            fprintf(file2, '\t\t\t<node id="%d" bc="z"/>\n', eleprinterboundary);
        end
    end
end

for eleprinterboundary = (floor(numx/2)*numy*numz3+1) : (floor(numx/2)*numy*numz3+numy)
    fprintf(file2, '\t\t\t<node id="%d" bc="xz"/>\n', eleprinterboundary);
end

fprintf(file2, '\t\t</fix>\n');

fprintf(file2, '\t<contact type="sliding-tension-compression">\n');
fprintf(file2, '\t\t<laugon>0</laugon>\n');
fprintf(file2, '\t\t<tolerance>0.2</tolerance>\n');
fprintf(file2, '\t\t<gaptol>0</gaptol>\n');
fprintf(file2, '\t\t<penalty>%e</penalty>\n', penalty);
fprintf(file2, '\t\t<auto_penalty>1</auto_penalty>\n');
fprintf(file2, '\t\t<two_pass>0</two_pass>\n');
fprintf(file2, '\t\t<search_tol>0.01</search_tol>\n');
fprintf(file2, '\t\t<symmetric_stiffness>0</symmetric_stiffness>\n');
fprintf(file2, '\t\t<search_radius>1</search_radius>\n');
fprintf(file2, '\t\t<seg_up>0</seg_up>\n');
fprintf(file2, '\t\t<tension>1</tension>\n');
fprintf(file2, '\t\t<minaug>0</minaug>\n');
fprintf(file2, '\t\t<maxaug>10</maxaug>\n');

fprintf(file2, '\t\t<surface type="master">\n');
fprintf(file2, '\t\t\t<quad4 id="1" mat="2">\t%d,\t%d,\t%d,\t%d,</quad4>\n', totalnumnodes+2, totalnumnodes+1, totalnumnodes+3, totalnumnodes+4);
fprintf(file2, '\t\t</surface>\n');

fprintf(file2, '\t\t<surface type="slave">\n');

elecounter1 = 0;
for elecountx = 1 : floor(numx/2)-1
    for elecounty = 1 : numy-1
        elecounter1 = elecounter1 +1;
        eleprinter1 = numy*(elecountx)+elecounty+numy*floor(numx/2)*(numz3-1)+1;
        eleprinter2 = eleprinter1-1; 
        eleprinter3 = numy*(elecountx-1)+elecounty+numy*floor(numx/2)*(numz3-1);
        eleprinter4 = eleprinter3+1;
        fprintf(file2, '\t\t\t<quad4 id="%d" mat="2">\t%d,\t%d,\t%d,\t%d,</quad4>\n', elecounter1, eleprinter1, eleprinter2, eleprinter3, eleprinter4);
    end
end

for elecountx = floor(numx/2) : floor(numx/2)
    for elecounty = 1 : numy-1
        elecounter1 = elecounter1 +1;
        eleprinter1 = numy*(elecountx)+elecounty+numy*floor(numx/2)*(numz3-1)+1;
        eleprinter2 = eleprinter1-1; 
        eleprinter3 = numy*(elecountx-1)+elecounty+numy*floor(numx/2)*(numz3-1);
        eleprinter4 = eleprinter3+1;
        fprintf(file2, '\t\t\t<quad4 id="%d" mat="2">\t%d,\t%d,\t%d,\t%d,</quad4>\n', elecounter1, eleprinter1, eleprinter2, eleprinter3, eleprinter4);
    end
end

fprintf(file2, '\t\t</surface>\n');
fprintf(file2, '\t</contact>\n');

fprintf(file2, '\t</Boundary>\n');
fprintf(file2, '\t<Loads>\n');
fprintf(file2, '\t\t<pressure>\n');

elecounter1 =0;

for elecountz = 1:(numz3-1)
    for elecountx = 1 : (floor(numx/2)-1)
        for elecounty = numy-1 : numy-1
            
            elecounter1 = elecounter1+1;
            
            eleprinter1 = numy*(elecountx)+elecounty+1+numy*floor(numx/2)*(elecountz-1);
            eleprinter2 = numy*(elecountx)+elecounty+1+numy*floor(numx/2)*(elecountz);
            eleprinter3 = numy*(elecountx-1)+elecounty+1+numy*floor(numx/2)*(elecountz);
            eleprinter4 = numy*(elecountx-1)+elecounty+1+numy*floor(numx/2)*(elecountz-1);
            
            fprintf(file2, '\t\t\t<quad4 id="%d" lc="1" scale="%d">\t%d,\t%d,\t%d,\t%d</quad4>\n', elecounter1, IOP, eleprinter1, eleprinter2, eleprinter3, eleprinter4);
            
        end
    end
end

for elecountz = 1:numz3-1
    for elecountx = floor(numx/2) : floor(numx/2)
        for elecounty = numy-1 : (numy-1)
            
            elecounter1 = elecounter1+1;
            
            eleprinter1 = numy*(elecountx-1)+elecounty+1+ numy*floor(numx/2)*(elecountz);
            eleprinter2 = numy*(elecountx-1)+elecounty+1+ numy*floor(numx/2)*(elecountz-1);
            eleprinter3 = numy*floor(numx/2)*numz3+elecounty+1;
            
            fprintf(file2, '\t\t\t<tri3 id="%d" lc="1" scale="%d">\t%d,\t%d,\t%d</tri3>\n', elecounter1, IOP, eleprinter1, eleprinter2, eleprinter3);
            
        end
    end
end

fprintf(file2, '\t\t</pressure>\n');
fprintf(file2, '\t</Loads>\n');

fprintf(file2, '\t<Constraints>\n');
fprintf(file2, '\t\t<rigid_body mat="2">\n');
fprintf(file2, '\t\t\t<trans_x type="fixed"/>\n');
fprintf(file2, '\t\t\t<trans_y type="fixed"/>\n');
fprintf(file2, '\t\t\t<trans_z type="fixed"/>\n');
fprintf(file2, '\t\t\t<rot_x type="fixed"/>\n');
fprintf(file2, '\t\t\t<rot_y type="fixed"/>\n');
fprintf(file2, '\t\t\t<rot_z type="fixed"/>\n');
fprintf(file2, '\t\t</rigid_body>\n');
fprintf(file2, '\t</Constraints>\n');

fprintf(file2, '\t<LoadData>\n');

fprintf(file2, '\t\t<loadcurve id="1" type="smooth">\n');

fprintf(file2, '\t\t\t<loadpoint>0,1</loadpoint>\n');
fprintf(file2, '\t\t\t<loadpoint>1,1</loadpoint>\n');

fprintf(file2, '\t\t</loadcurve>\n');
fprintf(file2, '\t\t<loadcurve id="2" type="smooth">\n');

for loadcount1 = 1 :size(loadfunction,1)
    fprintf(file2, '\t\t\t<loadpoint>%f,%f</loadpoint>\n', loadfunction(loadcount1,1),loadfunction(loadcount1,2));
end

fprintf(file2, '\t\t</loadcurve>\n');

fprintf(file2, '\t\t<loadcurve id="3" type="step">\n');

for loadcount2 = 1 :timestepmodel
    fprintf(file2, '\t\t\t<loadpoint>%f,1</loadpoint>\n', loadcount2*timestepsize);
end

fprintf(file2, '\t\t</loadcurve>\n');

fprintf(file2, '\t</LoadData>\n');
fprintf(file2, '\t<Output>\n');
fprintf(file2, '\t\t<plotfile type="febio">\n');
fprintf(file2, '\t\t\t<var type="acceleration"/>\n');
fprintf(file2, '\t\t\t<var type="displacement"/>\n');
fprintf(file2, '\t\t\t<var type="stress"/>\n');
fprintf(file2, '\t\t\t<var type="velocity"/>\n');
fprintf(file2, '\t\t</plotfile>\n');
fprintf(file2, '\t\t<logfile>\n');
fprintf(file2, '\t\t\t<node_data data="x;y;z" delim=", " file="%sPosition_%s.txt"/>\n', filelocation1, fname1);
fprintf(file2, '\t\t\t<node_data data="ux;uy;uz" delim=", " file="%sDisplacement_%s.txt"/>\n', filelocation1, fname1);
fprintf(file2, '\t\t</logfile>\n');
fprintf(file2, '\t</Output>\n');
fprintf(file2, '\t<Step name="Step01">\n');
fprintf(file2, '\t\t<Module type="solid"/>\n');
fprintf(file2, '\t\t<Control>\n');
fprintf(file2, '\t\t\t<time_steps>%d</time_steps>\n', timestepmodel);
fprintf(file2, '\t\t\t<step_size>%f</step_size>\n', timestepsize);
fprintf(file2, '\t\t\t<max_refs>15</max_refs>\n');
fprintf(file2, '\t\t\t<max_ups>10</max_ups>\n');
fprintf(file2, '\t\t\t<dtol>0.001</dtol>\n');
fprintf(file2, '\t\t\t<etol>0.01</etol>\n');
fprintf(file2, '\t\t\t<rtol>0</rtol>\n');
fprintf(file2, '\t\t\t<lstol>0.9</lstol>\n');

fprintf(file2, '\t\t\t<time_stepper>\n');
fprintf(file2, '\t\t\t\t<dtmin>0.01</dtmin>\n');
fprintf(file2, '\t\t\t\t<dtmax lc="3"></dtmax>\n');
fprintf(file2, '\t\t\t\t<max_retries>5</max_retries>\n');
fprintf(file2, '\t\t\t\t<opt_iter>10</opt_iter>\n');
fprintf(file2, '\t\t\t</time_stepper>\n');

fprintf(file2, '\t\t\t<analysis type="%s"/>\n', analysismethod);
fprintf(file2, '\t\t</Control>\n');
fprintf(file2, '\t\t<Loads>\n');
fprintf(file2, '\t\t\t<pressure>\n');

elecounter1 = 0;

for elecountz = 1:(numz3-1)
    for elecountx = 1 : (floor(numx/2)-1) %1 : (floor(numx/2)-1)
        for elecounty = 1 : 1
            
            elecounter1 = elecounter1+1;
            
            eleprinter1 = numy*(elecountx)+elecounty+numy*floor(numx/2)*(elecountz);
            eleprinter2 = numy*(elecountx)+elecounty+numy*floor(numx/2)*(elecountz-1);
            eleprinter3 = numy*(elecountx-1)+elecounty+numy*floor(numx/2)*(elecountz-1);
            eleprinter4 = numy*(elecountx-1)+elecounty+numy*floor(numx/2)*(elecountz);
            
            fprintf(file2, '\t\t\t\t<quad4 id="%d" lc="2" scale="%d">\t%d,\t%d,\t%d,\t%d</quad4>\n', elecounter1, airjetpmat1y(elecountx), eleprinter1, eleprinter2, eleprinter3, eleprinter4);
            
        end
    end
end

for elecountz = 1:numz3-1
    for elecountx = floor(numx/2) : floor(numx/2)
        for elecounty = 1 : 1
            
            elecounter1 = elecounter1+1;
            
            eleprinter1 = numy*(elecountx-1)+elecounty+ numy*floor(numx/2)*(elecountz);
            eleprinter2 = numy*floor(numx/2)*numz3+elecounty;
            eleprinter3 = numy*(elecountx-1)+elecounty+ numy*floor(numx/2)*(elecountz-1);
            
            fprintf(file2, '\t\t\t\t<tri3 id="%d" lc="2" scale="%d">\t%d,\t%d,\t%d</tri3>\n', elecounter1, airjetpmat1y(floor(numx/2)), eleprinter1, eleprinter2, eleprinter3);
            
        end
    end
end

fprintf(file2, '\t\t\t</pressure>\n');
fprintf(file2, '\t\t</Loads>\n');
fprintf(file2, '\t</Step>\n');
fprintf(file2, '</febio_spec>\n');

[status1 cmdout1]=system(fname4);

fclose(file1);
fclose(fileinitial);
fclose(filetimestepchoice);
fclose(fileload);
fclose(filevid);
fclose(file2);

end

