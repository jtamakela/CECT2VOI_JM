function [RESULT_PROFILES50, RESULT_PROFILES90] = CECTVOI_Analysis
%% m-file for analysing CECT image VOIs
%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, June / 2019

% Takes defined rectangular region from the center, finds the interfaces, saves the profiles. 

clear all, close all, clc;

% Which measurement
whichfile = 1;
% Not working: 2, 5, 9




% For saving
foldername = pwd;
foldername = foldername(max(strfind(foldername,'/'))+1:max(strfind(foldername,'/'))+4); %Checking the folder name

%If saving == 1, saves the output
saving = 0;

% Analysis
% Rectangle VOI - 1 mm3



% Loading the data  -----------------------------------------------------------------------------------
filename = dir('*RotatedVOI_data*.mat');
filename = filename(whichfile).name; %Reads the last mat file
load(filename);

names = DATA{1,1};

CECT_all = DATA{1,2};

Coordinates = DATA{1,3};

for location = 1:size(CECT_all,1) %How many measured locations
% location = 1;


counter = 1;

%How many measure dlocations
for measuredpoints = 1:size(CECT_all,1)
    
    for datalength = 1:2:length(names)
        names50{counter} = DATA{1,1}{1,datalength};
        CECT50{measuredpoints,counter} = DATA{1,2}{measuredpoints,datalength};
        
        names90{counter} = DATA{1,1}{1,1+datalength};
        CECT90{measuredpoints,counter} = DATA{1,2}{measuredpoints,1+datalength};
        
        counter = counter+1;
    end
    counter = 1;
end


% Displaying the images
% Can be driven after VOI_check_JM
timepoint = 1;


% % % % SEPARATE ENERGIES PLOT
% % % figure(1); 
% % % subplot(1,2,1)
% % % imagesc(squeeze(CECT50{location,timepoint}(:,floor(size(CECT50{location,timepoint},2)/2),:)), [-1500 5000]); %Can't flip because it'll mess up the indexes
% % % view(90, 90)
% % % axis equal
% % % title('50 kV')
% % % ylim([0 size(CECT50{location,timepoint},1)]);
% % % %caxis([0 5000])
% % % 
% % % subplot(1,2,2)
% % % imagesc(squeeze(CECT90{location,timepoint}(:,floor(size(CECT90{location,timepoint},2)/2),:)), [-1500 5000]);
% % % view(90, 90)
% % % axis equal
% % % title('90 kV')
% % % ylim([0 size(CECT90{location,timepoint},1)]);

%caxis([0 100])


% % % CHECK THESE! % % % %
resolution = 40;
windowsize = 1000/40; %um/resolution. Default 1mm, typically 1040um because of uneven 
% % % % % % % % % % % % % % 


% LIMITS  -----------------------------------------------------------------------------------------------------------------
%X and Y just from the middle
x1 = floor(size(CECT50{location,timepoint},1)/2)-ceil(windowsize/2); %Moves the window left if windowsize is odd number
x2 = floor(size(CECT50{location,timepoint},1)/2)+floor(windowsize/2);

y1 = floor(size(CECT50{location,timepoint},2)/2)-ceil(windowsize/2); %Moves the window left if windowsize is odd number
y2 = floor(size(CECT50{location,timepoint},2)/2)+floor(windowsize/2);




% ----------------------------------------------------------------------------------


clearvars profile50_x profile90_x
if ishandle(98)
close(98)
end
figure(98);


% Going through all the timepoints -------------
for i = 1:8
keke = CECT50{location,i};
keke(keke==-1000) = nan;
keke(keke==0) = nan;

profile50_x(:,i) = reshape(nanmean(keke(x1:x2,y1:y2,:),[1 2]),[],1); % Here (x1:x2,y1:y2,:) makes the code inspect only the center


plot(profile50_x)
hold on;

end
title('50 kV profiles');

if ishandle(99)
close(99)
end
figure(99);
for i = 1:8
keke = CECT90{location,i};
keke(keke==-1000) = nan;
keke(keke==0) = nan;

profile90_x(:,i) = reshape(nanmean(keke(x1:x2,y1:y2,:),[1 2]),[],1);

plot(profile90_x)
hold on;

end
title('90 kV profiles');

% FINDING THE LAYERS ------------------------------------------------------------------------------------------

% You can find the cartilage from the peaks of ka50 (first, last). Two points should be the same between 50 and 90

sumfor50 = sum(profile50_x,2);
sumfor90 = sum(profile90_x,2);


difference = diff(sumfor50);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% UNCOMMENT IF YOU WANT TO SEE THE PROFILE OR DIFFERENCE
% figure; plot(sumfor50)
% figure; plot(difference)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% [pks, locs] = findpeaks(difference(1:floor(length(difference)./2)),'NPeaks',1,'MinPeakWidth',3); %MinPeakWidth can exclude or
% includewron peaks. Using xax instead
[pks_temp1, locs_temp1] = findpeaks(difference(1:floor(length(difference)./2))); 
[pks, in] = max(pks_temp1); %finding the max
locs = locs_temp1(in); %Finding the location

% [pks(2), locs(2)] = findpeaks(difference(floor(length(difference)./2)+1:end),'NPeaks',1,'MinPeakWidth',2);
% locs(2) = locs(2) + floor(length(difference)./2)+1; %Correcting the index (missing the beginning)

[pks_temp2, locs_temp2] = findpeaks(difference(floor(length(difference)./2)+1:end));
[pks(2), in(2)] = max(pks_temp2); %finding the max
locs(2) = locs_temp2(in(2)) + floor(length(difference)./2); %Finding the location (and adding the missing half to the index)



% Renaming 
% Using the names z1,z2,x1,x2,y1,y2 to display the edges of ROI from now on. 
z1 = locs(1);
z2 = locs(2);



% disp(['Bone interface signal at ', num2str(location), ' point = ', num2str(pks(2))]);

% Displaying the lines in the image -----------------------------------------------------------------------------------------

figure(98);
line([locs(1) locs(1)], [min(profile50_x(:)) max(profile50_x(:))],'color','g');
line([locs(2) locs(2)], [min(profile50_x(:)) max(profile50_x(:))],'color','r');

figure(99);
line([locs(1) locs(1)], [min(profile50_x(:)) max(profile90_x(:))],'color','g');
line([locs(2) locs(2)], [min(profile50_x(:)) max(profile90_x(:))],'color','r');

% % IF YOU WANT TO SEE THE LIMITS ALSO WITH 90 kV
% figure(1)
% subplot(1,2,1)
% line([locs(1) locs(1)],[0, size(CECT50{location,timepoint},1)], 'Color', 'g')
% line([locs(2) locs(2)],[0, size(CECT50{location,timepoint},1)], 'Color', 'r')
% 
% subplot(1,2,2)
% line([locs(1) locs(1)],[0, size(CECT50{location,timepoint},1)], 'Color', 'g')
% line([locs(2) locs(2)],[0, size(CECT50{location,timepoint},1)], 'Color', 'r')



% Checking also the 2nd angle (y) ------------------------------------------------------------------------------------------

figure(2); 
subplot(1,2,1)
imagesc(squeeze(CECT50{location,timepoint}(:,floor(size(CECT50{location,timepoint},2)/2),:)), [-1500 5000]); %Can't flip because it'll mess up the indexes
view(90, 90)
axis equal
title('X-angle')
ylim([0 size(CECT50{location,timepoint},1)]);
%caxis([0 5000])

subplot(1,2,2)
imagesc(squeeze(CECT50{location,timepoint}(floor(size(CECT50{location,timepoint},1)/2),:,:)), [-1500 5000]); %Can't flip because it'll mess up the indexes
view(90, 90)
axis equal
title('Y-angle')
ylim([0 size(CECT50{location,timepoint},2)]);


% Cropping the 1mm rectangle region  -----------------------------------------------------------------------------------------------------
% Not doing this manually



% Taking a 1 mm rectangle inside the VOI
figure(2)

%Vertical
subplot(1,2,1)
line([locs(1) locs(1)],[0, size(CECT50{location,timepoint},1)], 'Color', 'g')
line([locs(2) locs(2)],[0, size(CECT50{location,timepoint},1)], 'Color', 'r')
subplot(1,2,2)
line([locs(1) locs(1)],[0, size(CECT50{location,timepoint},2)], 'Color', 'g')
line([locs(2) locs(2)],[0, size(CECT50{location,timepoint},2)], 'Color', 'r')

%Horizontal
subplot(1,2,1)
line([0, size(CECT50{location,timepoint},3)],[x1, x1], 'Color', 'g')
line([0, size(CECT50{location,timepoint},3)],[x2, x2], 'Color', 'r')
subplot(1,2,2)
line([0, size(CECT50{location,timepoint},3)],[y1, y1], 'Color', 'g')
line([0, size(CECT50{location,timepoint},3)],[y2, y2], 'Color', 'r')





% Saving the profiles

% Analyzing again all the timepoints
RESULT_PROFILES50{location} = profile50_x(z1:z2,:);
RESULT_PROFILES90{location} = profile90_x(z1:z2,:);


pause(1)
end


% Save, but don't overwrite
if saving == 1;
    savename = filename(1:end-14)
save([savename, '_', 'RESULT_PROFILES', num2str(length(dir('*RESULT_PROFILES*.mat'))+1), '.mat'],'RESULT_PROFILES50', 'RESULT_PROFILES90')
end

end % function