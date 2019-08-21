% function DATA = CECT2VOI_JM
%% m-file for analysing CECT image VOIs
%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, June / 2019

clear all,%close all, clc;


% Analysis
% Rectangle VOI - 1 mm3
% for mittauspiste = 1:6
mittauspiste = 1;



% Loading the data  -----------------------------------------------------------------------------------
filename = dir('kuusi*VOI_data*.mat');
filename = filename(end).name; %Reads the last mat file
load(filename);

names = DATA{1,1};

CECT_all = DATA{1,2};

Coordinates = DATA{1,3};


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

figure(1); 
subplot(1,2,1)
imagesc(squeeze(CECT50{mittauspiste,timepoint}(:,floor(size(CECT50{mittauspiste,timepoint},2)/2),:)), [-1500 5000]); %Can't flip because it'll mess up the indexes
view(90, 90)
axis equal
ylim([0 size(CECT50{mittauspiste,timepoint},1)]);
%caxis([0 5000])

subplot(1,2,2)
imagesc(squeeze(CECT90{mittauspiste,timepoint}(:,floor(size(CECT90{mittauspiste,timepoint},2)/2),:)), [-1500 5000]);
view(90, 90)
axis equal
ylim([0 size(CECT90{mittauspiste,timepoint},1)]);

%caxis([0 100])

% Cropping the region

% region_z = sum(

% ----------------------------------------------------------------------------------


clearvars profile50_x profile90_x
if ishandle(98)
close(98)
end
figure(98);
for i = 1:8
keke = CECT50{mittauspiste,i};
keke(keke==-1000) = nan;
keke(keke==0) = nan;

profile50_x(:,i) = reshape(nanmean(keke,[1 2]),[],1);
% profile50_y = mean(keke,[2 3]);


plot(profile50_x)
hold on;
% plot(profile50_y)
%ylim([0 2500])

end
title('50 kV profiles');

if ishandle(99)
close(99)
end
figure(99);
for i = 1:8
keke = CECT90{mittauspiste,i};
keke(keke==-1000) = nan;
keke(keke==0) = nan;

profile90_x(:,i) = reshape(nanmean(keke,[1 2]),[],1);
% profile50_y = mean(keke,[2 3]);


plot(profile90_x)
hold on;
% plot(profile50_y)
%ylim([0 2500])

end
title('90 kV profiles');

% FINDING THE LAYERS ------------------------------------------------------------------------------------------

% You can find the cartilage from the peaks of ka50 (first, last). Two points should be the same between 50 and 90

sumfor50 = sum(profile50_x,2);
sumfor90 = sum(profile90_x,2);

% figure; plot(sumfor50)

difference = diff(sumfor50);

% figure; plot(difference)

[pks, locs] = findpeaks(difference(1:floor(length(difference)./2)),'NPeaks',1,'MinPeakWidth',2); 
% Can't take the first as it might be the background interface



[pks(2), locs(2)] = findpeaks(difference(floor(length(difference)./2)+1:end),'NPeaks',1,'MinPeakWidth',2);
locs(2) = locs(2) + floor(length(difference)./2)+1; %Correcting the index (missing the beginning)

disp(['Bone interface signal at ', num2str(mittauspiste), ' point = ', num2str(pks(2))]);



% Displaying the lines in the image -----------------------------------------------------------------------------------------

figure(98);
line([locs(1) locs(1)], [min(profile50_x(:)) max(profile50_x(:))],'color','g');
line([locs(2) locs(2)], [min(profile50_x(:)) max(profile50_x(:))],'color','r');

figure(99);
line([locs(1) locs(1)], [min(profile50_x(:)) max(profile90_x(:))],'color','g');
line([locs(2) locs(2)], [min(profile50_x(:)) max(profile90_x(:))],'color','r');


figure(1)
subplot(1,2,1)
line([locs(1) locs(1)],[0, 100], 'Color', 'g')
line([locs(2) locs(2)],[0, 100], 'Color', 'r')

subplot(1,2,2)
line([locs(1) locs(1)],[0, 100], 'Color', 'g')
line([locs(2) locs(2)],[0, 100], 'Color', 'r')


pause(1);

% Choosing the window  -----------------------------------------------------------------------------------------------------
% Not doing this manually
resolution = 40;
windowsize = 2000/40;















% end

% Cropping the image ( 1 mm ) rectangle



