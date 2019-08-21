% For checking the created volumes of interest VOI

% clear all, close all, clc; 
close all

%cd /media/janne/Elements/uCT registrations/2Ri3

% This just for checking the registration

filename = dir('*_VOI_data*.mat');
filename = filename(end).name; %Reads the last mat file
DATA = load(filename);


h2 = waitbar(0,'Loading the files, please wait...'); %Display waitbar

whichpoint = 1; %Determines the point of measurement

for ii = [1:2:length(DATA.DATA{1,2}), 2:2:length(DATA.DATA{1,2})]

% figure;
% imagesc(DATA.DATA{1,2}{whichpoint,ii}(:,:,floor(size(DATA.DATA{1,2}{whichpoint,ii},3)/2))); %Showing mid slice
    waitbar(ii/length(DATA.DATA{1,2}));



% dicom_slider(DATA.DATA{1,2}{1,6})


%  Turning
% whichtime = 6;

Dicoms = DATA.DATA{1,2}{whichpoint,ii};
%Preallocating for efficiancy
SUBIM_x = ones(size(Dicoms,3), size(Dicoms,1), size(Dicoms,2));
SUBIM_y = ones(size(Dicoms,3), size(Dicoms,2), size(Dicoms,1));

% h = waitbar(0,'X-direction, please wait...'); %Display waitbar

%Showing projections
for i = 1:size(Dicoms,2)
    for j = 1:size(Dicoms,3)
        SUBIM_x(j,:,i) = Dicoms(:,i,j);
    end
%     waitbar(i/size(Dicoms,2));
end

% close(h)


% h = waitbar(0,'Y-direction, please wait...'); %Display waitbar

for i = 1:size(Dicoms,2)
    for j = 1:size(Dicoms,3)
        SUBIM_y(j,i,:) = Dicoms(:,i,j);
    end
%     waitbar(i/size(Dicoms,2));
end

% close(h)

figure; 
subplot(1,2,1)
imagesc(SUBIM_x(:,:,floor(size(SUBIM_x,3)/2)), [-1500 5000]);
title(DATA.DATA{1,1}{1,ii},'interpreter', 'none')
% caxis([0 5000])
subplot(1,2,2)
imagesc(SUBIM_y(:,:,floor(size(SUBIM_y,3)/2)), [-1500 5000]);
% caxis([0 5000])

pause(0.5)
end
close(h2)
