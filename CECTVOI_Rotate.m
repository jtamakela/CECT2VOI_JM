function [RESULT_PROFILES50, RESULT_PROFILES90] = CECTVOI_Analysis


% This could be just a a standalone that rotates matrixes

%% m-file for analysing CECT image VOIs
%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, June / 2019

% Takes defined rectangular region from the center, finds the interfaces, saves the profiles.

clear all, close all, clc;

% % % % % % % % % % % % For saving
% % % % % % % % % % % foldername = pwd;
% % % % % % % % % % % foldername = foldername(max(strfind(foldername,'/'))+1:max(strfind(foldername,'/'))+4); %Checking the folder name
% % % % % % % % % % % 
% % % % % % % % % % % %If saving == 1, saves the output
% % % % % % % % % % % saving = 0;

% Analysis
% Rectangle VOI - 1 mm3

% Loading the data  -----------------------------------------------------------------------------------
filename = dir('*VOI_data*.mat');

whichfile = 7; %Reads the numbered last mat file
whichlocation = 3; %There should be a six of these for the ponies


filename = filename(whichfile).name; 
load(filename);

names = DATA{1,1};

CECT_all = DATA{1,2};

Coordinates = DATA{1,3};

% ----------------------------------------------------------------------------------------------------


% Choosing the location and the sample  -----------------------------------------------------------------------------------

% % % % % % % % % % for location = 1:size(CECT_all,1) %How many measured locations
    
    %6,2 on hyvä vino esimerkki
    %6,1
    
    counter = 1;
    
    
% ----------------------------------------------------------------------------------------------------


% Importing the data ------------------------------------------------------------------------------------------------------
    %How many measured locations
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
    
% ----------------------------------------------------------------------------------------------------

    % Displaying the images
    % Which timepoint is used (default 1 = baseline)
    timepoint = 1;
    
    
    
% Finding the orientation -----------------------------------------------------------------------------------------------
fit_coefficients = find_demarcation(CECT90{whichlocation,timepoint})
    

fit_coefficients = find_demarcation(CECT90{whichlocation,timepoint})
% -----------------------------------------------------------------------------------------------------------------------



    
    
    
    
    
    
    %%
    
% % % %     % SEPARATE ENERGIES PLOT
% % % %     figure(1);
% % % %     subplot(1,2,1)
% % % %     %Taking the halfway slice
% % % %     imagesc(squeeze(CECT50{location,timepoint}(:,floor(size(CECT50{location,timepoint},2)/2),:)), [-1500 5000]); %Can't flip because it'll mess up the indexes
% % % %     view(90, 90)
% % % %     axis equal
% % % %     title('50 kV')
% % % %     ylim([0 size(CECT50{location,timepoint},1)]);
% % % %     caxis([0 5000])
% % % %     
% % % %     subplot(1,2,2)
% % % %     imagesc(squeeze(CECT90{location,timepoint}(:,floor(size(CECT90{location,timepoint},2)/2),:)), [-1500 5000]);
% % % %     view(90, 90)
% % % %     axis equal
% % % %     title('90 kV')
% % % %     ylim([0 size(CECT90{location,timepoint},1)]);
% % % %     
% % % %     caxis([0 1000])
    
    
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
    
    % ---------------------------------------------------------------------------------------------------------------------------------
    % Here we rotate the stack if needed
    
% % % % % % %     Let's go first with the 50 kV % % % % % % % % % % % % % % 
    
    % Creates point-of-views from both directions
    SUBIM = CECT50{location,timepoint};
    [dicom_swmask, dicom_swmask_y] = maskcreator(SUBIM);
    
    
    slider_question = menu('Does the figure need to be aligned:','1) Yes','2) No');
    
    while slider_question < 3
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        %Creates a image stack from x-angle (ROT_SUBIM) and if needed, orientates the image
        % First, checking the angle
        [ROT_SUBIM, direction_question, angle] = orientation(SUBIM, dicom_swmask, dicom_swmask_y,slider_question);
        
        
       %% 
        
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - 
%         Tää ei nyt oikein toimi % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%         ELI KÄÄNNÖN JÄLKEEN % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - 

        keyboard
        %Then the rotation
        if direction_question == 1
            ROT_SUBIM = imrotate(SUBIM_x,angle); %Actual imagestack
        else
            %This is only from the y-direction
            ROT_SUBIM_Y = imrotate(SUBIM_y,angle); %Actual imagestack
            %             dicom_slider(ROT_SUBIM_Y,98) %Using dicom_slider.m function for viewing
            
            figure(3);
            %imshow(sneakpeak);
            imagesc(sneakpeak);
            axis equal;
            
            %Correcting the X-direction
            clear ROT_SUBIM %Needs to be cleared because of the new size
            for i = 1:size(ROT_SUBIM_Y,2)
                for j = 1:size(ROT_SUBIM_Y,3)
                    ROT_SUBIM(:,j,i) = ROT_SUBIM_Y(:,i,j);
                end
            end
        end
        
        
        
        
        
        
        
        if slider_question == 1 %Make a new SUBIM and ask for a new run
            
            %For re-orienting, create new original image stack SUBIM
            clear SUBIM
            h = waitbar(0,'Creating a new image stack, please wait...'); %Display waitbar
            for i = 1:size(ROT_SUBIM,3)
                for j = 1:size(ROT_SUBIM,1)
                    SUBIM(:,i,j) = ROT_SUBIM(j,:,i);
                end
                waitbar(i/size(ROT_SUBIM,3));
            end
            close(h)
            close figure 100;
            dicom_slider(SUBIM,100) %Using dicom_slider.m function for viewing
            
            %Make new masks
            [dicom_swmask, dicom_swmask_y] = maskcreator(SUBIM);
            
            new_question = menu('Does the figure still need to be aligned:','1) Yes','2) No');
        else
            slider_question = 3;
            new_question = 2;
        end
        
        if new_question == 2
            slider_question = 3; %moving on from while
        end
        
    end
    
    
    
    
    % ---------------------------------------------------------------------------------------------------------------------------------
    
    
    
    % ----------------------------------------------------------------------------------
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     clearvars profile50_x profile90_x
% % % % % % % % % %     if ishandle(98)
% % % % % % % % % %         close(98)
% % % % % % % % % %     end
% % % % % % % % % %     figure(98);
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     % Going through all the timepoints -------------
% % % % % % % % % %     for i = 1:8
% % % % % % % % % %         keke = CECT50{location,i};
% % % % % % % % % %         
% % % % % % % % % %         keke(keke==-1000) = nan;
% % % % % % % % % %         keke(keke==0) = nan;
% % % % % % % % % %         
% % % % % % % % % %         profile50_x(:,i) = reshape(nanmean(keke(x1:x2,y1:y2,:),[1 2]),[],1); % Here (x1:x2,y1:y2,:) makes the code inspect only the center
% % % % % % % % % %         
% % % % % % % % % %         
% % % % % % % % % %         plot(profile50_x)
% % % % % % % % % %         hold on;
% % % % % % % % % %         
% % % % % % % % % %     end
% % % % % % % % % %     title('50 kV profiles');
% % % % % % % % % %     
% % % % % % % % % %     if ishandle(99)
% % % % % % % % % %         close(99)
% % % % % % % % % %     end
% % % % % % % % % %     figure(99);
% % % % % % % % % %     for i = 1:8
% % % % % % % % % %         keke = CECT90{location,i};
% % % % % % % % % %         keke(keke==-1000) = nan;
% % % % % % % % % %         keke(keke==0) = nan;
% % % % % % % % % %         
% % % % % % % % % %         profile90_x(:,i) = reshape(nanmean(keke(x1:x2,y1:y2,:),[1 2]),[],1);
% % % % % % % % % %         
% % % % % % % % % %         plot(profile90_x)
% % % % % % % % % %         hold on;
% % % % % % % % % %         
% % % % % % % % % %     end
% % % % % % % % % %     title('90 kV profiles');
% % % % % % % % % %     
% % % % % % % % % %     % FINDING THE LAYERS ------------------------------------------------------------------------------------------
% % % % % % % % % %     
% % % % % % % % % %     % You can find the cartilage from the peaks of ka50 (first, last). Two points should be the same between 50 and 90
% % % % % % % % % %     
% % % % % % % % % %     sumfor50 = sum(profile50_x,2);
% % % % % % % % % %     sumfor90 = sum(profile90_x,2);
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     difference = diff(sumfor50);
% % % % % % % % % %     
% % % % % % % % % %     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % %     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % %     % UNCOMMENT IF YOU WANT TO SEE THE PROFILE OR DIFFERENCE
% % % % % % % % % %     % figure; plot(sumfor50)
% % % % % % % % % %     % figure; plot(difference)
% % % % % % % % % %     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % %     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % %     
% % % % % % % % % %     % [pks, locs] = findpeaks(difference(1:floor(length(difference)./2)),'NPeaks',1,'MinPeakWidth',3); %MinPeakWidth can exclude or
% % % % % % % % % %     % includewron peaks. Using xax instead
% % % % % % % % % %     [pks_temp1, locs_temp1] = findpeaks(difference(1:floor(length(difference)./2)));
% % % % % % % % % %     [pks, in] = max(pks_temp1); %finding the max
% % % % % % % % % %     locs = locs_temp1(in); %Finding the location
% % % % % % % % % %     
% % % % % % % % % %     % [pks(2), locs(2)] = findpeaks(difference(floor(length(difference)./2)+1:end),'NPeaks',1,'MinPeakWidth',2);
% % % % % % % % % %     % locs(2) = locs(2) + floor(length(difference)./2)+1; %Correcting the index (missing the beginning)
% % % % % % % % % %     
% % % % % % % % % %     [pks_temp2, locs_temp2] = findpeaks(difference(floor(length(difference)./2)+1:end));
% % % % % % % % % %     [pks(2), in(2)] = max(pks_temp2); %finding the max
% % % % % % % % % %     locs(2) = locs_temp2(in(2)) + floor(length(difference)./2); %Finding the location (and adding the missing half to the index)
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     % Renaming
% % % % % % % % % %     % Using the names z1,z2,x1,x2,y1,y2 to display the edges of ROI from now on.
% % % % % % % % % %     z1 = locs(1);
% % % % % % % % % %     z2 = locs(2);
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     % disp(['Bone interface signal at ', num2str(location), ' point = ', num2str(pks(2))]);
% % % % % % % % % %     
% % % % % % % % % %     % Displaying the lines in the image -----------------------------------------------------------------------------------------
% % % % % % % % % %     
% % % % % % % % % %     figure(98);
% % % % % % % % % %     line([locs(1) locs(1)], [min(profile50_x(:)) max(profile50_x(:))],'color','g');
% % % % % % % % % %     line([locs(2) locs(2)], [min(profile50_x(:)) max(profile50_x(:))],'color','r');
% % % % % % % % % %     
% % % % % % % % % %     figure(99);
% % % % % % % % % %     line([locs(1) locs(1)], [min(profile50_x(:)) max(profile90_x(:))],'color','g');
% % % % % % % % % %     line([locs(2) locs(2)], [min(profile50_x(:)) max(profile90_x(:))],'color','r');
% % % % % % % % % %     
% % % % % % % % % %     % % IF YOU WANT TO SEE THE LIMITS ALSO WITH 90 kV
% % % % % % % % % %     % figure(1)
% % % % % % % % % %     % subplot(1,2,1)
% % % % % % % % % %     % line([locs(1) locs(1)],[0, size(CECT50{location,timepoint},1)], 'Color', 'g')
% % % % % % % % % %     % line([locs(2) locs(2)],[0, size(CECT50{location,timepoint},1)], 'Color', 'r')
% % % % % % % % % %     %
% % % % % % % % % %     % subplot(1,2,2)
% % % % % % % % % %     % line([locs(1) locs(1)],[0, size(CECT50{location,timepoint},1)], 'Color', 'g')
% % % % % % % % % %     % line([locs(2) locs(2)],[0, size(CECT50{location,timepoint},1)], 'Color', 'r')
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     % Checking also the 2nd angle (y) ------------------------------------------------------------------------------------------
% % % % % % % % % %     
% % % % % % % % % %     figure(2);
% % % % % % % % % %     subplot(1,2,1)
% % % % % % % % % %     imagesc(squeeze(CECT50{location,timepoint}(:,floor(size(CECT50{location,timepoint},2)/2),:)), [-1500 5000]); %Can't flip because it'll mess up the indexes
% % % % % % % % % %     view(90, 90)
% % % % % % % % % %     axis equal
% % % % % % % % % %     title('X-angle')
% % % % % % % % % %     ylim([0 size(CECT50{location,timepoint},1)]);
% % % % % % % % % %     %caxis([0 5000])
% % % % % % % % % %     
% % % % % % % % % %     subplot(1,2,2)
% % % % % % % % % %     imagesc(squeeze(CECT50{location,timepoint}(floor(size(CECT50{location,timepoint},1)/2),:,:)), [-1500 5000]); %Can't flip because it'll mess up the indexes
% % % % % % % % % %     view(90, 90)
% % % % % % % % % %     axis equal
% % % % % % % % % %     title('Y-angle')
% % % % % % % % % %     ylim([0 size(CECT50{location,timepoint},2)]);
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     % Cropping the 1mm rectangle region  -----------------------------------------------------------------------------------------------------
% % % % % % % % % %     % Not doing this manually
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     % Taking a 1 mm rectangle inside the VOI
% % % % % % % % % %     figure(2)
% % % % % % % % % %     
% % % % % % % % % %     %Vertical
% % % % % % % % % %     subplot(1,2,1)
% % % % % % % % % %     line([locs(1) locs(1)],[0, size(CECT50{location,timepoint},1)], 'Color', 'g')
% % % % % % % % % %     line([locs(2) locs(2)],[0, size(CECT50{location,timepoint},1)], 'Color', 'r')
% % % % % % % % % %     subplot(1,2,2)
% % % % % % % % % %     line([locs(1) locs(1)],[0, size(CECT50{location,timepoint},2)], 'Color', 'g')
% % % % % % % % % %     line([locs(2) locs(2)],[0, size(CECT50{location,timepoint},2)], 'Color', 'r')
% % % % % % % % % %     
% % % % % % % % % %     %Horizontal
% % % % % % % % % %     subplot(1,2,1)
% % % % % % % % % %     line([0, size(CECT50{location,timepoint},3)],[x1, x1], 'Color', 'g')
% % % % % % % % % %     line([0, size(CECT50{location,timepoint},3)],[x2, x2], 'Color', 'r')
% % % % % % % % % %     subplot(1,2,2)
% % % % % % % % % %     line([0, size(CECT50{location,timepoint},3)],[y1, y1], 'Color', 'g')
% % % % % % % % % %     line([0, size(CECT50{location,timepoint},3)],[y2, y2], 'Color', 'r')
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     % Saving the profiles
% % % % % % % % % %     
% % % % % % % % % %     % Analyzing again all the timepoints
% % % % % % % % % %     RESULT_PROFILES50{location} = profile50_x(z1:z2,:);
% % % % % % % % % %     RESULT_PROFILES90{location} = profile90_x(z1:z2,:);
% % % % % % % % % %     
% % % % % % % % % %     
% % % % % % % % % %     pause(1)
% % % % % % % % % % end
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % % % Save, but don't overwrite
% % % % % % % % % % if saving == 1;
% % % % % % % % % %     save([foldername, '_', 'RESULT_PROFILES', num2str(length(dir('*RESULT_PROFILES*.mat'))+1), '.mat'],'RESULT_PROFILES50', 'RESULT_PROFILES90')
% % % % % % % % % % end

end % function

% To automatically orient
function treshold_fit = find_orientation(thestack)

    % Let's check if the figure needs to be automatically aligned
    % 90 keV data goes great for trying to find bone cartilage interface
    orientationtesting = thestack;
    
    % ############################################################################################################################
    %The following needs to be in its own function
    
    figure(2);
    
    
    for which_slice= 30;%1:floor(size(orientationtesting,2))
        midslice = squeeze(orientationtesting(:,which_slice,:));
        %
        %         imagesc(midslice);
        %         view(90, 90)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fororientation = midslice>400; %THIS THRESHOLD HAS NOT BEEN OPTIMIZED
        % Though doesn't matter if not excatly on point.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Just taking the bone-cartilage interface
        interface = zeros(1,size(fororientation,1)); %Preallocating
        for orientation_i = 1:size(fororientation,1) %Going through rows
            if find(fororientation(orientation_i,:) == 1) %If not empty
                interface(orientation_i) = min(find(fororientation(orientation_i,:) == 1)); %Just the interface
            end
        end
        
        imagesc(midslice)
        view(90, 90)
        hold on;
        % Showing the interface with red crosses
        plot(interface, [1:length(interface)], 'rx');
        pause(0.1);
        
        %Not taking zeroes into account
        % Make this into an own function
        function_for_fit = interface;
        
        %Won't do any fitting if empty
%         if function_for_fit %Not empty
            removable = (function_for_fit == 0);
            function_for_fit(removable) = [];
            points = 1:length(interface);
            %removing zeroes
            points(removable) = [];
            
            plot(function_for_fit, points, 'go');
            
            
            %     figure;
            if ~isempty(points)
            treshold_fit = fit(points',function_for_fit','poly1');
            hold on;
            plot([0:100].*treshold_fit.p1+treshold_fit.p2,[0:100], 'g--', 'linewidth' ,2)
            %     view(90, 90) %If you want to make a separate figure
            %     axis ij
                    pause(0.1);

        end
        
        
    end

    
    % ############################################################################################################################
end





function [dicom_swmask, dicom_swmask_y] = maskcreator(SUBIM)

clear dicom_swmask dicom_swmask_y

h = waitbar(0,'Creating Masks, please wait...'); %Display waitbar
for i = 1:size(SUBIM,1) %y-direction
    for j = 1:size(SUBIM,3) %x-direction
        %Commented can be used to create black and white figures
        % % %         if find(SUBIM(i,:,j) > 0) %Look for values above background
        % % %             dicom_swmask(j,i,:) = 1;
        % % %         else
        % % %             dicom_swmask(j,i,:) = 0;
        % % %         end
        %       Or just take mean
        dicom_swmask(j,i,:) = mean(SUBIM(i,:,j));
    end
end

for i = 1:size(SUBIM,2) %y-direction
    for j = 1:size(SUBIM,3) %x-direction
        %Commented can be used to create black and white figures
        % % %         if find(SUBIM(:,i,j) > 0) %Look for values above background
        % % %             dicom_swmask_y(j,i,:) = 1;
        % % %         else
        % % %             dicom_swmask_y(j,i,:) = 0;
        % % %         end
        dicom_swmask_y(j,i,:) = mean(SUBIM(:,i,j));
    end
    waitbar(i/(size(SUBIM,1)));
end
close(h);

% figure(2); %Mask image
% imshow(dicom_swmask)
figure(200)
clf(figure(200))
set(200,'position', [400 200 1200 500]);
subplot(1,2,1);
% imshow(dicom_swmask);
imagesc(dicom_swmask);
axis equal;
title('From x-direction (Scanner door?)'); % Haven't confirmed this
subplot(1,2,2);
% imshow(dicom_swmask_y);
imagesc(dicom_swmask_y);
axis equal;
title('From y-direction');
end

%%

function [ROT_SUBIM, direction_question, rotangle] = orientation(SUBIM, dicom_swmask, dicom_swmask_y, slider_question_fororientation)
% Rotating the imagestack  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(200); %Mask image

%No orientation done
if slider_question_fororientation == 2
    
    for i = 1:size(SUBIM,2)
        for j = 1:size(SUBIM,3)
            SUBIM_x(j,:,i) = SUBIM(:,i,j);
        end
    end
    
    ROT_SUBIM = SUBIM_x;
    rotangle = 0;
    dicom_slider(ROT_SUBIM,98)
end

%Orientation is done
if slider_question_fororientation == 1
    
    
    q = 2; % For the while command
    direction_question = menu('Which Direction:','1) X','2) Y');
    
    
    %Rotating the figure
    while q == 2
        
        % EDIT X-DIRECTION % % % % % % % % % % % % % % % % % % % % % % % % %
        if direction_question == 1
            figure(200); %Mask image
            subplot(1,2,1);
            pause(0.6); %Gets rid off random freezing in linux
            title('Please pick using your mouse two points in the center, left and right');
            
            %imshow(dicom_swmask)
            imagesc(dicom_swmask)
            axis equal;
            hold on;
            
            [xrot, yrot] = ginput(2);
            line(xrot,yrot)
            rotangle = atan( (yrot(2)-yrot(1)) / (xrot(2)-xrot(1)) ) * (180/pi); %The angle in degrees
            
            sneakpeak = imrotate(dicom_swmask,rotangle); %Display the mask
            %Draw lines for comparison
            for i = 1:10:size(sneakpeak,2)
                sneakpeak(i,:) = 0;
            end
            
% % %             ROT_SUBIM = imrotate(SUBIM_x,rotangle); %Actual imagestack
%             dicom_slider(ROT_SUBIM,98) %Using dicom_slider.m function for viewing
            
            figure(3);
            %imshow(sneakpeak);
            imagesc(sneakpeak);
            axis equal;
            
        end
        
        
        
        % EDIT Y-DIRECTION % % % % % % % % % % % % % % % % % % % % % % % % %
        if direction_question == 2
            
            
            figure(200); %Mask image
            subplot(1,2,2);
            pause(0.6); %Gets rid off random freezing in linux
            title('Please pick using your mouse two points in the center, left and right');
            
            %imshow(dicom_swmask_y)
            imagesc(dicom_swmask_y)
            axis equal;
            hold on;
            
            [xrot, yrot] = ginput(2);
            line(xrot,yrot)
            rotangle = atan( (yrot(2)-yrot(1)) / (xrot(2)-xrot(1)) ) * (180/pi); %The angle in degrees
            
            sneakpeak = imrotate(dicom_swmask_y,rotangle); %Display the mask
            %Draw lines for comparison
            for i = 1:10:size(sneakpeak,2)
                sneakpeak(i,:) = 0;
            end
            
            %This is only from the y-direction
% % %             ROT_SUBIM_Y = imrotate(SUBIM_y,rotangle); %Actual imagestack
%             dicom_slider(ROT_SUBIM_Y,98) %Using dicom_slider.m function for viewing
            
            figure(3);
            %imshow(sneakpeak);
            imagesc(sneakpeak);
            axis equal;
            
            %Correcting the X-direction
% % %             clear ROT_SUBIM %Needs to be cleared because of the new size
% % %             for i = 1:size(ROT_SUBIM_Y,2)
% % %                 for j = 1:size(ROT_SUBIM_Y,3)
% % %                     ROT_SUBIM(:,j,i) = ROT_SUBIM_Y(:,i,j);
% % %                 end
% % %             end
            
        end
        
        q = menu('Are you satisfied with the angle:','1) Yes','2) No');
        
    end
end
end

