function DATA = CECT2VOI_JM
%% m-file for analysing registered CECT images (.nii)
%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/CECT2VOI_JM
%% (c) Janne Mäkelä June / 2019
% Click on the measurement location and analyze
% Creates a cubical VOI from the chosen location which is then saved for later analysis


% SAVES 'FOLDERNAME'_VOI_DATA.MAT -FILE
% With names, matrixes, and coordinates for backup
% Coordinates [x1, y1, z1; x2, y2, z2];


clear all, close all, %clc;

% cd /media/janne/Data/UEF/Measurements/Ponies/'uCT registrations'/2Ri3/

foldername = pwd;
foldername = foldername(max(strfind(foldername,'/'))+1:max(strfind(foldername,'/'))+4); %Checking the folder name

filesavename = 'VOI_data.mat'; %Saved under this name




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% CHECK THIS ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
% RESOLUTION OF THE CT STACK
resolution = [40 40 40]; %[Z X Y]. Voxel size in micrometers. Defines also the aspect ratio in figures
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %








aspectratio = resolution./min(resolution); %Drawing the figures based on the given resolution

%Preallocating the final parameters
% Thicknesses = cell([]);
info = [];

filuname = dir('0h_50registration.nii');
filuname = filuname.name;
% LOAD IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dicoms = niftiread(filuname);
info = niftiinfo(filuname);


% % This is for rescaling the voxel values
% % Not necessary if imagesc is used
%Dicoms = Dicoms.*info.RescaleSlope+info.RescaleIntercept;
% % Otherwise handles data using native pixel values (original, short integer value)


%Orienting the figures
[Dicoms_x, Dicoms_y] = orientation(Dicoms);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % Options to display all the slices using dicom_slider function
% dicom_slider(Dicoms,100)
% dicom_slider(Dicoms_x,100) %Using dicom_slider.m function for viewing
% dicom_slider(Dicoms_y,100) %Using dicom_slider.m function for viewing
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%Mean image for picking the measurement point
dicom_mask = mean(Dicoms,3);
figure(1);
colormap jet
imagesc(dicom_mask)
axis equal;
xlabel(['x']); ylabel(['y']);
hold on;


%Choose the location
figure(1);
pause(1) %Reduces crashing
title('Pick the location. Press enter to quit'); %

% question = menu('Satisfied?','1) Yes','2) No'); %Option to fall back should be added
question = 1; %"Satisfied?"

%Pick the first location
[xcoord, ycoord] = ginput(1);

% Preallocation
location_i = 1;
% Marking the locations ----------------------------------------------------------------
while ~isempty(xcoord) %If enter is not pressed
    
    plot(xcoord,ycoord,'+','markersize', 40, 'Linewidth', 1.5)
    text(xcoord+20,ycoord-20,num2str(location_i),'HorizontalAlignment','center','fontsize', 20);
    
    %Displaying the chosen location from two angles
    %Calculating the thickness
    
    
    slice_x = Dicoms_x(:,:,round(xcoord));
    slice_y = Dicoms_y(:,:, round(ycoord));
    
    Coordinates{location_i} = imdistancecalculator(slice_x,slice_y, xcoord, ycoord, aspectratio); %[Z X Y]
    
    if Coordinates{location_i} ~= [1, 2, 3, 4, 5, 6] %If user is not satisfied with the location, imdistancecalulator returns [1, 2, 3, 4, 5, 6];
        Thicknesses(location_i) = sqrt( (resolution(1)*(Coordinates{location_i}(4)-Coordinates{location_i}(1)))^2 + (resolution(2)*(Coordinates{location_i}(5)-Coordinates{location_i}(2)))^2 + (resolution(3)*(Coordinates{location_i}(6)-Coordinates{location_i}(3)))^2);
        
        disp(['Measured thickness in #', num2str(location_i), ' is ', num2str(Thicknesses(location_i)), ' um'])
        
        XCOORD(location_i) = xcoord; %Saving
        YCOORD(location_i) = ycoord;
        
        save('THICKNESS_temp.mat','Thicknesses', 'XCOORD', 'YCOORD') % In case the code crashes
        
        
        location_i = location_i+1;
    end
    
    figure(1);
    pause(1) %Reduces crashing
    title('Pick the location. Press enter to quit'); %
    [xcoord, ycoord] = ginput(1);
    
end


%Sorting -----------------------------------------------------------------------------------------------

fileorder = {'Baseline_50','Baseline_90','0h_50','0h_90', '30min_50','30min_90','1h_50','1h_90','2h_50','2h_90','6h_50','6h_90','10h_50','10h_90','23h_50','23h_90'};

niifiles = dir('*.nii');


for i = 1: length(niifiles)
    niiname{i} = niifiles(i).name;
end

order = ones(1,length(niifiles));
%Finding the string order
for i = 1: length(niifiles)
    if length(find(contains(niiname, fileorder{i}) == 1) == 1) == 1 %I know, this is awful. Checks if there is only one instance (length) of found strings
        order(i) = find(contains(niiname, fileorder{i}) == 1);
        % '0h_50','0h_90' need to be separately searched because the 10h strings include these
    else
        temptemp =  find(contains(niiname, fileorder{i}) == 1);
        order(i) = min(temptemp);
    end
end


% Creating masks -----------------------------------------------------------------------------------------------
% keyboard
    % VOIs are analyzed for all the locations length(Coordinates)   
    h2 = waitbar(0,'Creating masks, please wait...'); %Display waitbar
    for numofpoints = 1:length(Coordinates) %The number of locations (6)
        
        %The two points that define ROI
        POINT1 = Coordinates{numofpoints}(1:3); % Above [X Y Z]
        POINT2 = Coordinates{numofpoints}(4:6); % Below [X Y Z]
        
        %Cropping [x, y, z]
        [MASK{numofpoints} lims{numofpoints}] = VOI_returner(Dicoms, POINT1, POINT2, numofpoints);
            waitbar(numofpoints/length(Coordinates));

    end
close(h2)
% -----------------------------------------------------------------------------------------------


% ANALYZING ALL THE FILES

h2 = waitbar(0,'Analyzing all timepoints, please wait...'); %Display waitbar
counter = 1;
for numoffiles = order %The order of measurements (16)
    
%     niiname{numoffiles} = niifiles(numoffiles).name;
    % LOAD IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dicoms = niftiread(niiname{numoffiles});
    for numofpoints = 1:length(Coordinates)                                         % POINT1                           % POINT2
        ORIENTEDVOI{numofpoints, counter} = VOIOrienter(Dicoms, MASK{numofpoints}, Coordinates{numofpoints}(1:3), Coordinates{numofpoints}(1:3), lims{numofpoints}); %Does the masking and orientation
    end
    waitbar(counter/length(niifiles));
    counter = counter + 1;
end
close(h2)

DATA = {niiname(order), ORIENTEDVOI, Coordinates'};

save([foldername, '_', filesavename], 'DATA')


end


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% This function takes the coordinates and returns the VOI

function [MASK, lims] = VOI_returner(Dicoms, POINT1, POINT2, point_i)

%For subplots
dualA = [POINT1(1) POINT1(3); POINT1(2:3)]; % Above point X,Z;Y,Z
dualB = [POINT2(1) POINT2(3); POINT2(2:3)]; % Below point X,Z;Y,Z

% SUBIMAGES: 3D to 2D
% This can be left undone for the following timepoints
for subimages_i = 1:2 % 1=X-direction, 2=Y-direction
    
    A = dualA(subimages_i,:);
    B = dualB(subimages_i,:);
    
    % Calculating normals
    % Length half the length of the vector
    normal = 0.5*norm(A-B).*(null(A-B)');
    
    %Coordinates of the normals (takes with if into account the orientations)
    if normal(1) < 0
        normalUL(subimages_i,:) = A + normal; %Upper left
        normalUR(subimages_i,:) = A - normal; %Upper right
        normalBL(subimages_i,:) = B + normal; %Bottom left
        normalBR(subimages_i,:) = B - normal; %Bottom right
    else
        normalUL(subimages_i,:) = A - normal; %Upper left
        normalUR(subimages_i,:) = A + normal; %Upper right
        normalBL(subimages_i,:) = B - normal; %Bottom left
        normalBR(subimages_i,:) = B + normal; %Bottom right
    end
    
        figure(point_i+1);
        subplot(1,2,3-subimages_i) %Right first
        
        line([A(1),normalUL(subimages_i,1)],[A(2),normalUL(subimages_i,2)],'color','r','LineWidth',2)
        line([A(1),normalUR(subimages_i,1)],[A(2),normalUR(subimages_i,2)],'color','r','LineWidth',2)
        
        line([B(1),normalBL(subimages_i,1)],[B(2),normalBL(subimages_i,2)],'color','r','LineWidth',2)
        line([B(1),normalBR(subimages_i,1)],[B(2),normalBR(subimages_i,2)],'color','r','LineWidth',2)
        
        %Rectangle
        line([normalUL(subimages_i,1),normalBL(subimages_i,1)],[normalUL(subimages_i,2),normalBL(subimages_i,2)],'color','r','LineWidth',2)
        line([normalUR(subimages_i,1),normalBR(subimages_i,1)],[normalUR(subimages_i,2),normalBR(subimages_i,2)],'color','r','LineWidth',2)
end

%
% 2D back to 3D
% 4 points to define cuboid sides
%          Looking 1 subplot
% X Y Z in original figure 1
% FROM ABOVE
%Origo - top left
%          subfig2       subfig1
ROI1 = floor([normalUL(1,1) normalUL(2,1) (normalUL(1,2)-A(2)+normalUL(2,2)-A(2))+A(2)]);
%X - top right
ROI2 = floor([normalUR(1,1) normalUL(2,1) (normalUR(1,2)-A(2)+normalUL(2,2)-A(2))+A(2)]);
%Y - bottom left
ROI3 = floor([normalUL(1,1) normalUR(2,1) (normalUL(1,2)-A(2)+normalUR(2,2)-A(2))+A(2)]);
%Z - top left (below)
ROI4 = floor([normalBL(1,1) normalBL(2,1) (normalBL(1,2)-B(2)+normalBL(2,2)-B(2))+B(2)]);
%Bottom right
ROI5 = floor([normalUR(1,1) normalUR(2,1) (normalUR(1,2)-A(2)+normalUR(2,2)-A(2))+A(2)]);
%Bottom right (below)
ROI6 = floor([normalBR(1,1) normalBR(2,1) (normalBR(1,2)-B(2)+normalBR(2,2)-B(2))+B(2)]);
%Bottom left (below)
ROI7 = floor([normalBL(1,1) normalBR(2,1) (normalBL(1,2)-B(2)+normalBR(2,2)-B(2))+B(2)]);
%Top right (below)
ROI8 = floor([normalBR(1,1) normalBL(2,1) (normalBR(1,2)-B(2)+normalBL(2,2)-B(2))+B(2)]);

ROI = [ROI1; ROI2; ROI3; ROI4; ROI5; ROI6; ROI7; ROI8];


%Checking and plotting the locations
% Only for the first timepoint
    
    figure(1);
    for i = 1:3
        eval(['plot(ROI',num2str(i),'(1), ROI',num2str(i),'(2),''ko'',''markersize'', 6, ''Linewidth'', 3)']);
    end
    plot(ROI4(1), ROI4(2),'ko','markersize', 6, 'Linewidth', 1); %Lower point plotted thinner
    
    for i = 2:4
        eval(['line([ROI1(1), ROI',num2str(i),'(1)],[ROI1(2), ROI',num2str(i),'(2)])']);
    end
    
% DRAWING THE CUBICLE IN A SEPARATE FIGURE
% % % % % % Separate figure of the axii
% % % % % figure;
% % % % %
% % % % % ax = axes;
% % % % % plot3(ROI1(1),ROI1(2),ROI1(3), '+','markersize', 3, 'Linewidth', 1.5); %First corner point
% % % % % hold on;
% % % % %
% % % % % xlabel(['x']); ylabel(['y']); zlabel(['z']);
% % % % % set(ax,'ZDir','reverse')
% % % % % % set(ax,'YDir','reverse')
% % % % % % set(ax,'XDir','reverse')
% % % % % axis equal
% % % % % % xlim([0 512]); ylim([0 512]); zlim([0 512])
% % % % %
% % % % % for i = 1:8
% % % % %     for ii = 1:8
% % % % %     line([ROI(ii,1) ROI(i,1)], [ROI(ii,2) ROI(i,2)], [ROI(ii,3) ROI(i,3)])
% % % % %     end
% % % % % end
% % % % % legend({'cp'; 'X'; 'Y'; 'Z'});
% % % % %
% % % % % plot3(POINT1(1), POINT1(2), POINT1(3), 'rx', 'markersize', 18)
% % % % % plot3(POINT2(1), POINT2(2), POINT2(3), 'rx', 'markersize', 18)
% % % % % line([POINT1(1) POINT2(1)], [POINT1(2) POINT2(2)], [POINT1(3) POINT2(3)], 'color', 'r', 'linewidth', 3)

% -------------------------------------------------------------------------------------------


% Outer limits to disregard completely the values outside these
x_lims = [min(ROI(:,1)) max(ROI(:,1))];
y_lims = [min(ROI(:,2)) max(ROI(:,2))];
z_lims = [min(ROI(:,3)) max(ROI(:,3))];


mask = logical(zeros(size(Dicoms)));
% VOI = (zeros(size(Dicoms)));

xyz = [ROI; [0 0 0 ]];
len = length(xyz(:,1));

% % % h = waitbar(0,'Counting the voxels of the cubical VOI, Please wait...'); %Display waitbar
% Checking only the values close to the mask
for x=x_lims(1):x_lims(2)
    for y=y_lims(1):y_lims(2)
        for z=z_lims(1):z_lims(2)
            
            xyz(end,:) = [x y z];
            % Inside the cube?
            % https://se.mathworks.com/matlabcentral/answers/21429-determine-if-a-point-is-inside-a-cube
            K = convhulln(xyz);
            if ~any(any(K==len))
                mask(y,x,z) = 1; %I dunno why x and y are in inverse order?
                %             VOI(x,y,z) = Dicoms(x,y,z);
            end
        end
        
    end
    % % %     waitbar((x-x_lims(1))/(x_lims(2)-x_lims(1)))
    
end
% % % close(h)


%Cropping [x, y, z]
MASK = mask(y_lims(1):y_lims(2),x_lims(1):x_lims(2),z_lims(1):z_lims(2));
lims = [x_lims; y_lims; z_lims];

end

 % %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%% %% %%

function VOIRotated = VOIOrienter(Dicoms, MASK, POINT1, POINT2, lims)

%Calculating angles -----------------------------------------------------------------------------------------
%
V = [POINT1(1)-POINT2(1); POINT1(2)-POINT2(2); POINT1(3)-POINT2(3)]; %Use these, these are exact
% % Taking between ROI 1 and ROI4 because these have been oriented correclty earlier
% V = [ROI(1,1)-ROI(4,1); ROI(1,2)-ROI(4,2); ROI(1,3)-ROI(4,3)]

% WORKS ONLY IF POINT 1 IS ABOVE
ax = -90+180/pi*atan2(norm(cross(V,[1;0;0])),dot(V,[1;0;0]));
ay = 90-180/pi*atan2(norm(cross(V,[0;1;0])),dot(V,[0;1;0]));
az = 180/pi*atan2(norm(cross(V,[0;0;1])),dot(V,[0;0;1])); %Not used, Most likely wrong atm.
th = [ax ay az];

% The limits for cropping
x_lims = [lims(1,1), lims(1,2)];
y_lims = [lims(2,1), lims(2,2)];
z_lims = [lims(3,1), lims(3,2)];
%

VOIRotated = MASK.*Dicoms(y_lims(1):y_lims(2),x_lims(1):x_lims(2),z_lims(1):z_lims(2));


VOIRotated = imrotate3(VOIRotated,th(2),[1 0 0],'nearest','loose','FillValues',-1000); %Nearest neighbor changes histogram very little
VOIRotated = imrotate3(VOIRotated,th(1),[0 1 0],'nearest','loose','FillValues',-1000);

end %FUNCTION

% % % % % % %% tHIS CAN BE USED FOR VALIDATION
% % % % % %
% % % % % % clear VOIRotated_x
% % % % % % clear VOIRotated_y
% % % % % %
% % % % % % clear VOIRotatedslice_x
% % % % % % clear VOIRotatedslice_y
% % % % % %
% % % % % % % Displaying the VOI
% % % % % % VOIRotated_y = ones(size(VOIRotated,3), size(VOIRotated,2), size(VOIRotated,1));
% % % % % % VOIRotated_x = ones(size(VOIRotated,3), size(VOIRotated,1), size(VOIRotated,2));
% % % % % %
% % % % % % h = waitbar(0,'X-direction, please wait...'); %Display waitbar
% % % % % %
% % % % % % for i = 1:size(VOIRotated,2)
% % % % % %     for j = 1:size(VOIRotated,3)
% % % % % %         VOIRotated_x(j,:,i) = VOIRotated(:,i,j);
% % % % % %     end
% % % % % %     waitbar(i/size(VOIRotated,2));
% % % % % % end
% % % % % %
% % % % % % close(h)
% % % % % %
% % % % % %
% % % % % % h = waitbar(0,'Y-direction, please wait...'); %Display waitbar
% % % % % %
% % % % % % for i = 1:size(VOIRotated,2)
% % % % % %     for j = 1:size(VOIRotated,3)
% % % % % %         VOIRotated_y(j,i,:) = VOIRotated(:,i,j);
% % % % % %     end
% % % % % %     waitbar(i/size(VOIRotated,2));
% % % % % % end
% % % % % %
% % % % % % close(h)
% % % % % %
% % % % % %
% % % % % %
% % % % % % VOIRotatedslice_x = VOIRotated_x(:,:,floor(mean(size(VOIRotated_x,3))/2));
% % % % % % VOIRotatedslice_y = VOIRotated_y(:,:,floor(mean(size(VOIRotated_y,3))/2));
% % % % % %
% % % % % % figure;
% % % % % % subplot(1,2,1)
% % % % % % imagesc(VOIRotatedslice_x);
% % % % % % title('X ->'); %Left hand figure
% % % % % % axis equal
% % % % % % subplot(1,2,2)
% % % % % % imagesc(VOIRotatedslice_y);%Right hand figure
% % % % % % title('Y -^', 'Interpreter', 'none');
% % % % % % axis equal
% % % % % %
% % % % % %
% % % % % %
% % % % % % % figure;
% % % % % % % slice(double(VOIRotated),floor(mean(x_lims)),floor(mean(y_lims)),floor(mean(z_lims)));
% % % % % % % grid on, shading interp, colormap jet
% % % % % % % xlabel(['x']); ylabel(['y']); zlabel(['z']);
% % % % % %
% % % % % % % figure;
% % % % % % % dicom_slider(VOIRotated)
% % % % % %
% % % % % % % % %
% % % % % % % % %
% % % % % % % % %
% % % % % % % dicom_slider(B)
% % % % % %
% % % % % % figure;
% % % % % % histogram(VOIRotated)
% % % % % % hold on;
% % % % % % histogram(VOI);
% % % % % % legend({'VOIRotated with nearest neighbor interp.'; 'VOI'})


%%
% function [Dicoms, info] = load_dicoms()
% %Loading first the dicoms
%
% path = uigetdir; %Choose the folder where the DICOMS are
%
% f = filesep; %Checks what's the file separator for current operating system (windows,unix,linux)
%
% dicomnames = dir([num2str(path) f '*.dcm*']); %Read dicoms.
% disp(['Folder: ', dicomnames(1).folder]); %display folder
% %Dicom info
% info = dicominfo([num2str(path) f dicomnames(1).name]);
%
% h = waitbar(0,'Loading dicoms, please wait...'); %Display waitbar
%
% %Import dicomsdicom_slider(SUBIM_x,100) %Using dicom_slider.m function for viewing
% % % % % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %Preallocating to save speed (With 2.08s, without, 2.56s on i5-6267U processor)
% temp = dicomread([num2str(path) f dicomnames(1).name]);
% Dicoms= int16(zeros(size(temp,1),size(temp,2), length(dicomnames)));
%
% for i = 1:length(dicomnames)
%     % Dicoms(:,:,i)= dicomread([num2str(path) f dicomnames(i).name]); %Using native stack
%     % Ali's figures are upside down
%     Dicoms(:,:,length(dicomnames)+1-i)= fliplr(dicomread([num2str(path) f dicomnames(i).name])); %Flipping the stack (up=down, left=right)
%     waitbar(i/length(dicomnames));
% end
% close(h);
%
% end

%---------------------------------------------------------------------

function dicom_slider(Dicoms,x) %(Stack, figure number)
% Function to use slider for image

switch nargin
    case 2
        fig=figure(x); %Uses the same figure
    case 1
        fig = figure;
end

Stack = Dicoms;

koko = size(Stack,3);

%fig=figure;
set(fig,'Name','Image','Toolbar','figure');%,...
%'NumberTitle','off')
% Create an axes to plot in
axes('Position',[.15 .05 .7 .9]);
% sliders for epsilon and lambda
slider1_handle=uicontrol(fig,'Style','slider','Max',koko,'Min',1,...
    'Value',2,'SliderStep',[1/(koko-1) 10/(koko-1)],...
    'Units','normalized','Position',[.02 .02 .14 .05]);
uicontrol(fig,'Style','text','Units','normalized','Position',[.02 .07 .14 .04],...
    'String','Choose frame');
% Set up callbacks
vars=struct('slider1_handle',slider1_handle,'Stack',Stack);
set(slider1_handle,'Callback',{@slider1_callback,vars});
plotterfcn(vars)
% End of main file

% Callback subfunctions to support UI actions
    function slider1_callback(~,~,vars)
        % Run slider1 which controls value of epsilon
        plotterfcn(vars)
    end

    function plotterfcn(vars)
        % Plots the image
        %imshow(vars.Stack(:,:,round(get(vars.slider1_handle,'Value'))));
        imagesc(vars.Stack(:,:,round(get(vars.slider1_handle,'Value'))));
        axis equal;
        title(num2str(get(vars.slider1_handle,'Value')));
        
    end
end

%---------------------------------------------------------------------

function [SUBIM_x, SUBIM_y] = orientation(Dicoms)
%Creates two image stacks from two different directions

%Preallocating for efficiancy
SUBIM_x = ones(size(Dicoms,3), size(Dicoms,1), size(Dicoms,2));
SUBIM_y = ones(size(Dicoms,3), size(Dicoms,2), size(Dicoms,1));

h = waitbar(0,'X-direction, please wait...'); %Display waitbar

for i = 1:size(Dicoms,2)
    for j = 1:size(Dicoms,3)
        SUBIM_x(j,:,i) = Dicoms(:,i,j);
    end
    waitbar(i/size(Dicoms,2));
end

close(h)


h = waitbar(0,'Y-direction, please wait...'); %Display waitbar

for i = 1:size(Dicoms,2)
    for j = 1:size(Dicoms,3)
        SUBIM_y(j,i,:) = Dicoms(:,i,j);
    end
    waitbar(i/size(Dicoms,2));
end

close(h)

end

%---------------------------------------------------------------------
function Coordinates = imdistancecalculator(slice_x, slice_y, xcoord, ycoord, aspectratio)
%Drawing and calculating

clear position

%Preallocating
thickness = 0;

t = figure('name', 'Click from top to bottom!', 'outerposition', [-317        1082        2560        1440]);

%t = figure('Name', 'Please press Enter when ready');

figchoice = 1000;
subplot(1,2,1)
imagesc(slice_x);
daspect(aspectratio)
% axis equal;

%title('X-direction ->', 'interpreter', 'none')
title('Choose your angle: X -> [1]');

line([ycoord,ycoord], [1, length(slice_x)],'Color','red','LineStyle','--') %Displays the y-slice

subplot(1,2,2)
imagesc(slice_y);
daspect(aspectratio)
% axis equal;

%title('Y-direction -^', 'interpreter', 'none')
title('Choose your angle: Y -^ [2]', 'Interpreter', 'none');

line([xcoord,xcoord], [1, length(slice_x)],'Color','red','LineStyle','--') %Displays the x-slice

%Allocating lines
%X-axis figure
leftpoint1 = ycoord;
rightpoint1 = ycoord;
%Y-axis figure
leftpoint2 = xcoord;
rightpoint2 = xcoord;

%Temporary
% % % % % question = menu('Satisfied?','1) Yes','2) No'); %Option to fall back should be added
question = 1;
if question == 1
    
    while figchoice ~= 13
        
        % This is just to assure that there exists no empty position
        if exist('position')
            if isempty(position.position1) || isempty(position.position2)
                clear position
            end
        end
        
        %Checking existence and displaying the current situation
        if exist('position')
            subplot(1,2,1)
            imagesc(slice_x);
            %title('X-direction ->', 'interpreter', 'none')
            title('Choose your angle: X -> [1]');
            xlabel('Y-axis')
            daspect(aspectratio)
            % axis equal;
            
            line([ycoord,ycoord], [1, length(slice_x)],'Color','red','LineStyle','--') %Displays the y-slice
            h1 = imline(gca,position.position1);
            subplot(1,2,2)
            imagesc(slice_y);
            title('Y-direction -^', 'interpreter', 'none')
            xlabel('X-axis')
            title('Choose your angle: Y -^ [2]', 'Interpreter', 'none');
            
            daspect(aspectratio)
            % axis equal;
            
            
            line([xcoord,xcoord], [1, length(slice_x)],'Color','red','LineStyle','--') %Displays the x-slice
            h2 = imline(gca,position.position2);
        end
        % drawsubplot(2)
        
        % % %
        pause %This pause waits for the user's input
        % % %
        %And returns what subplot is used
        figchoice = double(get(t,'CurrentCharacter'))
        
        
        
        % % % Playing with the lines
        if figchoice == 49 %1-button
            
            subplot(1,2,1)
            imagesc(slice_x);
            title('X-direction ->', 'interpreter', 'none')
            % title('Choose your angle: X -> [1]');
            daspect(aspectratio)
            % axis equal;
            
            line([ycoord,ycoord], [1, length(slice_x)],'Color','red','LineStyle','--') %Displays the y-slice
            if exist('position')
                h1 = imline(gca,position.position1);
            else
                h1 = imline;
            end
            position.position1 = wait(h1);
            subplot(1,2,2)
            imagesc(slice_y);
            title('Choose your angle: Y -^ [2]', 'Interpreter', 'none');
            daspect(aspectratio)
            % axis equal;
            
            
            line([xcoord,xcoord], [1, length(slice_x)],'Color','red','LineStyle','--') %Displays the x-slice
            
            %Calulating the new location
            % Z-axis
            highpoint = position.position1(1,2);
            lowpoint = position.position1(2,2);
            
            %These are used in the X-axis figure
            % With these the angle can change
            leftpoint1 = position.position1(1,1);
            rightpoint1 = position.position1(2,1);
            
            
            position.position2 = [leftpoint2 highpoint;  rightpoint2 lowpoint];
            
            %And displaying
            h2 = imline(gca,position.position2);
        elseif figchoice == 50
            
            subplot(1,2,2)
            imagesc(slice_y);
            title('Y-direction -^', 'interpreter', 'none')
            daspect(aspectratio)
            % axis equal;
            
            
            line([xcoord,xcoord], [1, length(slice_x)],'Color','red','LineStyle','--') %Displays the x-slice
            if exist('position')
                h2 = imline(gca,position.position2);
            else
                h2 = imline;
            end
            position.position2 = wait(h2);
            subplot(1,2,1)
            imagesc(slice_x);
            title('Choose your angle: X -> [1]');
            daspect(aspectratio)
            % axis equal;
            
            
            line([ycoord,ycoord], [1, length(slice_x)],'Color','red','LineStyle','--') %Displays the y-slice
            %Calculating the new location
            % Z-axis
            highpoint = position.position2(1,2);
            lowpoint = position.position2(2,2);
            %These are used in the X-axis figure
            % With these the angle can change
            leftpoint2 = position.position2(1,1);
            rightpoint2 = position.position2(2,1);
            
            position.position1 = [leftpoint1 highpoint;  rightpoint1 lowpoint];
            %And displaying
            h1 = imline(gca,position.position1);
        end
        
        
    end %while
    
    %     Dimensions = [(highpoint-lowpoint), (rightpoint1 - leftpoint1), (rightpoint2 - leftpoint2)]; %[Z X Y]
    Coordinates = [leftpoint2, leftpoint1, highpoint, rightpoint2, rightpoint1, lowpoint]; %[Z X Y]
    
elseif question == 2
    
    Coordinates = [1, 2, 3, 4, 5, 6];
    figchoice = 13; %Stops the execution
    
end %if question

end %function