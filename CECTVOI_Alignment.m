function angles = CECTVOI_Alignment(whichfile, whichlocation)
%% Returns the orientation of a given location based on the baseline
%% This needs to be run first so that the orientation can be found. 


%% m-file for further rotating the earlier created CECT image VOIs
%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, October / 2019

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
% 6,2 ei toimi!!!!!!!!!!!!!!!!!!!!!!!!

% Added to the function call
% % % % % whichfile = 7; %Reads the numbered last mat file
% % % % % whichlocation = 1; %There should be a six of these for the ponies


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
    
   
    
    % This is what we rotate
    [ROTATED_RESULT_PROFILES90, angles] = Rotate_me_a_matrix(CECT90{whichlocation,timepoint});
%     EI KÄÄNNY




end

function [backto_CECT90, angles] = Rotate_me_a_matrix(thematrix);
% Finding the orientation -----------------------------------------------------------------------------------------------
fit_coefficients = find_demarcation(thematrix)



%The angle
rotangle_x = atan( fit_coefficients.p1 ) * (180/pi)

% Rotating a slice -----------------------------------------------------------------------------------------------------

%
% Rotating The matrix from the other angle
% -> Turning the angle 90 deg and redoing the previous -----------------------------------------------------------------------------------------------------

new_CECT90 = imrotate3(thematrix, rotangle_x, [1 0 0]);

% Testing a slice
% sizeIn = size(new_CECT90);
% hFigOriginal = figure;
% hAxOriginal  = axes;
% slice(double(new_CECT90),sizeIn(2)/2,sizeIn(1)/2,sizeIn(3)/2);
% grid on, shading interp, %colormap gray



showmeasneakpeak(new_CECT90);



% Turning 90 deg
new_rot_CECT90 = rot90(new_CECT90);
% Finding the angle
fit_coefficients = find_demarcation(new_rot_CECT90)
rotangle_y = atan( fit_coefficients.p1 ) * (180/pi)

% Rotating
new_new_rot_CECT90 = imrotate3(new_rot_CECT90, rotangle_y, [1 0 0]);
showmeasneakpeak(new_new_rot_CECT90);

% Turning back to original
backto_CECT90 = rot90(new_new_rot_CECT90,-1);
% THIS IS SUPPOSED TO BE THE SAME AS new_CECT90 BECAUSE MID SLICE IS NOT ROTATED

% % Display
sizeIn = size(backto_CECT90);
hFigOriginal = figure;
hAxOriginal  = axes;
slice(double(fliplr(backto_CECT90)),sizeIn(2)/2,sizeIn(1)/2,sizeIn(3)/2);
grid on, shading interp, %colormap gray


angles = [rotangle_x, rotangle_y];


% -----------------------------------------------------------------------------------------------------------------------
end

function showmeasneakpeak(matrixtoshow)

thehalf = floor(size(matrixtoshow,2)/2);

sneakpeak = squeeze(matrixtoshow(:,thehalf,:)); %midslice

figure;
imagesc(sneakpeak)
view(90, 90)


end


