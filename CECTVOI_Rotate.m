function  NEW_DATA = CECTVOI_Rotate % Doesn't need to be a function

%% m-file for further rotating the earlier created CECT image VOIs
% Uses CECTVOI_Alignment and find_demarcation to turn

%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, October / 2019


% clear all, close all, clc;

% Are we saving?
saving = 1; % 1 for saving


%Testing if this returns the angles

% Don't change the file name, won't work with CECTVOI_Alignment
filelist = dir('*_VOI_data*.mat');


for whichfile = 1:length(filelist)  %Which sample
    clearvars -except saving filelist whichfile
    close all
    
    
    % Going through all the timepoints
    filename = filelist(whichfile).name %Reads the last mat file
    
    load(filename); %Imports DATA
    
    % DATA{1,2} has the VOIs stored
    CECT_original = DATA{1,2}; %2 energies, 6 locations, 8 timepoints
    
    
    for whichlocation = 1:size(CECT_original,1) % %Which of the measured locations
        
        %Check the angle
        angles = CECTVOI_Alignment(whichfile, whichlocation);
        
        
        % Going through all timepoints & locations
        h2 = waitbar(0,'Rotating the files, please wait...'); %Display waitbar
        
        for whichtimepoint = 1:16
            waitbar(whichtimepoint/length(DATA{1,2}));
            
            % Rotate X
            new_CECT90 = imrotate3(CECT_original{whichlocation, whichtimepoint}, angles(1), [1 0 0]);
            % Turning 90 deg
            new_rot_CECT90 = rot90(new_CECT90);
            % Rotate Y
            new_new_rot_CECT90 = imrotate3(new_rot_CECT90, angles(2), [1 0 0]);
            % Turning back to original
            NEW_CECT{whichlocation, whichtimepoint} = rot90(new_new_rot_CECT90,-1);
            
        end
        close(h2)
        
    end
    
    savename = filename(1:end-14);
    
    NEWDATA = {DATA{1,1}, NEW_CECT, DATA{1,3}};
    
    % Don't overwrite
    if saving == 1
        save([savename, '_RotatedVOI_data', num2str(length(dir(['*', num2str(savename),'_Rotated*']))+1), '.mat'], 'NEWDATA');
    end
    
    
    
end


