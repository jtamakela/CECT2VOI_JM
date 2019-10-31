function angles = CECTVOI_Rotate(whichfile, whichlocation)

%% m-file for further rotating the earlier created CECT image VOIs
% Uses CECTVOI_Alignment and find_demarcation to turn

%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, October / 2019


clear all, close all, clc; 

%Testing if this returns the angles

whichfile = 1;  %Which sample
whichlocation = 1; %Which of the measured locations

test = CECTVOI_Alignment(whichfile, whichlocation)






end


