clear all, close all,% clc; 

%cd /media/janne/Elements/uCT registrations/2Ri3

% This just for checking the registration

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

% -----------------------------------------------------------------------------------------------------------

niinames = dir('*.nii');
h2 = waitbar(0,'Loading the files, please wait...'); %Display waitbar

counter = 1;

% V = single(zeros(512,512,512));

tic 
for i = order
   
filuname = niinames(i).name;

V = niftiread(filuname);
info = niftiinfo(filuname);

figure;
imagesc(V(:,:,300), [-1500 5000]);
title(filuname, 'interpreter', 'none')
    waitbar(counter/length(niinames));

counter = counter +1 ;
end
close(h2)
toc;


