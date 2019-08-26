%% m-file for analysing earlier created contrast agent partitions
%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, August / 2019

clear all, close all, clc;

filename = dir('*RESULT_PROFILES*.mat');

filename = filename(end).name; %Reads the last/latest mat file
load(filename);

fileorder = {'Baseline_50','Baseline_90','0h_50','0h_90', '30min_50','30min_90','1h_50','1h_90','2h_50','2h_90','6h_50','6h_90','10h_50','10h_90','23h_50','23h_90'};
timepoints = [0, 0.5, 1, 2, 6, 10, 23]; %hours
interpolation_n = 20; %You can determine how many points the data is interpolated into. 
depths = linspace(1,100, interpolation_n);

% Inspecting first the 50 kV
h2 = waitbar(0,'Loading the files, please wait...'); %Display waitbar

for location = 1:6

profiles50 = RESULT_PROFILES50{location};
profiles90 = RESULT_PROFILES90{location};

% Interpolation to 100 points (so that can be combined)
P100_profile50(:,:,location) = interp1(linspace(1,100,length(profiles50)), profiles50, depths,'pchip');
% % Checking
% figure;
% plot(linspace(1,100,length(profiles50)), profiles50, 'o-');
% hold on;
% plot(p100_profile50,'*');
P100_profile90(:,:,location) = interp1(linspace(1,100,length(profiles90)), profiles90, depths,'pchip');


% % 50kV
% figure;
% subplot(1,2,1)
% % plot(profiles50);
% plot(P100_profile50(:,:,location));
% ylabel('Attenuation (AU)');
% xlabel('thickness (px)')
% % legend(fileorder{1:2:end}, 'location', 'se', 'interpreter', 'none')
% title('50 kV');
% ylim([-500 5000])
% 
% % 90 kV
% % figure;
% subplot(1,2,2)
% % plot(profiles90);
% plot(P100_profile90(:,:,location));
% ylabel('Attenuation (AU)');
% xlabel('thickness (px)')
% % legend(fileorder{1:2:end}, 'location', 'se', 'interpreter', 'none')
% title('90 kV');
% ylim([-1000 2000])


% waitbar(location/6);

end
close(h2) %Waitbar


%Playing with averages


%%%     Incorrect_partitions = 100*(Incorrect_partitions -Incorrect_partitions(1,:))./( Baselines(Baseline_i)-mean(Water)' ); 
%%%     Averages = Incorrect_partitions./((Tube_volume-0.01*Incorrect_partitions.*Plug_volume)./Tube_volume);

% Subtracting cartilage
P100_profile50 = P100_profile50-P100_profile50(:,1,:);
P100_profile90 = P100_profile90-P100_profile90(:,1,:);

P100_profile50(:,1,:) = []; %Removing zeroes
P100_profile90(:,1,:) = [];


% Plotting

meanprofiles50 = squeeze(mean(P100_profile50,1));
stdprofiles50 = squeeze(std(P100_profile50,0,1));

meanprofiles90 = squeeze(mean(P100_profile90,1));
stdprofiles90 = squeeze(std(P100_profile90,0,1));

%FIGURES
figure;
plot(timepoints,meanprofiles50, 'linewidth', 3)
xlabel('Time (h)');
ylabel('Attenuation (AU)')
% hold on;
% plot(timepoints,meanprofiles50+stdprofiles50, 'k--')
% plot(timepoints,meanprofiles50-stdprofiles50, 'r--')
title('50 kV Profile');

figure;
plot(timepoints,meanprofiles90, 'linewidth', 3)
xlabel('Time (h)');
ylabel('Attenuation (AU)')
% hold on;
% plot(timepoints,meanprofiles90+stdprofiles90, 'k--')
% plot(timepoints,meanprofiles90-stdprofiles90, 'r--')
title('90 kV Profile');
%% --

% % Diffusion rate
% 
% Diffusion_50 = meanprofiles50';
% Diffusion_90 = meanprofiles90';
% 
% figure;
% subplot(2,1,1)
% plot(timepoints, Diffusion_50, 'o-', 'linewidth', 2)
% xlabel('Time (h)');
% ylabel('Attenuation (AU)')
% title('50 kV')
% legend(num2str(floor(depths)'), 'location' ,'northeastoutside')
% 
% subplot(2,1,2)
% plot(timepoints, Diffusion_90, 'o-', 'linewidth', 2)
% xlabel('Time (h)');
% ylabel('Attenuation (AU)')
% title('90 kV')
% legend([''], 'location', 'southeastoutside');
% legend('boxoff')

finalmeanprofiles50 = mean(meanprofiles50,1);
finalmeanprofiles90 = mean(meanprofiles90,1);

% finalstdprofiles50 = std(meanprofiles50,0,1);
% finalstdprofiles90 = std(meanprofiles90,0,1);

% The proper way is to calculate deviation between the sample means, like this: 

finalstdprofiles50 = std( mean(P100_profile50,1),0,3); %First the sample, then the group
finalstdprofiles90 = std( mean(P100_profile90,1),0,3);

% % % figure; 
% % % plot(finalstdprofiles50);
% % % hold on;
% % % plot(finalstdprofiles50_2)

figure;
errorbar(timepoints,finalmeanprofiles50,finalstdprofiles50)
















