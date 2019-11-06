function [name, tempGd, temp_I] = Partition_Analysis

%% m-file for analysing earlier created contrast agent partitions
%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, August / 2019

% clear all, close all, clc;

% You want extra figures plotted?
% 1 for yes
doyouwantfigures = 0;

% You can separately decide if you want all the profiles plotted
doyouwantprofilefigures = 0;


files = dir('*Rotated_RESULT_PROFILES*.mat');

for sample_i = 1:length(files)
%     for sample_i = [13,39,45] %FAILED
    %     sample_i = 1
    if doyouwantprofilefigures == 1 %If plotting, adding a pause
%         pause(2)
    end
    
    clearvars -except files sample_i name tempGd temp_I Eq_Gd Eq_I doyouwantfigures doyouwantprofilefigures
    
    filename = files(sample_i).name; %Reads the last/latest mat file
    load(filename);
    
    % 10 mgI/ml, 20 mgGd/ml
    % CORRECTED COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % By Ali % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % As a function of concentration
    % Calibration/Phantom_measurement.xlsx
    Gd50 = [36.32, 49.64];
    Gd90 = [37.34, 20.50];
    CA50 = [73.87, 89.58];
    CA90 = [29.53, 40.46];
    
    % Constants
    Water50 = 543.166; %Above values are normalized
    Water90 = -234.254;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fileorder = {'Baseline_50','Baseline_90','0h_50','0h_90', '30min_50','30min_90','1h_50','1h_90','2h_50','2h_90','6h_50','6h_90','10h_50','10h_90','23h_50','23h_90'};
    timepoints = [0, 0.5, 1, 2, 6, 10, 23]; %hours
    interpolation_n = 50; %You can determine how many points the data is interpolated into.
    depths = linspace(1,100, interpolation_n);
    
    % Inspecting first the 50 kV
    % h2 = waitbar(0,'Loading the files, please wait...'); %Display waitbar
    if doyouwantfigures == 1
        figure;
    end
    
    for location = 1:length(RESULT_PROFILES50)
        
        profiles50 = RESULT_PROFILES50{location};
        profiles90 = RESULT_PROFILES90{location};
        
        % Interpolation to 100 points (so that can be combined)
        P100_profile50(:,:,location) = interp1(linspace(1,100,size(profiles50,1)), profiles50, depths,'pchip');
        % % Checking
        % figure;
        % plot(linspace(1,100,length(profiles50)), profiles50, 'o-');
        % hold on;
        % plot(p100_profile50,'*');
        P100_profile90(:,:,location) = interp1(linspace(1,100,size(profiles90,1)), profiles90, depths,'pchip');
        
        
        if doyouwantfigures == 1
            % Depth-dependent attenuation profiles for all the measurement locations at each timepoints
            % 50kV
            subplot(1,2,1)
            % plot(profiles50);
            plot(P100_profile50(:,:,location));
            hold on;
            ylabel('Attenuation (AU)');
            xlabel('thickness (px)')
            % legend(fileorder{1:2:end}, 'location', 'se', 'interpreter', 'none')
            title('Raw 50 kV');
            ylim([-500 6000])
            
            % 90 kV
            % figure;
            subplot(1,2,2)
            % plot(profiles90);
            plot(P100_profile90(:,:,location));
            hold on;
            ylabel('Attenuation (AU)');
            xlabel('thickness (px)')
            % legend(fileorder{1:2:end}, 'location', 'se', 'interpreter', 'none')
            title('Raw 90 kV');
            ylim([-1000 2500])
            
        end
        % waitbar(location/6);
        
    end
    % close(h2) %Waitbar
    
    %
    %Playing with averages
    %
    
    %%%     Incorrect_partitions = 100*(Incorrect_partitions -Incorrect_partitions(1,:))./( Baselines(Baseline_i)-mean(Water)' );
    %%%     Averages = Incorrect_partitions./((Tube_volume-0.01*Incorrect_partitions.*Plug_volume)./Tube_volume);
    
    
    % % % % % Subtracting cartilage
    P100_profile50 = P100_profile50-(P100_profile50(:,1,:)); %Water doesn't need to be subtracted here as it is already done in calibration
    P100_profile90 = P100_profile90-(P100_profile90(:,1,:));
    
    
    P100_profile50(:,1,:) = []; %Removing zeroes
    P100_profile90(:,1,:) = [];
    
    % Normalized plots % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    if doyouwantfigures == 1
        % Plotting
        figure
        for location = 1:length(RESULT_PROFILES50)
            % 50kV
            subplot(1,2,1)
            % plot(profiles50);
            plot(P100_profile50(:,:,location));
            hold on;
            ylabel('Attenuation (AU)');
            xlabel('thickness (px)')
            % legend(fileorder{1:2:end}, 'location', 'se', 'interpreter', 'none')
            title('Normalized 50 kV');
            ylim([-500 6000])
            
            % 90 kV
            % figure;
            subplot(1,2,2)
            % plot(profiles90);
            plot(P100_profile90(:,:,location));
            hold on;
            ylabel('Attenuation (AU)');
            xlabel('thickness (px)')
            % legend(fileorder{1:2:end}, 'location', 'se', 'interpreter', 'none')
            title('Normalized 90 kV');
            ylim([-500 2500])
            
        end
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    
    %
    % Plotting
    
    meanprofiles50 = squeeze(mean(P100_profile50,1)); %Averaging the thickness
    stdprofiles50 = squeeze(std(P100_profile50,0,1));
    meanprofiles90 = squeeze(mean(P100_profile90,1));
    stdprofiles90 = squeeze(std(P100_profile90,0,1));
    
    
    % % If you want to look the depth-dependent profiles
    % meanprofiles50 = mean(P100_profile50,3)'; %Averaging the measured locations
    % stdprofiles50 = squeeze(std(P100_profile50,0,3)');
    % plot(depths,meanprofiles50, 'linewidth', 3) %Depth-dependent profiles at different timepoints
    
    
    
    % -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    
    if doyouwantprofilefigures == 1
        
        % % % % %FIGURES
        % Diffusion profiles for each point of interest
        Diffusion_50 = meanprofiles50';
        Diffusion_90 = meanprofiles90';
        
        hFig(sample_i) = figure;
        set(hFig(sample_i), 'Position', [256 356 1280 1280])
        subplot(2,1,1)
        plot(timepoints, Diffusion_50, 'o-', 'linewidth', 2)
        xlabel('Time (h)');
        ylabel('Attenuation (AU)')
        title('Locations, 50 kV Profile');
        lgd = legend(num2str(floor(depths)'), 'location' ,'southeast');
        title(lgd,'Location')
        
        subplot(2,1,2)
        plot(timepoints, Diffusion_90, 'o-', 'linewidth', 2)
        xlabel('Time (h)');
        ylabel('Attenuation (AU)')
        title('Locations, 90 kV Profile');
    end
    
    
    
    
    % % % % % % % % % % % % % % % % AVERAGES FOR THE LOCATIONS --------------------------------------------------
    % % % % % % % % % % % % % % % finalmeanprofiles50 = mean(meanprofiles50,2);
    % % % % % % % % % % % % % % % finalmeanprofiles90 = mean(meanprofiles90,2);
    % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % finalstdprofiles50 = std(meanprofiles50,0,1);
    % % % % % % % % % % % % % % % % finalstdprofiles90 = std(meanprofiles90,0,1);
    % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % The proper way is to calculate deviation between the sample means, like this:
    % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % finalstdprofiles50 = std( mean(P100_profile50,1),0,3); %First the sample, then the group
    % % % % % % % % % % % % % % % finalstdprofiles90 = std( mean(P100_profile90,1),0,3);
    % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % figure;
    % % % % % % % % % % % % % % % % % % plot(finalstdprofiles50);
    % % % % % % % % % % % % % % % % % % hold on;
    % % % % % % % % % % % % % % % % % % plot(finalstdprofiles50_2)
    % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % figure;
    % % % % % % % % % % % % % % % errorbar(timepoints,finalmeanprofiles50,finalstdprofiles50)
    % % % % % % % % % % % % % % % hold on;
    % % % % % % % % % % % % % % % errorbar(timepoints,finalmeanprofiles90,finalstdprofiles90)
    % % % % % % % % % % % % % % % legend(['50 kV'; '90 kV'], 'location', 'northwest')
    % % % % % % % % % % % % % % % title('Location''s average');
    
    %
    
    % Annin koodia
    C_CA = 10;% Change if you want partition: 10;
    C_Gd = 20;%20;
    
    
    for location = 1:length(RESULT_PROFILES50)
        if doyouwantfigures == 1 
            figure
        end
        
        for time = 1:7
            
            A = [CA50(1)*C_CA+CA50(2) Gd50(1)*C_Gd+Gd50(2); %Keeping the residue in results (zero points)
                CA90(1)*C_CA+CA90(2) Gd90(1)*C_Gd+CA90(2)];
            
            
            % b = [meanprofiles50(:,location) meanprofiles90(:,location)]';
            b = [P100_profile50(:,time,location) P100_profile90(:,time,location)]';
            
            
            concentration_profile = A\b;
            
            % yyaxis left
            
            if doyouwantfigures == 1
                % plot(timepoints, concentration_profile(1,:))
                plot(depths, concentration_profile(1,:),'b')
                % ylabel('Iodine (%)')
                
                hold on;
                % yyaxis right
                % plot(timepoints, concentration_profile(2,:), 'r')
                plot(depths, concentration_profile(2,:), 'r')
                title('Enzymatically + Mechanically - Partition in Cartilage','fontsize',14)
                legend('I', 'Gd');
                set(gca,'fontsize',14)
                % ylabel('Gadolinium (%)');
                % ylim([0 0.5]);
            end
            
            GADOLINIUM(:,time,location) = concentration_profile(2,:)';
            IODINE(:,time,location) = concentration_profile(1,:)';
            
        end
        
    end
    
    name{sample_i} = filename(1:end-29)
    tempGd{sample_i} = squeeze(mean(GADOLINIUM));
    temp_I{sample_i} = squeeze(mean(IODINE));
    
    
    Eq_Gd(:,sample_i) = tempGd{sample_i}(end,:)';
    Eq_I(:,sample_i)= temp_I{sample_i}(end,:)';
    
    
    
    
end %for sample_i
%%












