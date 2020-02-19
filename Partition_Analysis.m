function [names_CECT, tempGd, temp_I] = Partition_Analysis

%% m-file for analysing earlier created contrast agent partitions
%% Developed for triple contrast images
%% The code is available at https://github.com/jtamakela/
%% (c) Janne T.A. Mäkelä, August / 2019

clear all, close all, clc;

% You want extra figures plotted?
% 1 for yes
doyouwantfigures = 0;

% You can separately decide if you want all the profiles plotted
doyouwantprofilefigures = 0;

% If you want to analyze only specific samples
% chosen = {'3LR1',	'4RR1',	'5RR1',	'6Li3',	'7Li3',	'8Ri3',	'9Ri3',	'10Ri3'};
% chosen = {'5LR2'};
plotcolors = ['yrcgbkyr'];

files = dir('*Rotated_RESULT_PROFILES*.mat');

for sample_i = 1:length(files)
    %     for sample_i = [13,39,45] %FAILED
    %     sample_i = 1
    
    clearvars -except files sample_i names_CECT tempGd temp_I Eq_Gd Eq_I doyouwantfigures doyouwantprofilefigures chosen plotcolors
    
    filename = files(sample_i).name; %Reads the last/latest mat file

    
    
    % Adding a rule that the execution is continued only if the file matches the chosen
    % Comment 'exist' out if you want this to execute
    if ~exist('chosen') || contains(num2str(cell2mat(chosen)),filename(1:4))
        %         if isempty( strfind(num2str(cell2mat(chosen)),filename(1:4)) ) % Same thing        
        
        
        load(filename);
        
        % CORRECTED COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % By Ali % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % As a function of concentration
        % Calibration/Phantom_measurement.xlsx
        % HOUNSFIELDS
%         From Ali
        Gd50 = [37.34, 20.50];
        Gd90 = [36.32, 49.64];
        CA50 = [73.87, 89.58];
        CA90 = [29.53, 40.46];
        % ARBITRARY UNITS (RAW DATA)

        
        % 2 months later. Gd perhaps gone bad?
%         Gd50 = [45.7 19.87];
%         Gd90 = [45.1 0.098];
%         CA50 = [73.98 16.22];
%         CA90 = [28.36 10.63];

        
        % Constants, not used
%         Water50 = 543.166; %Above values are normalized
%         Water90 = -234.254;


% FROM ABHISEK (HU)
% αI(50 kV) = 54.75 CI + 65.41
% αI(90 kV) = 35.21 CI + 59.91
% b) gadolinium
% αGd(50 kV) = 34.28 CGd + 40.31
% αGd(90 kV) = 54.81 CGd + 53.77

% FROM ABHISEK (AU)
% αI(50 kV) = 119 CI + 142 
% αI(90 kV) =  46 CI + 78
% b) gadolinium
% αGd(50 kV) = 75 CGd + 76
% αGd(90 kV) = 72 CGd + 100 




        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fileorder = {'Baseline_50','Baseline_90','0h_50','0h_90', '30min_50','30min_90','1h_50','1h_90','2h_50','2h_90','6h_50','6h_90','10h_50','10h_90','23h_50','23h_90'};
        timepoints = [0, 0.5, 1, 2, 6, 10, 23]; %hours
        interpolation_n = 50; %You can determine how many points the data is interpolated into.
        depths = linspace(1,100, interpolation_n);
        
        % Inspecting first the 50 kV
        % h2 = waitbar(0,'Loading the files, please wait...'); %Display waitbar

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
   
            % waitbar(location/6);
            
        end
        % close(h2) %Waitbar
        
        %
        %Playing with averages
        %
        
        %%%     Incorrect_partitions = 100*(Incorrect_partitions -Incorrect_partitions(1,:))./( Baselines(Baseline_i)-mean(Water)' );
        %%%     Averages = Incorrect_partitions./((Tube_volume-0.01*Incorrect_partitions.*Plug_volume)./Tube_volume);
        
        
        % % % % % Subtracting cartilage
        P100_profile50_normal = P100_profile50-(P100_profile50(:,1,:)); %Water doesn't need to be subtracted here as it is already done in calibration
        P100_profile90_normal = P100_profile90-(P100_profile90(:,1,:));
        
        
        P100_profile50_normal(:,1,:) = []; %Removing baseline
        P100_profile90_normal(:,1,:) = [];
               
        
        % Normalized plots % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        if doyouwantfigures == 1
            % Plotting
            figure
            for time = 1:size(P100_profile50,2)
                
                                % Depth-dependent attenuation profiles for all the measurement locations at each timepoints
                % 50kV
                subplot(2,2,1)
                % plot(profiles50);
                plot(squeeze(P100_profile50(:,time,:)), 'color', num2str(plotcolors(time)));
                hold on;
                ylabel('Attenuation (AU)');
                xlabel('thickness (px)')
                % legend(fileorder{1:2:end}, 'location', 'se', 'interpreter', 'none')
                title('Raw 50 kV');
                ylim([-500 6000])
                
                % 90 kV
                % figure;
                subplot(2,2,2)
                % plot(profiles90);
                plot(squeeze(P100_profile90(:,time,:)), 'color', num2str(plotcolors(time)));
                hold on;
                ylabel('Attenuation (AU)');
                xlabel('thickness (px)')
                % legend(fileorder{1:2:end}, 'location', 'se', 'interpreter', 'none')
                title('Raw 90 kV');
                ylim([-1000 2500])
                
            end
            for time = 1:size(P100_profile50_normal,2)
                
                % 50kV
                subplot(2,2,3)
                % plot(profiles50);
                plot(squeeze(P100_profile50_normal(:,time,:)), 'color', num2str(plotcolors(time+1)));
                hold on;
                ylabel('Attenuation (AU)');
                xlabel('thickness (px)')
                % legend(fileorder{1:2:end}, 'time', 'se', 'interpreter', 'none')
                title('Normalized 50 kV');
                ylim([-500 6000])
                
                % 90 kV
                % figure;
                subplot(2,2,4)
                % plot(profiles90);
                plot(squeeze(P100_profile90_normal(:,time,:)), 'color', num2str(plotcolors(time+1)));
                hold on;
                ylabel('Attenuation (AU)');
                xlabel('thickness (px)')
                % legend(fileorder{1:2:end}, 'time', 'se', 'interpreter', 'none')
                title('Normalized 90 kV');
                ylim([-500 2500])
                
            end
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        
        %
        % Plotting
        
        meanprofiles50 = squeeze(mean(P100_profile50_normal,1)); %Averaging the thickness
        stdprofiles50 = squeeze(std(P100_profile50_normal,0,1));
        meanprofiles90 = squeeze(mean(P100_profile90_normal,1));
        stdprofiles90 = squeeze(std(P100_profile90_normal,0,1));
        
        
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
            lgd = legend(num2str([1:6]'), 'location' ,'southeast');
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
        
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % 10 mgI/ml, 20 mgGd/ml

        C_CA = 10;% Change if you want partition: 10;
        C_Gd = 20;%20;
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        
        for location = 1:length(RESULT_PROFILES50)
            if doyouwantfigures == 1
                figure
            end
            
            
            for time = 1:size(P100_profile50_normal,2)
                
% % % % %                 
                A = [CA50(1)*C_CA Gd50(1)*C_Gd ; %Keeping the residue in results (zero points)
                    CA90(1)*C_CA Gd90(1)*C_Gd];
                
                % A*x = b
                
                % b = [meanprofiles50(:,location) meanprofiles90(:,location)]';
                b = [P100_profile50_normal(:,time,location) P100_profile90_normal(:,time,location)]';
%                 b = [845 570]';%A coctail, 5 I, 10 Gd
                
                
                concentration_profile = A\b;
                
                % yyaxis left
                
                if doyouwantfigures == 1
                    % plot(timepoints, concentration_profile(1,:))
                    yyaxis left
                    plot(depths, concentration_profile(1,:), '-', 'color', num2str(plotcolors(time+1)))
                    % ylabel('Iodine (%)')
                    
                    hold on;
                    % yyaxis right
                    % plot(timepoints, concentration_profile(2,:), 'r')
                    yyaxis right
                    plot(depths, concentration_profile(2,:), '--', 'color', num2str(plotcolors(time+1)))
                    title('Partition in Cartilage','fontsize',14)
                    legend(num2str(timepoints'), 'location', 'southeast');
                    set(gca,'fontsize',14)
                    % ylabel('Gadolinium (%)');
                    % ylim([0 0.5]);
                end
                
                %Taking into account only the chosen depth
                GADOLINIUM(:,time,location) = concentration_profile(2,31:40)';
                IODINE(:,time,location) = concentration_profile(1,31:40)';
                
            end
            
        end
        
        names_CECT{sample_i} = filename(1:end-29);
        tempGd{sample_i} = squeeze(mean(GADOLINIUM));
        temp_I{sample_i} = squeeze(mean(IODINE));
        
        Gd_profile_ave(:,time,location,sample_i) = concentration_profile(2,:)';
        I_profile_ave(:,time,location,sample_i) = concentration_profile(1,:)';
        
        Eq_Gd(:,sample_i) = tempGd{sample_i}(end,:)';
        Eq_I(:,sample_i)= temp_I{sample_i}(end,:)';
        
%         break % Adding this rule for debugging. Inspects only the first 
    end %if 'chosen' exists and filename contains the chosen
    
    
end %for sample_i

% save('CECT_data_60_80.mat', 'names_CECT', 'tempGd', 'temp_I');







