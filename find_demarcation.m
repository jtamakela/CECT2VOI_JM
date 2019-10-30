
function treshold_fit = find_orientation(thestack)
% Takes an image stack and finds cartilage-bone interface
% Outputs a fit of the demarcation


    % Let's check if the figure needs to be automatically aligned
    % 90 keV data goes great for trying to find bone cartilage interface
    orientationtesting = thestack;
    
    % ############################################################################################################################
    %The following needs to be in its own function
    
    figure;
    
    
    for which_slice= floor(size(orientationtesting,2)/2) %1:floor(size(orientationtesting,2))
        midslice = squeeze(orientationtesting(:,which_slice,:));
        %
        %         imagesc(midslice);
        %         view(90, 90)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fororientation = midslice>400; %THIS THRESHOLD HAS NOT BEEN OPTIMIZED
        % However, doesn't have to point to bone. The trends are still there
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
