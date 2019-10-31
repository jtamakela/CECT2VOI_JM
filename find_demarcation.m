
function treshold_fit = find_orientation(thestack)
% Takes an image stack and finds cartilage-bone interface
% Outputs a fit of the demarcation


% ############################################################################################################################
% These pauses keep the plots from wondering to wrong windows. 
pause(0.1);

figure;
pause(0.1);

%Use this if you want all slices displayed
% for which_slice= 1:floor(size(orientationtesting,2))
% Otherwise (Mid slice):
for which_slice= floor(size(thestack,2)/2) %1:floor(size(orientationtesting,2))
    
    midslice = squeeze(thestack(:,which_slice,:));
    %
    %         imagesc(midslice);
    %         view(90, 90)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fororientation = midslice>400; %THIS THRESHOLD HAS NOT BEEN OPTIMIZED
    % However, doesn't have to point exactly to bone. The trends are still there
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
    pause(0.1);

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
    
    pause(0.1);
    plot(function_for_fit, points, 'go');
    pause(0.2);
    
    
    %     figure;
    if ~isempty(points)
        treshold_fit = fit(points',function_for_fit','poly1');
        
        hold on;
        plot([0:100].*treshold_fit.p1+treshold_fit.p2,[0:100], 'g--', 'linewidth' ,2)
        %     view(90, 90) %If you want to make a separate figure
        %     axis ij
        
    end
    
    
end


% ############################################################################################################################
end
