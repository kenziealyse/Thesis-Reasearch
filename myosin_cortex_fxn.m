function [max_MhcA_vec, T] = myosin_cortex_fxn (frames, bleb_frame, bleb_MhcA_data, plotFrames)


frameRate = 2.14; %Time interval between frames


bleb_MhcA = bleb_MhcA_data(:,bleb_frame:bleb_frame+frames); % To isolote bleb or boundary of interest

% Remove nan cells from data

first_nan_idx = find(isnan(bleb_MhcA(:,1)));

if first_nan_idx >1 
     bleb_MhcA_gd = bleb_MhcA(1:first_nan_idx(1)-1,:);
else
     bleb_MhcA_gd = bleb_MhcA;
end


% SCALE INTENSITY DATA

[n,m] = size(bleb_MhcA_gd);
bleb_MhcA = zeros(n,m);
bleb_MhcA(:,1)= bleb_MhcA_gd(:,1);
max_MhcA_vec = zeros(1,frames); % I think this needs to be frames and not m (or m-1)

for i = 1:frames
    cyto_MhcA_mean = mean(bleb_MhcA_gd(1:10,i+1));
    bleb_MhcA(:,i+1)= bleb_MhcA_gd(:,i+1)/cyto_MhcA_mean;
    max_MhcA_vec(i)= max(bleb_MhcA(:,i+1));   % For generating plots later on   
end


% PLOT SCALED INTENSITY DATA

max_MhcA = max(max_MhcA_vec);

figure

if plotFrames == 1

for i = 1:m-1
    subplot(3,5,i)
    plot(bleb_MhcA(:,1),bleb_MhcA(:,i+1),'g-','LineWidth',3)
    grid on
    ylim([0 max(max_MhcA)])
    yticks(0:0.5:max(max_MhcA))
    xticks(0:0.5:max(bleb_MhcA(:,1)))
    xlim([0 max(bleb_MhcA(:,1))])
    xlabel('Distance (in microns)')
    ylabel('Fluoresence (AU)')
    title([num2str((i-1)*frameRate) ' Secs'])
    legend('MhcA')
    hold off
end

end

for i = 1:m-1
    T(i) = (i-1)*frameRate;
end

end
