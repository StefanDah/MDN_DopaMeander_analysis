%% Load data

load('\postAnalysis.mat')

filepath='placeholder';

smoothing_factor = 2;

sampling_rate = 20000;
data_reduction_for_plotting = 400;
binsize_factor = 0.2; %200 ms
numPokes = 10;
%t_ramp = 0.5; % 0.5s ramp time
t_poke = 2; % 2s poke hold
interstim_interval = 10; %10s between pokes
% Length of analysis window
windowSizeinS = 0.5;
% Longer analysis window for behavior analysis
windowSizeinSBehavior = 1;

binsize = 4000;
bin = 20;

for k = 1 : length(analysis)
    analysis(k).xveloc_in_mm = analysis(k).VelocX(:,1)*8.79;
    analysis(k).zveloc_in_degree_per_s = analysis(k).VelocZ(:,1)*158.9;
    analysis(k).VM_medfilt = medfilt1(analysis(k).VM,1000);
end

%% trajectory

for k = 1 : length(analysis)
    analysis(k).motion = [analysis(k).VelocX, analysis(k).VelocY, analysis(k).VelocZ];
end
figure('Name', 'Trajectory')
%tiledlayout(round(length(analysis)/4), 4)
%tiledlayout(1,2)

for k = 1 %: length(analysis)

    currentMotion  = analysis(k).motion(1:data_reduction_for_plotting:end, :);

    allTrajectories = ball2trajectory(currentMotion);


    nexttile %%Generate new tile for each trajectory
    title(analysis(k).ID,'Interpreter','none') %% Interpreter 'none' required for preseting underscores (_) properly

    % figure('Name','Trajecory');
    hold on
    xline(0,'Color',[0.8,0.8,0.8]);
    yline(0,'Color',[0.8,0.8,0.8]);
    caxis('manual')
    caxis([0 6])
    T = allTrajectories;
    x = T(1:end,2);
    y = T(1:end,3);
    z = zeros(size(x));

    scatter(x(1), y(1), [3000], '.')

    scatter(x(1:end-2), y(1:end-2), [], 'c',  'filled', 'o');

    colormap(flipud(cool))
    colorbar;
    caxis([0 1])


    plot(0,0,...
        'Marker','o', 'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    %Total distance walked
    d = hypot(diff(x), diff(y));                            % Distance Of Each Segment
    d_tot = sum(d);

    xlabel('X dimension [mm]');
    ylabel('Y dimension [mm]');
    axis equal;

end

mean_velocity = d_tot/((length(currentMotion)*400)/sampling_rate); %%mean velocity in mm/s

%% spike binning
clearvars binstart binstarts binsize

binstarts = nan;
binsize = binsize_factor*sampling_rate;

for k = 1 : length(analysis)

    binstarts = 1:binsize:length(analysis(k).spikes);
    spikesbinned = nan;
    for bin = 1 : length(binstarts)-1
        spikesbinned(:,bin) = sum(analysis(k).spikes(binstarts(bin):binstarts(bin+1)),'omitnan');
    end
    analysis(k).spikesbinned = spikesbinned/binsize_factor;

end

%% x-velocity binning

clearvars binstart binstarts binsize

binstarts = nan;
binsize = binsize_factor*sampling_rate;

for k = 1 : length(analysis)

    binstarts = 1:binsize:length(analysis(k).VelocX);
    xvelocbinned=nan;
    for bin = 1 : length(binstarts)-1
        xvelocbinned(:,bin) = sum(analysis(k).VelocX(binstarts(bin):binstarts(bin+1)),'omitnan');
    end
    xvelocbinned = xvelocbinned/binsize;
    analysis(k).xvelocbinned = xvelocbinned;

end

%% z-velocity binning

clearvars binstart binstarts binsize

binstarts = nan;
binsize = binsize_factor*sampling_rate;

for k = 1 : length(analysis)

    binstarts = 1:binsize:length(analysis(k).VelocZ);
    zvelocbinned=nan;
    for bin = 1 : length(binstarts)-1
        zvelocbinned(:,bin) = sum(analysis(k).VelocZ(binstarts(bin):binstarts(bin+1)),'omitnan');
    end
    zvelocbinned = zvelocbinned/binsize;
    analysis(k).zvelocbinned = zvelocbinned;

end
%% spike binning during stimulation window

clearvars binstart binstarts binsize spikesbinned binned_spikerate_right

binstarts = nan;
binsize = binsize_factor*sampling_rate;
binned_spikerate= cell(1,length(analysis));

for k = 1 : length(analysis)
    temp_spikes = [];


    for m = 1 : numPokes
        window_on = analysis(k).trigger_on(m)-1*sampling_rate; %1s before poking
        window_off = analysis(k).trigger_off(m)+1*sampling_rate; %1s afte poking

        binstarts = window_on:binsize:window_off;
        for bin = 1 : length(binstarts)-1
            temp_spikes{m,bin} = sum(analysis(k).spikes(binstarts(bin):binstarts(bin+1)),'omitnan');
        end
    end
    analysis(k).binnedSpikerate = cellfun(@(x)x/binsize_factor, temp_spikes);
    analysis(k).binnedSpikerateMean = mean(analysis(k).binnedSpikerate);
end

%% x-velocity in mm binning

clearvars binstart binstarts binsize

binstarts=nan;
binsize = binsize_factor*sampling_rate;

for k = 1 : length(analysis)

    binstarts = 1:binsize:length(analysis(k).xveloc_in_mm);
    xvelocbinned=nan;
    for bin = 1 : length(binstarts)-1
        xvelocbinned(:,bin) = sum(analysis(k).xveloc_in_mm(binstarts(bin):binstarts(bin+1)),'omitnan');
    end
    xvelocbinned = xvelocbinned/binsize;
    analysis(k).xveloc_in_mm_binned = xvelocbinned;

end
%% X velocity binning during stimulation window

clearvars binstart binstarts binsize binned_xveloc_right

binstarts = nan;
binsize = binsize_factor*sampling_rate;

%binned_xveloc_right = cell(1,length(analysis));

for k = 1 : length(analysis)

    temp_xveloc = [];

    for m = 1 : numPokes
        window_on = analysis(k).trigger_on(m)-1*sampling_rate; %1s before poking
        window_off = analysis(k).trigger_off(m)+1*sampling_rate; %1s after poking

        binstarts = window_on:binsize:window_off;
        for bin = 1 : length(binstarts)-1
            temp_xveloc{m,bin} = sum(analysis(k).xveloc_in_mm...
                (binstarts(bin):binstarts(bin+1)),'omitnan');
        end
    end
    analysis(k).binnedXveloc = cellfun(@(x)x/binsize, temp_xveloc);
    analysis(k).binnedXvelocMean = mean(analysis(k).binnedXveloc);

end

%% Z velocity binning during stimulation window
for k = 1 : length(analysis)
    analysis(k).zveloc_in_mm = analysis(k).VelocZ(:,1)*8.79;
end

%%%%% Right Antenna
clearvars binstart binstarts binsize

binstarts = nan;
binsize = binsize_factor*sampling_rate;

%binned_zveloc_right = cell(1,length(analysis));

for k = 1 : length(analysis)

    temp_zveloc = [];
    temp_zveloc_degree = [];

    for m = 1 : numPokes

        window_on = analysis(k).trigger_on(m)-1*sampling_rate; %1s before poking
        window_off = analysis(k).trigger_off(m)+1*sampling_rate; %1s after poking

        binstarts = window_on:binsize:window_off;
        for bin = 1 : length(binstarts)-1
            temp_zveloc{m,bin} = sum(analysis(k).zveloc_in_mm(binstarts(bin):binstarts(bin+1)),'omitnan');
            temp_zveloc_degree{m,bin} = sum(analysis(k).zveloc_in_degree_per_s(binstarts(bin):binstarts(bin+1)),'omitnan');

        end

        % if to_be_removed_ID(k,m) == 1
        %     temp_zveloc(k,m) = NaN;
        % else
        % end



    end
    analysis(k).binnedZveloc = cellfun(@(x)x/binsize, temp_zveloc);
    analysis(k).binnedZvelocMean = mean(analysis(k).binnedZveloc);
    analysis(k).binnedZvelocDegree = cellfun(@(x)x/binsize, temp_zveloc_degree);
    analysis(k).binnedZvelocDegreeMean = mean(analysis(k).binnedZvelocDegree);
end

%% VM binning during stimulation window

%%%%% Right Antenna
clearvars binstart binstarts binsize

binstarts = nan;
binsize = binsize_factor*sampling_rate;

%binned_zveloc_right = cell(1,length(analysis));

for k = 1 : length(analysis)

    temp_VM = [];
    temp_medVM = [];
    for m = 1 : numPokes

        window_on = analysis(k).trigger_on(m)-1*sampling_rate; %1s before poking
        window_off = analysis(k).trigger_off(m)+1*sampling_rate; %1s after poking

        binstarts = window_on:binsize:window_off;
        for bin = 1 : length(binstarts)-1
            temp_VM{m,bin} = sum(analysis(k).VM(binstarts(bin):binstarts(bin+1)),'omitnan');
            temp_medVM{m,bin} = sum(analysis(k). VM_medfilt(binstarts(bin):binstarts(bin+1)),'omitnan');
        end

        % if to_be_removed_ID(k,m) == 1
        %     temp_zveloc(k,m) = NaN;
        % else
        % end

    end
    analysis(k).binnedVM = cellfun(@(x)x/binsize, temp_VM);
    analysis(k).binnedVMMean = mean(analysis(k).binnedVM);
    analysis(k).binnedVM_medfilt = cellfun(@(x)x/binsize, temp_medVM);
    analysis(k).binnedVMMean_medfilt = mean(analysis(k).binnedVM_medfilt);

end
%% Binning for cross correlation

clearvars binstart binstarts binsize temp_spikes temp_xveloc temp_zveloc_degree

%Use a smaller binsize to increase the resolution of the cross correlation
binsize_factor_xcorr = 0.05; %50 ms
binstarts = nan;
binsize_xcorr = binsize_factor_xcorr*sampling_rate;

for k = 1 : length(analysis)
    temp_spikes = [];
    temp_xveloc = [];
    temp_zveloc_degree = [];

    binstarts = 1:binsize_xcorr:length(analysis(k).spikes);

    for bin = 1 : length(binstarts)-1
        temp_spikes(:,bin) = sum(analysis(k).spikes...
            (binstarts(bin):binstarts(bin+1)),'omitnan');
        temp_xveloc(:,bin) = sum(analysis(k).xveloc_in_mm...
            (binstarts(bin):binstarts(bin+1)),'omitnan');
        temp_zveloc_degree(:,bin) = sum(analysis(k).zveloc_in_degree_per_s...
            (binstarts(bin):binstarts(bin+1)),'omitnan');

    end
    temp_zveloc_degree = temp_zveloc_degree/binsize_xcorr;
    analysis(k).zvelocbinned_degree_xcorr = temp_zveloc_degree;
    temp_xveloc = temp_xveloc/binsize_xcorr;
    analysis(k).xveloc_in_mm_binned_xcorr = temp_xveloc;
    analysis(k).spikesbinned_xcorr = temp_spikes/binsize_factor;
end

%% Analyze poker effect

clearvars right_antenna left_antenna right_antenna_currentInjection...
    left_antenna_currentInjection

for i = 1 : length(analysis)
    if isnan(analysis(i).currentInjection) == true

        if strcmp(analysis(i).antenna, 'right')
            right_antenna(i)= analysis(i);

        else
            left_antenna(i) = analysis(i);

        end
    else
        if strcmp(analysis(i).antenna, 'right')
            right_antenna_currentInjection(i)= analysis(i);

        else
            left_antenna_currentInjection(i) = analysis(i);

        end
    end
end

right_antenna = right_antenna(~cellfun(@isempty,{right_antenna.VM}));
left_antenna = left_antenna(~cellfun(@isempty,{left_antenna.VM}));

% right_antenna_currentInjection = right_antenna_currentInjection(~cellfun...
%     (@isempty,{right_antenna_currentInjection.VM}));
% left_antenna_currentInjection = left_antenna_currentInjection(~cellfun...
%     (@isempty,{left_antenna_currentInjection.VM}));

%save('postPoker_right_currentInjection','right_antenna_currentInjection','-v7.3');
%save('postPoker_left_currentInjection','left_antenna_currentInjection','-v7.3');

%save('postAnalysis', 'analysis', '-v7.3');



%% Quality control


right_antenna = right_antenna(~cellfun(@(x)x>-25,{right_antenna.medianVM}));
left_antenna = left_antenna(~cellfun(@(x)x>-25,{left_antenna.medianVM}));
right_antenna = right_antenna(~cellfun(@(x)x>25,{right_antenna.meanSpikeHz}));
left_antenna = left_antenna(~cellfun(@(x)x>25,{left_antenna.meanSpikeHz}));

% save('postPoker_right','right_antenna','-v7.3');
% save('postPoker_left','left_antenna','-v7.3');

% right_antenna = right_antenna(~cellfun(@(x)x<5,{right_antenna.meanSpikeHz}));
% left_antenna = left_antenna(~cellfun(@(x)x<5,{left_antenna.meanSpikeHz}));
%
% right_antenna = right_antenna(~cellfun(@(x)x>20,{right_antenna.meanSpikeHz}));
% left_antenna = left_antenna(~cellfun(@(x)x>20,{left_antenna.meanSpikeHz}));

numCells = length(right_antenna);


%% Area under the curve (AUC) testing to filter out trials with walking before poking
clearvars t v timeSpan velocitySpan startIndex endIndex data to_be_removed_ID

% Preallocate array for IDs of traces to be excluded
to_be_removed_ID = {zeros(length(right_antenna), numPokes), zeros(length(right_antenna), numPokes)};
% Assign datasets
datasets = {right_antenna, left_antenna};
datasets_names = {'right_antenna', 'left_antenna'};

for n = 1 : length(datasets)

    clearvars t v timeSpan velocitySpan startIndex endIndex data
    data = datasets{n};

    figure
    for i = 1 : length(data)

        % velcostiy data
        v = (data(i).xveloc_in_mm); 
        % Timeseries of velocity data
        t = (1:length(v))/sampling_rate; 

        nexttile
        plot(t, v, 'LineWidth', 1);
        for k = 1 : numPokes

            % Start 1 second (sampling_rate/2) before poking
            startIndexPre = data(i).trigger_on(k)-(sampling_rate/2); 
            % End at start of poking
            endIndexPre = data(i).trigger_on(k); 
            % Start at poking
            startIndexPoke = data(i).trigger_on(k);  
            % End at the end of poke (2s long)
            endIndexPoke = data(i).trigger_on(k)+(sampling_rate); 
            % Start at the end of poking
            startIndexPost = data(i).trigger_off(k); 
            % End 1 s after end of poke
            endIndexPost = data(i).trigger_off(k)+(sampling_rate/2); 


            % Extract the relevant portion of the time and velocity vectors
            timeSpanPre = t(startIndexPre:endIndexPre);
            velocitySpanPre = v(startIndexPre:endIndexPre);
            timeSpanPoke = t(startIndexPoke:endIndexPoke);
            velocitySpanPoke= v(startIndexPoke:endIndexPoke);
            timeSpanPost = t(startIndexPost:endIndexPost);
            velocitySpanPost= v(startIndexPost:endIndexPost);

            % Calculate the area under the curve using the trapz function
            areaUnderCurvePre = trapz(timeSpanPre, velocitySpanPre);
            areaUnderCurvePoke = trapz(timeSpanPoke, velocitySpanPoke);
            areaUnderCurvePost = trapz(timeSpanPost, velocitySpanPost);

            % Plot the time series
            hold on;
            %Highlight the specified time span
            highlight_span_pre = plot([startIndexPre/sampling_rate, endIndexPre/sampling_rate],...
                [0, 0], 'm', 'LineWidth', 5);
            highlight_span_poke = plot([startIndexPoke/sampling_rate, endIndexPoke/sampling_rate],...
                [0, 0], 'r', 'LineWidth', 5);
            highlight_span_post = plot([startIndexPost/sampling_rate, endIndexPost/sampling_rate],...
                [0, 0], 'g', 'LineWidth', 5);

            % Fill the region corresponding to the AUC
            areaPre = area(timeSpanPre, velocitySpanPre, 'FaceColor', 'm', 'EdgeColor', 'm');
            areaPoke = area(timeSpanPoke, velocitySpanPoke, 'FaceColor', 'r', 'EdgeColor', 'r');
            areaPost = area(timeSpanPost, velocitySpanPost, 'FaceColor', 'g', 'EdgeColor', 'g');

            %If AUC 0.5 s before the poke (Pre) exceeds threshold, mark the AUC
            %with an asterisk and save the ID of the poke (fly,poke)
            if areaUnderCurvePre > 0.1
                text(startIndexPre/sampling_rate, max(areaPre.YData)+2, '*','Color','red','FontSize',14)
                to_be_removed_ID{n}(i,k) = 1;
            else
            end

        end

        % Add labels and title
        xlabel('Time (s)');
        ylabel('Forward Velocity in mm/s');
        title('X Velocity with AUC');
        subtitle(['No.', num2str(i) , datasets_names{n}])

        % Add legend
        % legend([highlight_span_pre], ['Pre Poke']);
        legend('','Pre', 'Poke', 'Post');

    end

end

% How many trials are exlucded
nTotalRight = numel(to_be_removed_ID{1});
nExcludedRight = sum(to_be_removed_ID{1},'all');
nTotalLeft = numel(to_be_removed_ID{2});
nExcludedLeft = sum(to_be_removed_ID{2},'all');
%% Plotting mean and binned spikerates during stimulation with 1 s before and 1 s after
clearvars binned_spikerate_right_mean binned_spikerate_left_mean binned_xveloc_right_mean...
    binned_xveloc_left_mean binned_zveloc_right_mean binned_zveloc_left_mean xbin

%Save means of all cells in one cell
meanSpikerateCellRight = arrayfun(@(s) s.binnedSpikerateMean, right_antenna, 'UniformOutput', false);
%Calculate mean from all cells
meanSpikerateOverallRight = mean(cat(1, meanSpikerateCellRight{:}));

meanSpikerateCellLeft = arrayfun(@(s) s.binnedSpikerateMean, left_antenna, 'UniformOutput', false);
meanSpikerateOverallLeft = mean(cat(1, meanSpikerateCellLeft{:}));

meanVMCellRight = arrayfun(@(s) s.binnedVMMean, right_antenna, 'UniformOutput', false);
meanVMOverallRight = mean(cat(1, meanVMCellRight{:}));

meanVMCellLeft = arrayfun(@(s) s.binnedVMMean, left_antenna, 'UniformOutput', false);
meanVMOverallLeft = mean(cat(1, meanVMCellLeft{:}));

meanVM_medCellRight = arrayfun(@(s) s.binnedVMMean_medfilt, right_antenna, 'UniformOutput', false);
meanVM_medOverallRight = mean(cat(1, meanVM_medCellRight{:}));

meanVM_medCellLeft = arrayfun(@(s) s.binnedVMMean_medfilt, left_antenna, 'UniformOutput', false);
meanVM_medOverallLeft = mean(cat(1, meanVM_medCellLeft{:}));

meanXvelocCellRight = arrayfun(@(s) s.binnedXvelocMean, right_antenna, 'UniformOutput', false);
meanXvelocOverallRight = mean(cat(1, meanXvelocCellRight{:}));

meanXvelocCellLeft = arrayfun(@(s) s.binnedXvelocMean, left_antenna, 'UniformOutput', false);
meanXvelocOverallLeft = mean(cat(1, meanXvelocCellLeft{:}));

meanZvelocCellRight = arrayfun(@(s) s.binnedZvelocDegreeMean, right_antenna, 'UniformOutput', false);
meanZvelocOverallRight = mean(cat(1, meanZvelocCellRight{:}));

meanZvelocCellLeft = arrayfun(@(s) s.binnedZvelocDegreeMean, left_antenna, 'UniformOutput', false);
meanZvelocOverallLeft = mean(cat(1, meanZvelocCellLeft{:}));


%xbin = [1:bin]/(sampling_rate/binsize);
xbin = ((binsize/2)/sampling_rate):binsize/sampling_rate:(bin*(binsize/sampling_rate));

% Create a colormap
cmap = hsv(numel(right_antenna));

%Plot right antenna firing rate and x velocity
figure('name', 'Mean per Cell -- 2s poke 0.2 bins')

tiledlayout(3,2)
nexttile;
for i = 1 : length(right_antenna)
    color = cmap(i, :);
    plot(xbin, cell2mat(meanSpikerateCellRight(i)), 'Color', color)
    hold on
end
plot(xbin, meanSpikerateOverallRight,'  k', 'Linewidth', 3)
%rectangle('Position',[2 0 0.5 30], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
title('Spike rate')
ylim([0 45])
%rectangle('Position',[2 0 0.5 30], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 0 2 45], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')

nexttile;
for i = 1 : length(left_antenna)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanSpikerateCellLeft(i)), 'Color', color)
    hold on
end
plot(xbin, meanSpikerateOverallLeft, 'k','Linewidth', 2)
%rectangle('Position',[2 0 0.5 30], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 0 2 45], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('Spike rate')
ylim([0 45])
ylabel('Spike Frequency (Hz)');


nexttile;
for i = 1 : length(right_antenna)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanXvelocCellRight(i)), 'Color', color)
    hold on
end
plot(xbin, meanXvelocOverallRight, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('X velocity')
ylim([-2 1])
ylabel('X Velocity (mm)');

nexttile;
for i = 1 : length(left_antenna)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanXvelocCellLeft(i)), 'Color', color)
    hold on
end
plot(xbin, meanXvelocOverallLeft, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('X velocity')
ylim([-2 1])
xlabel('X Velocity');

nexttile;
for i = 1 : length(right_antenna)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanZvelocCellRight(i)), 'Color', color)
    hold on
end
plot(xbin, meanZvelocOverallRight, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('Z velocity')
% ylim([-2.5 2.5])
ylabel('Z Velocity (mm)');
xlabel('Time (s)')

nexttile;
for i = 1 : length(right_antenna)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanZvelocCellLeft(i)), 'Color', color)
    hold on
end
plot(xbin, meanZvelocOverallLeft, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('Z velocity')
% ylim([-2.5 2.5])
xlabel('Time (s)')


%Calculate standard deviation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdev_right = std(cat(1, meanSpikerateCellRight{:}));
stdev_left = std(cat(1, meanSpikerateCellLeft{:}));

y_upper_limit_X_Right = meanSpikerateOverallRight + stdev_right;
y_lower_limit_X_Right = meanSpikerateOverallRight - stdev_right;
y_upper_limit_X_Left = meanSpikerateOverallLeft + stdev_left;
y_lower_limit_X_Left = meanSpikerateOverallLeft - stdev_left;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('name', 'Overall mean -- 2s poke 0.2 bins')

t = tiledlayout('flow');
title(t, "Overview" + newline + "N =" + num2str(length(right_antenna)))

nexttile

%plot Left antenna x velocity with SD as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_X_Right fliplr(y_lower_limit_X_Right)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanSpikerateOverallRight, 'Linewidth', 3)
fill([xbin fliplr(xbin)], [y_upper_limit_X_Left fliplr(y_lower_limit_X_Left)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanSpikerateOverallLeft, 'Linewidth', 2)

title('Spike rate')
legend( '', 'left antenna', '', 'right antenna','Location','northwest')
rectangle('Position',[1 0 2 45], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
ylim([0 (max(y_upper_limit_X_Right)+5)])
ylabel('Spike Frequency (Hz)');


%Calculate standard deviation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdev_right = std(cat(1, meanVM_medCellRight{:}));
stdev_left = std(cat(1, meanVM_medCellLeft{:}));

y_upper_limit_X_Right = meanVM_medOverallRight + stdev_right;
y_lower_limit_X_Right = meanVM_medOverallRight - stdev_right;
y_upper_limit_X_Left = meanVM_medOverallLeft + stdev_left;
y_lower_limit_X_Left = meanVM_medOverallLeft - stdev_left;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nexttile

fill([xbin fliplr(xbin)], [y_upper_limit_X_Right fliplr(y_lower_limit_X_Right)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanVM_medOverallRight, 'Linewidth', 3)
fill([xbin fliplr(xbin)], [y_upper_limit_X_Left fliplr(y_lower_limit_X_Left)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanVM_medOverallLeft, 'Linewidth', 3)

title('median VM velocity')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
ylim([min(y_lower_limit_X_Right)-0.25 max(y_upper_limit_X_Right)+0.25])
ylabel('median VM velocity (VM))');


%Calculate standard deviation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdev_right = std(cat(1, meanXvelocCellRight{:}));
stdev_left = std(cat(1, meanXvelocCellLeft{:}));

y_upper_limit_X_Right = meanXvelocOverallRight + stdev_right;
y_lower_limit_X_Right = meanXvelocOverallRight - stdev_right;
y_upper_limit_X_Left = meanXvelocOverallLeft + stdev_left;
y_lower_limit_X_Left = meanXvelocOverallLeft - stdev_left;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile

fill([xbin fliplr(xbin)], [y_upper_limit_X_Right fliplr(y_lower_limit_X_Right)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallRight, 'Linewidth', 3)
fill([xbin fliplr(xbin)], [y_upper_limit_X_Left fliplr(y_lower_limit_X_Left)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanXvelocOverallLeft, 'Linewidth', 3)

title('X velocity')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
ylim([min(y_lower_limit_X_Right)-0.25 max(y_upper_limit_X_Right)+0.25])
ylabel('X Velocity (mm))');

%Calculate standard deviation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdev_right = std(cat(1, meanZvelocCellRight{:}));
stdev_left = std(cat(1, meanZvelocCellLeft{:}));

y_upper_limit_X_Right = meanZvelocOverallRight + stdev_right;
y_lower_limit_X_Right = meanZvelocOverallRight - stdev_right;
y_upper_limit_X_Left = meanZvelocOverallLeft + stdev_left;
y_lower_limit_X_Left = meanZvelocOverallLeft - stdev_left;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile

fill([xbin fliplr(xbin)], [y_upper_limit_X_Right fliplr(y_lower_limit_X_Right)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallRight, 'Linewidth', 3)
fill([xbin fliplr(xbin)], [y_upper_limit_X_Left fliplr(y_lower_limit_X_Left)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanZvelocOverallLeft, 'Linewidth', 3)

title('Z velocity')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
ylim([min(y_lower_limit_X_Right)-0.25 max(y_upper_limit_X_Left)+0.25])
ylabel('Z Velocity (mm))');
xlabel('Time (s)');

print(['Overview' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
print(['Overview' '.png'], '-dpng','-r300', '-vector')


%Plot all cells with firing rate, x and z velocities next to each other
figure('name', 'Cell by cell comparison -- 2s poke 0.2 bins')

t = tiledlayout(3,length(right_antenna),"TileSpacing","compact");

for i = 1 : length(right_antenna)
    nexttile;
    plot(xbin, meanSpikerateCellRight{i})
    hold on
    plot(xbin, meanSpikerateCellLeft{i})
    set(gca, 'XTickLabel', [])
    ylim([0 35])
end

for i = 1 : length(right_antenna)
    nexttile;
    plot(xbin, meanXvelocCellRight{i})
    hold on
    plot(xbin, meanXvelocCellLeft{i})
    set(gca, 'XTickLabel', [])
    ylim([-1.75 3.5])
end

for i = 1 : length(right_antenna)
    nexttile;
    plot(xbin, meanZvelocCellRight{i})
    hold on
    plot(xbin, meanZvelocCellLeft{i})
    % ylim([-4.25 3.25])
end

title(t,'All cells')
xlabel(t,'Time (s)')


% print(['individual_cells' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['individual_cells' '.png'], '-dpng','-r300', '-vector')

%% Calculate mean firing rate and create boxplots
clearvars spikesPrePokeRight spikesDuringPokeRight spikesPostPokeRight spikesPrePokeLeft...
    spikesDuringPokeLeft spikesPostPokeLeft

%Size of poking analysis window in s
analysis_window = windowSizeinS * sampling_rate;


for i = 1 : length(right_antenna)
    for k = 1 : numPokes
        %
        % if to_be_removed_ID{1}(i,k) == 1
        %     spikesPrePokeRight{i}(k) = NaN;
        %     spikesDuringPokeRight{i}(k) = NaN;
        %     spikesPostPokeRight{i}(k) = NaN;
        %     continue
        % else
        windowOnRight = right_antenna(i).trigger_on(k);
        windowOffRight = right_antenna(i).trigger_off(k);
        %Sum spikes pre poke
        spikesPrePokeRight{i}(k) = sum(right_antenna(i).spikes...
            (windowOnRight-analysis_window:windowOnRight),'omitnan')*2;
        %Sum spikes during poke
        spikesDuringPokeRight{i}(k) = sum(right_antenna(i).spikes...
            (windowOnRight:windowOnRight+analysis_window),'omitnan')*2;
        %Sum spikes post poke
        spikesPostPokeRight{i}(k) = sum(right_antenna(i).spikes...
            (windowOffRight:windowOffRight+analysis_window),'omitnan')*2;
    end
end
% end


for i = 1 : length(left_antenna)
    for k = 1 : numPokes

        % if to_be_removed_ID{2}(i,k) == 1
        %     spikesPrePokeLeft{i}(k) = NaN;
        %     spikesDuringPokeLeft{i}(k) = NaN;
        %     spikesPostPokeLeft{i}(k) = NaN;
        %     continue
        % else
        windowOnLeft = left_antenna(i).trigger_on(k);
        windowOffLeft = left_antenna(i).trigger_off(k);
        %Sum spikes pre poke
        spikesPrePokeLeft{i}(k) = sum(left_antenna(i).spikes...
            (windowOnLeft-analysis_window:windowOnLeft),'omitnan')*2;
        %Sum spikes during poke
        spikesDuringPokeLeft{i}(k) = sum(left_antenna(i).spikes...
            (windowOnLeft:windowOnLeft+analysis_window),'omitnan')*2;
        %Sum spikes post poke
        spikesPostPokeLeft{i}(k) = sum(left_antenna(i).spikes...
            (windowOffLeft:windowOffLeft+analysis_window),'omitnan')*2;
    end
end
% end

%Calculate means of spikes
% meanPreRight = cellfun(@mean, spikesPrePokeRight);
meanPreRight = cellfun(@(x) mean(x, 'omitnan'), spikesPrePokeRight)
%meanPreLeft = cellfun(@mean, spikesPrePokeLeft);
meanPreLeft = cellfun(@(x) mean(x, 'omitnan'), spikesPrePokeLeft)
% meanPokeRight = cellfun(@mean, spikesDuringPokeRight);
meanPokeRight = cellfun(@(x) mean(x, 'omitnan'), spikesDuringPokeRight)
% meanPokeLeft = cellfun(@mean, spikesDuringPokeLeft);
meanPokeLeft = cellfun(@(x) mean(x, 'omitnan'), spikesDuringPokeLeft)
% meanPostRight = cellfun(@mean, spikesPostPokeRight);
meanPostRight = cellfun(@(x) mean(x, 'omitnan'), spikesPostPokeRight)
% meanPostLeft = cellfun(@mean, spikesPostPokeLeft);
meanPostLeft = cellfun(@(x) mean(x, 'omitnan'), spikesPostPokeLeft)

%Statistics
[pWilPreVsPokeSpikesRight] = signrank(meanPreRight, meanPokeRight);
[pWilPostVsPokeSpikesRight] = signrank(meanPostRight, meanPokeRight);

[pWilPreVsPokeSpikesLeft] = signrank(meanPreLeft, meanPokeLeft);
[pWilPostVsPokeSpikesLeft] = signrank(meanPostLeft, meanPokeLeft);

%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean firing rates pre, during 1st and 2nd s and post stimulus -- RIGHT')
t = tiledlayout('flow');
title(t, "Spike Rate Pre, During and Post Poke")
%Labels for groups
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna)))
groupNames = ["Pre"; "Poke " + num2str(windowSizeinS) + "sec"; "Post"];
groups = [1; 2; 3;];
plot_y = [meanPreRight; meanPokeRight; meanPostRight];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels([groupNames]);

hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
text(1,mean(meanPreRight), ['p=' num2str(round(pWilPreVsPokeSpikesRight,4))])
text(3,mean(meanPostRight), ['p=' num2str(round(pWilPostVsPokeSpikesRight,4))])
ylabel('Spike Rate (Hz)')
ylim([0 25])
%print median
medianFiringRateRightPre = h.med(1).YData(1)
medianFiringRateRightPoke = h.med(2).YData(1)
medianFiringRateRightPost = h.med(3).YData(1)

% subtitle("Excluded trials: " + num2str(nExcludedRight)+"/"+ num2str(nTotalRight));

nexttile
title("LEFT" + newline + "N =" + num2str(length(right_antenna)))
groupNames = ["Pre"; "Poke " + num2str(windowSizeinS) + "sec"; "Post"];
groups = [1; 2; 3;];
plot_y = [meanPreLeft; meanPokeLeft; meanPostLeft];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
text(1,mean(meanPreLeft), ['p=' num2str(round(pWilPreVsPokeSpikesLeft,4))])
text(3,mean(meanPostLeft), ['p=' num2str(round(pWilPostVsPokeSpikesLeft,4))])
ylabel('Spike Rate (Hz)')
ylim([0 25])
%print median
medianFiringRateLeftPre = h.med(1).YData(1)
medianFiringRateLeftPoke = h.med(2).YData(1)
medianFiringRateLeftPost = h.med(3).YData(1)

% subtitle("Excluded trials: " + num2str(nExcludedLeft)+"/"+ num2str(nTotalLeft));

%Pool spikes from left and right antenna
meanPrePooled = [meanPreLeft meanPreRight];
meanPokePooled = [meanPokeLeft meanPokeRight];
meanPostPooled = [meanPostLeft meanPostRight];

[pWilPreVsPokePooled] = signrank(meanPrePooled, meanPokePooled);
[pWilPostVsPokePooled] = signrank(meanPostPooled, meanPokePooled);

% nexttile
% title("Pooled" + newline + "N =" + num2str(length(meanPrePooled)))
% groupNames = ["Pre"; "Poke 1st sec"; "Post"];
% groups = [1; 2; 3;];
% plot_y = [meanPrePooled; meanPokePooled; meanPostPooled];
% h = boxplot2(plot_y',groups,'whisker', 0);
% set(h.out, 'Marker', 'none')
% xticks([1,2,3]);
% xticklabels([groupNames]);
% hold on
% %Use scatter to plot filled circles of each mean
% s = scatter(groups, plot_y);
% set(s, 'Marker', {'none'});
% line(groups, plot_y)
% text(1,mean(meanPreLeft), ['p=' num2str(round(pWilPreVsPokePooled,4))])
% text(3,mean(meanPostLeft), ['p=' num2str(round(pWilPostVsPokePooled,4))])
% ylabel('Spike Rate (Hz)')
% ylim([0 25])
set(gcf,'position',[400, 100, 400,1000])

%print median
medianFiringRatePooledPre = h.med(1).YData(1);
medianFiringRatePooledPoke = h.med(2).YData(1);
medianFiringRatePooledPost = h.med(3).YData(1);
% 
 % print(['spikerate_boxplot_500ms' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['spikerate_boxplot_500ms' '.png'], '-dpng','-r300', '-vector')


%% pool pre and post from left and right and baseline subtract for Spike Rate

preLeft = [];
preRight = [];
prePooled = [];
postPooled = [];
PokeRightSpikeBaselineSub = [];
PokeLeftSpikeBaselineSub = [];
PrePokeBaselineSub = [];
PostPokeBaselineSub = [];
PrePostPooled = [];
PostRightSpikeBaselineSub = [];
PostLeftSpikeBaselineSub = [];

for i = 1 : numCells

% Pool baseline from pre and post poke     
preLeft(i) = mean(spikesPrePokeLeft{i},'omitnan');
preRight(i) = mean(spikesPrePokeRight{i},'omitnan');
prePooled(i) = mean([spikesPrePokeRight{i} spikesPrePokeLeft{i}],'omitnan');
postPooled(i) = mean([spikesPostPokeRight{i} spikesPostPokeLeft{i}],'omitnan');

% Pool baseline from pre and post poke     
PrePostPooled(i) = mean([spikesPrePokeRight{i} spikesPrePokeLeft{i}...
     spikesPostPokeRight{i} spikesPostPokeLeft{i}],'omitnan');

end

% Subtract left baseline from left and right from right
PokeRightSpikeBaselineSub = meanPokeRight -  preRight;
PokeLeftSpikeBaselineSub = meanPokeLeft -  preLeft;
PostRightSpikeBaselineSub = meanPostRight - preRight;
PostLeftSpikeBaselineSub = meanPostLeft - preLeft;

[pWilPreVsPokeLeftBaselineSub] = signrank(preLeft-preLeft, PokeLeftSpikeBaselineSub);
[pWilPreVsPokeRightBaselineSub] = signrank(preRight-preRight, PokeRightSpikeBaselineSub);
[pWilPostVsPokeLeftBaselineSub] = signrank(PostLeftSpikeBaselineSub, PokeLeftSpikeBaselineSub);
[pWilPostVsPokeRightBaselineSub] = signrank(PostRightSpikeBaselineSub, PokeRightSpikeBaselineSub);

figure
t = tiledlayout('flow');
title(t, "Spike Rate Pre, During and Post Poke" + newline + "Baseline Subtracted")
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna)))
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [preRight-preRight; PokeRightSpikeBaselineSub; PostRightSpikeBaselineSub];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels(groupNames);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
ylabel('Spike Rate (Hz)')
ylim([-15 15])
text(1,13, ['Pre vs Right' newline 'p=' num2str(round(pWilPreVsPokeRightBaselineSub,4))])
text(3,13, ['Post vs Right' newline 'p=' num2str(round(pWilPostVsPokeRightBaselineSub,4))])
set(gcf,'position',[400, 100, 400,1000])

nexttile
title("LEFT" + newline + "N =" + num2str(length(right_antenna)))
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [preLeft-preLeft; PokeLeftSpikeBaselineSub; PostLeftSpikeBaselineSub];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels(groupNames);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
ylabel('Spike Rate (Hz)')
ylim([-15 15])
text(1,13, ['Pre vs Left' newline 'p=' num2str(round(pWilPreVsPokeLeftBaselineSub,4))])
text(3,13, ['Post vs Left' newline 'p=' num2str(round(pWilPostVsPokeLeftBaselineSub,4))])


PokeRightSpikeBaselineSub = [];
PokeLeftSpikeBaselineSub = [];
PrePokeBaselineSub = [];
PostPokeBaselineSub = [];

% Subtract baseline from pre, poke and post
PokeRightSpikeBaselineSub = meanPokeRight- PrePostPooled;
PokeLeftSpikeBaselineSub = meanPokeLeft - PrePostPooled;
PrePokeBaselineSub = prePooled - PrePostPooled;
PostPokeBaselineSub = postPooled - PrePostPooled;

figure
title("Spike Rate" + newline + "Baseline Pre+Post subtracted")
groupNames = ["Pre"; "Left"; "Right"; "Post"];
groups = [1; 2; 3; 4;];
plot_y = [PrePokeBaselineSub; PokeLeftSpikeBaselineSub; PokeRightSpikeBaselineSub; PostPokeBaselineSub];
h = boxplot2(plot_y',groups);
xticks([1,2,3,4]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
line([1 2], [12.5 12.5])
line([1 2], [14 14])
line([3 4], [12.5 12.5])
line([3 4], [14 14])

text(1,13, ['Pre vs Left' newline 'p=' num2str(round(pWilPreVsPokeLeftBaselineSub,4))])
text(3,13, ['Post vs Left' newline 'p=' num2str(round(pWilPostVsPokeLeftBaselineSub,4))])
text(1,14.5, ['Pre vs Right' newline 'p=' num2str(round(pWilPreVsPokeRightBaselineSub,4))])
text(3,14.5, ['Post vs Right' newline 'p=' num2str(round(pWilPostVsPokeRightBaselineSub,4))])
ylabel('Spike Rate (Hz)')
ylim([-5 20])

set(gcf,'position',[400, 100, 400,1000])

%% Sort left and right response

% Based on the median change in response to antennal touch sort left and
% right antenna per fly in groups. Group 1 contains the stronger response
% while Group 2 has the lower one

% Compare median change in spike rate

groupOne = zeros(2,numCells);
groupTwo = zeros(2,numCells);

for i = 1 : numCells
    % Check if the response in left is higher than add it to group 1 and
    % right value to group 2, if not vice versa 
    if meanPokeLeft(i) > meanPokeRight(i)
        groupOne(1,i) = meanPokeLeft(i);
        % A one indicated left antenna, 0 is by default right antenna
        groupOne(2,i) = 1;
        groupTwo(1,i) = meanPokeRight(i);

    else
        groupTwo(1,i) = meanPokeLeft(i);
         % A one indicated left antenna, 0 is by default right antenna
        groupTwo(2,i) = 1;
        groupOne(1,i) = meanPokeRight(i);
    end
end
% figure
% plot_y = [groupOne(1,:); groupTwo(1,:)];
% groups = [1;2];
% groupNames = ["1","2"];
% boxplot2(plot_y',groups)
% xticks([1,2]);
% xticklabels(groupNames);

hold on
colors = 'mk';

h = gscatter(ones(1,9), groupOne(1,:), groupOne(2,:), colors, [],25,'filled');
h = gscatter(ones(1,9)+1, groupTwo(1,:), groupTwo(2,:), colors,[], 25,'filled');
legend("Right Antenna", "Left Antenna",'Location','northeast')
hold off
xticks([1,2]);
xlim([0.5 2.5])
ylabel("Spike Rate during Poking")
xlabel("Grouped Responses")

contrast = [];

% Conrast analysis using c = (x-y)/(x+y)
for i = 1 : numCells
contrast(i) = (groupOne(1,i) - groupTwo(1,i)) / (groupOne(1,i) + groupTwo(1,i))
end

figure
plot_y = contrast';
h = boxplot2(plot_y, 1,'whisker', 0);
set(h.out, 'Marker', 'none')
hold on
scatter(ones(1,9), contrast, "filled", 'jitter','on', 'jitterAmount',0.05)
xticks(1);
ylabel("Contrast - Group 1 vs 2")
ylim([0 1])

 % print(['Contrast_500ms' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['Contrast_500ms' '.png'], '-dpng','-r300', '-vector')

%% Calculate mean VM  and create boxplots
clearvars medVMPrePokeRight medVMDuringPokeRight medVMPostPokeRight

%Size of poking analysis window in s
analysis_window = windowSizeinS * sampling_rate;

for i = 1 : length(right_antenna)
    for k = 1 : numPokes

        % if to_be_removed_ID{1}(i,k) == 1
        %     medVMPrePokeRight{i}(k) = NaN;
        %     medVMDuringPokeRight{i}(k) = NaN;
        %     medVMPostPokeRight{i}(k) = NaN;
        %     continue
        % else

            windowOnRight = right_antenna(i).trigger_on(k);
            windowOffRight = right_antenna(i).trigger_off(k);
            %Sum spikes 1 s pre poke
            medVMPrePokeRight{i}(k) = mean(right_antenna(i).VM_medfilt...
                (windowOnRight-analysis_window:windowOnRight),'omitnan');
            %Sum spikes 1st second of poke
            medVMDuringPokeRight{i}(k) = mean(right_antenna(i).VM_medfilt...
                (windowOnRight:windowOnRight+analysis_window),'omitnan');
            %Sum spikes after poke
            medVMPostPokeRight{i}(k)  = mean(right_antenna(i).VM_medfilt...
                (windowOffRight:windowOffRight+analysis_window),'omitnan');

        end
    end
% end

for i = 1 : length(left_antenna)
    for k = 1 : numPokes

        % if to_be_removed_ID{2}(i,k) == 1
        %     medVMPrePokeLeft{i}(k) = NaN;
        %     medVMDuringPokeLeft{i}(k) = NaN;
        %     medVMPostPokeLeft{i}(k) = NaN;
        %     continue
        % else

            windowOnLeft = left_antenna(i).trigger_on(k);
            windowOffLeft = left_antenna(i).trigger_off(k);
            %Sum spikes 1 s pre poke
            medVMPrePokeLeft{i}(k) = mean(left_antenna(i).VM_medfilt...
                (windowOnLeft-analysis_window:windowOnLeft),'omitnan');
            %Sum spikes 1st second of poke
            medVMDuringPokeLeft{i}(k) = mean(left_antenna(i).VM_medfilt...
                (windowOnLeft:windowOnLeft+analysis_window),'omitnan');
            %Sum spikes after poke
            medVMPostPokeLeft{i}(k)  = mean(left_antenna(i).VM_medfilt...
                (windowOffLeft:windowOffLeft+analysis_window),'omitnan');
        end
    end
% end

%Calculate means of X velocity
% meanPreVMRight = cellfun(@mean, medVMPrePokeRight);
meanPreVMRight = cellfun(@(x) mean(x, 'omitnan'), medVMPrePokeRight)
% meanPokeVMRight = cellfun(@mean, medVMDuringPokeRight);
meanPokeVMRight = cellfun(@(x) mean(x, 'omitnan'), medVMDuringPokeRight)
% meanPostVMRight = cellfun(@mean, medVMPostPokeRight);
meanPostVMRight = cellfun(@(x) mean(x, 'omitnan'), medVMPostPokeRight)
% meanPreVMLeft = cellfun(@mean, medVMPrePokeLeft);
meanPreVMLeft = cellfun(@(x) mean(x, 'omitnan'), medVMPrePokeLeft)
% meanPokeVMLeft = cellfun(@mean, medVMDuringPokeLeft);
meanPokeVMLeft = cellfun(@(x) mean(x, 'omitnan'), medVMDuringPokeLeft)
% meanPostVMLeft = cellfun(@mean, medVMPostPokeLeft);
meanPostVMLeft = cellfun(@(x) mean(x, 'omitnan'), medVMPostPokeLeft)

%Statistics
[pWilPreVsPokeVMRight] = signrank(meanPreVMRight, meanPokeVMRight);
[pWilPostVsPokeVMRight] = signrank(meanPostVMRight, meanPokeVMRight);

[pWilPreVsPokeVMLeft] = signrank(meanPreVMLeft, meanPokeVMLeft);
[pWilPostVsPokeVMLeft] = signrank(meanPostVMLeft, meanPokeVMLeft);

%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean medVM pre, during 1st  500 ms and post stimulus -- RIGHT')
t = tiledlayout('flow');
title(t, 'median VM')
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreVMRight; meanPokeVMRight; meanPostVMRight];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,-38.5, ['p=' num2str(pWilPreVsPokeVMRight)])
text(3,-38.5, ['p=' num2str(pWilPostVsPokeVMRight)])
ylabel('median VM (mV)')

% subtitle("Excluded trials: " + num2str(nExcludedRight)+"/"+ num2str(nTotalRight));


nexttile
title("LEFT" + newline + "N =" + num2str(length(left_antenna)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreVMLeft; meanPokeVMLeft; meanPostVMLeft];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,-38.5, ['p=' num2str(pWilPreVsPokeVMLeft)])
text(3,-38.5, ['p=' num2str(pWilPostVsPokeVMLeft)])
ylabel('median VM (mV)')

% subtitle("Excluded trials: " + num2str(nExcludedLeft)+"/"+ num2str(nTotalLeft));


%Pool data from right and left antenna
meanPreVMPooled = [meanPreVMRight meanPreVMLeft];
meanPokeVMPooled = [meanPokeVMRight meanPokeVMLeft];
meanPostVMPooled = [meanPostVMRight meanPostVMLeft];

[pWilPreVsPokeVMPooled] = signrank(meanPreVMPooled, meanPokeVMPooled);
[pWilPostVsPokeVMPooled] = signrank(meanPostVMPooled, meanPokeVMPooled);

nexttile
title("POOLED" + newline + "N =" + num2str(length(meanPreVMPooled)))

groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreVMPooled; meanPokeVMPooled; meanPostVMPooled];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
groups = [1,2,3];
%Use scatter to plot filled circles of each mean

scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,-38.5, ['p=' num2str(pWilPreVsPokeVMPooled)])
text(3,-38.5, ['p=' num2str(pWilPostVsPokeVMPooled)])
ylabel('median VM (mV)')

% print(['medVM_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['medVM_boxplot' '.png'], '-dpng','-r300', '-vector')

set(gcf,'position',[400, 100, 400,1000])
%% Calculate mean X velocity and create boxplots

%Size of poking analysis window in s
analysis_window = windowSizeinSBehavior * sampling_rate;

xVelocityPrePokeRight = [];
xVelocityDuringPokeRight = [];
xVelocityPostPokeRight = [];

for i = 1 : length(right_antenna)
    for k = 1 : numPokes

        % if to_be_removed_ID{1}(i,k) == 1
        %     xVelocityPrePokeRight{i}(k) = NaN;
        %     xVelocityDuringPokeRight{i}(k) = NaN;
        %     xVelocityPostPokeRight{i}(k) = NaN;
        %     continue
        % else

            windowOnRight = right_antenna(i).trigger_on(k);
            windowOffRight = right_antenna(i).trigger_off(k);
            %Sum spikes 1 s pre poke
            xVelocityPrePokeRight{i}(k) = mean(right_antenna(i).xveloc_in_mm...
                (windowOnRight-analysis_window:windowOnRight),'omitnan');
            %Sum spikes 1st second of poke
            xVelocityDuringPokeRight{i}(k) = mean(right_antenna(i).xveloc_in_mm...
                (windowOnRight:windowOnRight+analysis_window),'omitnan');
            %Sum spikes after poke
            xVelocityPostPokeRight{i}(k)  = mean(right_antenna(i).xveloc_in_mm...
                (windowOffRight:windowOffRight+analysis_window),'omitnan');

        end
    end
% end

xVelocityPrePokeLeft = [];
xVelocityDuringPokeLeft = [];
xVelocityPostPokeLeft = [];

for i = 1 : length(left_antenna)
    for k = 1 : numPokes

        % if to_be_removed_ID{2}(i,k) == 1
        %     xVelocityPrePokeLeft{i}(k) = NaN;
        %     xVelocityDuringPokeLeft{i}(k) = NaN;
        %     xVelocityPostPokeLeft{i}(k) = NaN;
        %     continue
        % else

            windowOnLeft = left_antenna(i).trigger_on(k);
            windowOffLeft = left_antenna(i).trigger_off(k);
            %Sum spikes 1 s pre poke
            xVelocityPrePokeLeft{i}(k) = mean(left_antenna(i).xveloc_in_mm...
                (windowOnLeft-analysis_window:windowOnLeft),'omitnan');
            %Sum spikes 1st second of poke
            xVelocityDuringPokeLeft{i}(k) = mean(left_antenna(i).xveloc_in_mm...
                (windowOnLeft:windowOnLeft+analysis_window),'omitnan');
            %Sum spikes after poke
            xVelocityPostPokeLeft{i}(k)  = mean(left_antenna(i).xveloc_in_mm...
                (windowOffLeft:windowOffLeft+analysis_window),'omitnan');
        end
    end
% end

%Calculate means of X velocity
% meanPreXRight = cellfun(@mean, xVelocityPrePokeRight);
meanPreXRight = cellfun(@(x) mean(x, 'omitnan'), xVelocityPrePokeRight)
% meanPokeXRight = cellfun(@mean, xVelocityDuringPokeRight);
meanPokeXRight = cellfun(@(x) mean(x, 'omitnan'), xVelocityDuringPokeRight)
% meanPostXRight = cellfun(@mean, xVelocityPostPokeRight);
meanPostXRight = cellfun(@(x) mean(x, 'omitnan'), xVelocityPostPokeRight)
% meanPreXLeft = cellfun(@mean, xVelocityPrePokeLeft);
meanPreXLeft = cellfun(@(x) mean(x, 'omitnan'), xVelocityPrePokeLeft)
% meanPokeXLeft = cellfun(@mean, xVelocityDuringPokeLeft);
meanPokeXLeft = cellfun(@(x) mean(x, 'omitnan'), xVelocityDuringPokeLeft)
% meanPostXLeft = cellfun(@mean, xVelocityPostPokeLeft);
meanPostXLeft = cellfun(@(x) mean(x, 'omitnan'), xVelocityPostPokeLeft)

%Statistics
[pWilPreVsPokeXRight] = signrank(meanPreXRight, meanPokeXRight);
[pWilPostVsPokeXRight] = signrank(meanPostXRight, meanPokeXRight);

[pWilPreVsPokeXLeft] = signrank(meanPreXLeft, meanPokeXLeft);
[pWilPostVsPokeXLeft] = signrank(meanPostXLeft, meanPokeXLeft);

%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean X Velocity pre, during 1st  500 ms and post stimulus -- RIGHT')
t = tiledlayout('flow');
title(t, 'X Velocity')
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna)))
%Labels for groups
groupNames = ["Pre"; "Poke " + num2str(windowSizeinSBehavior) + "sec"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXRight; meanPokeXRight; meanPostXRight];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
text(1,0.4, ['p=' num2str(round(pWilPreVsPokeXRight,4))])
text(3,0.4, ['p=' num2str(round(pWilPostVsPokeXRight,4))])
ylabel('Forward Velocity (mm/s)')
ylim([-1 1])

%print median
medianXVelocityRightPre = h.med(1).YData(1)
medianXVelocityRightPoke = h.med(2).YData(1)
medianXVelocityRightPost = h.med(3).YData(1)

% subtitle("Excluded trials: " + num2str(nExcludedRight)+"/"+ num2str(nTotalRight));

nexttile
title("LEFT" + newline + "N =" + num2str(length(left_antenna)))
%Labels for groups
groupNames = ["Pre"; "Poke " + num2str(windowSizeinSBehavior) + "sec"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXLeft; meanPokeXLeft; meanPostXLeft];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
text(1,0.4, ['p=' num2str(round(pWilPreVsPokeXLeft,4))])
text(3,0.4, ['p=' num2str(round(pWilPostVsPokeXLeft,4))])
ylabel('Forward Velocity (mm/s)')
ylim([-1 1])

%print median
medianXVelocityLeftPre = h.med(1).YData(1)
medianXVelocityLeftPoke = h.med(2).YData(1)
medianXVelocityLeftPost = h.med(3).YData(1)

%Pool data from right and left antenna
meanPreXPooled = [meanPreXRight meanPreXLeft];
meanPokeXPooled = [meanPokeXRight meanPokeXLeft];
meanPostXPooled = [meanPostXRight meanPostXLeft];

[pWilPreVsPokeXPooled] = signrank(meanPreXPooled, meanPokeXPooled)
[pWilPostVsPokeXPooled] = signrank(meanPostXPooled, meanPokeXPooled)

% subtitle("Excluded trials: " + num2str(nExcludedLeft)+"/"+ num2str(nTotalLeft));



nexttile
title("POOLED" + newline + "N =" + num2str(length(meanPreXPooled)))

groups = [1; 2; 3];
plot_y = [meanPreXPooled; meanPokeXPooled; meanPostXPooled];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
text(1,0.4, ['p=' num2str(round(pWilPreVsPokeXPooled,4))])
text(3,0.4, ['p=' num2str(round(pWilPostVsPokeXPooled,4))])
ylabel('Forward Velocity (mm/s)')
ylim([-1 1])

%print median
medianXVelocityPooledPre = h.med(1).YData(1)
medianXVelocityPooledPoke = h.med(2).YData(1)
medianXVelocityPooledPost = h.med(3).YData(1)

set(gcf,'position',[400, 100, 400,1000])

% print(['X_boxplot_500ms' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['X_boxplot_500ms' '.png'], '-dpng','-r300', '-vector')
%% pool pre and post from left and right and baseline subtract for Forward velocity

preRight = [];
preLeft = [];
prePooled = [];
postPooled = [];
PokeRightXBaselineSub = [];
PokeLeftXBaselineSub = [];
PrePokeBaselineSub = [];
PostPokeBaselineSub = [];
PrePostPooled = [];
PostRightXBaselineSub = [];
PostLeftXBaselineSub = [];

for i = 1 : numCells

preLeft(i) = mean(xVelocityPrePokeLeft{i},'omitnan');
preRight(i) = mean(xVelocityPrePokeRight{i},'omitnan');
prePooled(i) = mean([xVelocityPrePokeRight{i} xVelocityPrePokeLeft{i}],'omitnan');
postPooled(i) = mean([xVelocityPostPokeRight{i} xVelocityPostPokeLeft{i}],'omitnan');
PrePostPooled(i) = mean([xVelocityPrePokeRight{i} xVelocityPrePokeLeft{i}...
    xVelocityPostPokeRight{i} xVelocityPostPokeLeft{i}],'omitnan');

end

% Subtract left baseline from left and right from right
PokeRightXBaselineSub = meanPokeXRight - preRight;
PokeLeftXBaselineSub = meanPokeXLeft - preLeft;
PostPokeBaselineSubRight = meanPostXRight - preRight;
PostPokeBaselineSubLeft = meanPostXLeft - preLeft;

[pWilPreVsPokeLeftBaselineSub] = signrank(preLeft-preLeft, PokeLeftXBaselineSub);
[pWilPreVsPokeRightBaselineSub] = signrank(preRight-preRight, PokeRightXBaselineSub);
[pWilPostVsPokeLeftBaselineSub] = signrank(PostPokeBaselineSubLeft, PokeLeftXBaselineSub);
[pWilPostVsPokeRightBaselineSub] = signrank(PostPokeBaselineSubRight, PokeRightXBaselineSub);

figure
t = tiledlayout('flow');
title(t, "Forward Velocity Pre, During and Post Poke" + newline + "Baseline Subtracted")
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna)))
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [preRight-preRight; PokeRightXBaselineSub; PostPokeBaselineSubRight];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels(groupNames);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
ylabel('Forward Velocity (mm/s)')
ylim([-1 0.4])
text(1,0.3, ['Pre vs Right' newline 'p=' num2str(round(pWilPreVsPokeRightBaselineSub,4))])
text(3,0.3, ['Post vs Right' newline 'p=' num2str(round(pWilPostVsPokeRightBaselineSub,4))])
set(gcf,'position',[400, 100, 400,1000])

nexttile
title("LEFT" + newline + "N =" + num2str(length(right_antenna)))
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [preLeft-preLeft; PokeLeftXBaselineSub; PostPokeBaselineSubLeft];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels(groupNames);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
ylabel('Forward Velocity (mm/s)')
ylim([-1 0.4])
text(1,0.3, ['Pre vs Left' newline 'p=' num2str(round(pWilPreVsPokeLeftBaselineSub,4))])
text(3,0.3, ['Post vs Left' newline 'p=' num2str(round(pWilPostVsPokeLeftBaselineSub,4))])


PokeRightXBaselineSub = [];
PokeLeftXBaselineSub = [];
PrePokeBaselineSub = [];
PostPokeBaselineSub = [];

% Subtract baseline from pre, poke and post
PokeRightXBaselineSub = meanPokeXRight- PrePostPooled;
PokeLeftXBaselineSub = meanPokeXLeft - PrePostPooled;
PrePokeBaselineSub = prePooled - PrePostPooled;
PostPokeBaselineSub = postPooled - PrePostPooled;

[pWilPreVsPokeLeftBaselineSub] = signrank(PrePokeBaselineSub, PokeLeftXBaselineSub);
[pWilPreVsPokeRightBaselineSub] = signrank(PrePokeBaselineSub, PokeRightXBaselineSub);
[pWilPostVsPokeLeftBaselineSub] = signrank(PostPokeBaselineSub, PokeLeftXBaselineSub);
[pWilPostVsPokeRightBaselineSub] = signrank(PostPokeBaselineSub, PokeRightXBaselineSub);

figure
title("Forward Velocity" + newline + "Baseline Pre+Post subtracted")
groupNames = ["Pre"; "Left"; "Right"; "Post"];
groups = [1; 2; 3; 4;];
plot_y = [PrePokeBaselineSub; PokeLeftXBaselineSub; PokeRightXBaselineSub; PostPokeBaselineSub];
h = boxplot2(plot_y',groups);
xticks([1,2,3,4]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
line([1 2], [0.2 0.2])
line([1 2], [0.3 0.3])
line([3 4], [0.2 0.2])
line([3 4], [0.3 0.3])
ylim([-0.9 0.5])

text(1,0.23, ['Pre vs Left' newline 'p=' num2str(round(pWilPreVsPokeLeftBaselineSub,3))])
text(3,0.23, ['Post vs Left' newline 'p=' num2str(round(pWilPostVsPokeLeftBaselineSub,3))])
text(1,0.33, ['Pre vs Right' newline 'p=' num2str(round(pWilPreVsPokeRightBaselineSub,3))])
text(3,0.33, ['Post vs Right' newline 'p=' num2str(round(pWilPostVsPokeRightBaselineSub,3))])
ylabel('Forward Velocity (mm/s)')

set(gcf,'position',[400, 100, 400,1000])

%% Calculate mean Z velocity and create boxplots

%Size of poking analysis window in s
analysis_window = windowSizeinSBehavior * sampling_rate;


for i = 1 : length(right_antenna)
    for k = 1 : numPokes

        % if to_be_removed_ID{1}(i,k) == 1
        %     zVelocityPrePokeRight{i}(k) = NaN;
        %     zVelocityDuringPokeRight{i}(k) = NaN;
        %     zVelocityPostPokeRight{i}(k) = NaN;
        %     continue
        % else

            windowOnRight = right_antenna(i).trigger_on(k);
            windowOffRight = right_antenna(i).trigger_off(k);
            windowOnLeft = left_antenna(i).trigger_on(k);
            windowOffLeft = left_antenna(i).trigger_off(k);
            %Sum spikes 1 s pre poke
            zVelocityPrePokeRight{i}(k) = mean(right_antenna(i).zveloc_in_degree_per_s...
                (windowOnRight-analysis_window:windowOnRight),'omitnan');
            %Sum spikes 1st second of poke
            zVelocityDuringPokeRight{i}(k) = mean(right_antenna(i).zveloc_in_degree_per_s...
                (windowOnRight:windowOnRight+analysis_window),'omitnan');
            %Sum spikes after poke
            zVelocityPostPokeRight{i}(k)  = mean(right_antenna(i).zveloc_in_degree_per_s...
                (windowOffRight:windowOffRight+analysis_window),'omitnan');

        end
    end
% end

for i = 1 : length(left_antenna)
    for k = 1 : numPokes

        % if to_be_removed_ID{2}(i,k) == 1
        %     zVelocityPrePokeLeft{i}(k) = NaN;
        %     zVelocityDuringPokeLeft{i}(k) = NaN;
        %     zVelocityPostPokeLeft{i}(k) = NaN;
        %     continue
        % else

            windowOnLeft = left_antenna(i).trigger_on(k);
            windowOffLeft = left_antenna(i).trigger_off(k);
            %Sum spikes 1 s pre poke
            zVelocityPrePokeLeft{i}(k) = mean(left_antenna(i).zveloc_in_degree_per_s...
                (windowOnLeft-analysis_window:windowOnLeft),'omitnan');
            %Sum spikes 1st second of poke;
            zVelocityDuringPokeLeft{i}(k) = mean(left_antenna(i).zveloc_in_degree_per_s...
                (windowOnLeft:windowOnLeft+analysis_window),'omitnan');
            %Sum spikes after poke
            zVelocityPostPokeLeft{i}(k)  = mean(left_antenna(i).zveloc_in_degree_per_s...
                (windowOffLeft:windowOffLeft+analysis_window),'omitnan');

        end
    end
% end

%Calculate means of Z velocity
% meanPreZRight = cellfun(@mean, zVelocityPrePokeRight);
meanPreZRight = cellfun(@(x) mean(x, 'omitnan'), zVelocityPrePokeRight);
% meanPokeZRight = cellfun(@mean, zVelocityDuringPokeRight);
meanPokeZRight = cellfun(@(x) mean(x, 'omitnan'), zVelocityDuringPokeRight);
% meanPostZRight = cellfun(@mean, zVelocityPostPokeRight);
meanPostZRight = cellfun(@(x) mean(x, 'omitnan'), zVelocityPostPokeRight);
% meanPreZLeft = cellfun(@mean, zVelocityPrePokeLeft);
meanPreZLeft = cellfun(@(x) mean(x, 'omitnan'), zVelocityPrePokeLeft);
% meanPokeZLeft = cellfun(@mean, zVelocityDuringPokeLeft);
meanPokeZLeft = cellfun(@(x) mean(x, 'omitnan'), zVelocityDuringPokeLeft);
% meanPostZLeft = cellfun(@mean, zVelocityPostPokeLeft);
meanPostZLeft = cellfun(@(x) mean(x, 'omitnan'), zVelocityPostPokeLeft);

%Statistics
[pWilPreVsPokeZRight] = signrank(meanPreZRight, meanPokeZRight);
[pWilPostVsPokeZRight] = signrank(meanPostZRight, meanPokeZRight);

[pWilPreVsPokeZLeft] = signrank(meanPreZLeft, meanPokeZLeft);
[pWilPostVsPokeZLeft] = signrank(meanPostZLeft, meanPokeZLeft);

%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean Z velocity pre, during 1st  500 ms and post stimulus -- RIGHT')
t = tiledlayout('flow');
title(t, 'Angular Velocity')
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna)))
%Labels for groups
groupNames = ["Pre"; "Poke " + num2str(windowSizeinSBehavior) + "sec"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreZRight; meanPokeZRight; meanPostZRight];
h = boxplot2(plot_y',groups);
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
text(1,25, ['p=' num2str(round(pWilPreVsPokeZRight,4))])
text(3,25, ['p=' num2str(round(pWilPostVsPokeZRight,4))])
ylabel('Angular Velocity (°/s)')
ylim([-30 30])
%Print median
medianZVelocityRightPre = h.med(1).YData(1)
medianZVelocityRightPoke = h.med(2).YData(1)
medianZVelocityRightPost = h.med(3).YData(1)

% subtitle("Excluded trials: " + num2str(nExcludedRight)+"/"+ num2str(nTotalRight));


nexttile
title("LEFT" + newline + "N =" + num2str(length(left_antenna)))
%Labels for groups
plot_y = [meanPreZLeft; meanPokeZLeft; meanPostZLeft];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
plot_y = [meanPreZLeft; meanPokeZLeft; meanPostZLeft];
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
text(1,25, ['p=' num2str(round(pWilPreVsPokeZLeft,4))])
text(3,25, ['p=' num2str(round(pWilPostVsPokeZLeft,4))])
ylabel('Angular Velocity (°/s)')
ylim([-30 30])
%Print median
medianZVelocityLeftPre = h.med(1).YData(1)
medianZVelocityLeftPoke = h.med(2).YData(1)
medianZVelocityLeftPost = h.med(3).YData(1)

% subtitle("Excluded trials: " + num2str(nExcludedRight)+"/"+ num2str(nTotalRight));

%-----------------------------------------------------
%---------RIGHT ANTENNA IS MIRRORED -> *(-1)----------
%-----------------------------------------------------

%Pool and mirror z velocity data for right and left antenna
meanPreZPooled = [meanPreZLeft meanPreZRight*(-1)];
meanPokeZPooled = [meanPokeZLeft meanPokeZRight*(-1)];
meanPostZPooled = [meanPostZLeft meanPostZRight*(-1)];

% meanPreZPooled = [meanPreZLeft abs(meanPreZRight)];
% meanPokeZPooled = [meanPokeZLeft abs(meanPokeZRight)];
% meanPostZPooled = [meanPostZLeft abs(meanPostZRight)];

[pWilPreVsPokeZPooled] = signrank(meanPreZPooled, meanPokeZPooled);
[pWilPostVsPokeZPooled] = signrank(meanPostZPooled, meanPokeZPooled);

nexttile
title("Pooled" + newline + "N =" + num2str(length(meanPreZPooled)))
%Labels for groups
plot_y = [meanPreZPooled; meanPokeZPooled; meanPostZPooled];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y)
set(s, 'Marker', {'none'});
line(groups, plot_y)
text(1,25, ['p=' num2str(round(pWilPreVsPokeZPooled,4))])
text(3,25, ['p=' num2str(round(pWilPostVsPokeZPooled,4))])
ylabel('Angular Velocity (°/s)')
ylim([-30 30])
%Print median
medianZVelocityPooledPre = h.med(1).YData(1)
medianZVelocityPooledPoke = h.med(2).YData(1)
medianZVelocityPooledPost = h.med(3).YData(1)

set(gcf,'position',[400, 100, 400,1000])
%
% print(['Z_boxplot_500ms' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['Z_boxplot_500ms' '.png'], '-dpng','-r300', '-vector')

%% pool pre and post from left and right and baseline subtract for Angular velocity
preRight = [];
preLeft = [];
prePooled = [];
postPooled = [];
PokeRightZBaselineSub = [];
PokeLeftZBaselineSub = [];
PrePokeBaselineSub = [];
PostPokeBaselineSub = [];
PrePostPooled = [];
PostRightZBaselineSub = [];
PostLeftZBaselineSub = [];

for i = 1 : numCells

preLeft(i) = mean(zVelocityPrePokeLeft{i},'omitnan');
preRight(i) = mean(zVelocityPrePokeRight{i},'omitnan');
prePooled(i) = mean([zVelocityPrePokeRight{i} zVelocityPrePokeLeft{i}],'omitnan');
postPooled(i) = mean([zVelocityPostPokeRight{i} zVelocityPostPokeLeft{i}],'omitnan');
PrePostPooled(i) = mean([zVelocityPrePokeRight{i} zVelocityPrePokeLeft{i}...
    zVelocityPostPokeRight{i} zVelocityPostPokeLeft{i}],'omitnan');

end

% Subtract left baseline from left and right from right
PokeRightZBaselineSub = meanPokeZRight - preRight;
PokeLeftZBaselineSub = meanPokeZLeft - preLeft;
PostPokeBaselineSubRight = meanPostZRight - preRight;
PostPokeBaselineSubLeft = meanPostZLeft - preLeft;

[pWilPreVsPokeLeftBaselineSub] = signrank(preLeft-preLeft, PokeLeftZBaselineSub);
[pWilPreVsPokeRightBaselineSub] = signrank(preRight-preRight, PokeRightZBaselineSub);
[pWilPostVsPokeLeftBaselineSub] = signrank(PostPokeBaselineSubLeft, PokeLeftZBaselineSub);
[pWilPostVsPokeRightBaselineSub] = signrank(PostPokeBaselineSubRight, PokeRightZBaselineSub);

figure
t = tiledlayout('flow');
title(t, "Angular Velocity Pre, During and Post Poke" + newline + "Baseline Subtracted")
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna)))
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [preRight-preRight; PokeRightZBaselineSub; PostPokeBaselineSubRight];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels(groupNames);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
ylabel('Angular Velocity (°/s)')
% ylim([-1 0.4])
text(1,0.3, ['Pre vs Right' newline 'p=' num2str(round(pWilPreVsPokeRightBaselineSub,4))])
text(3,0.3, ['Post vs Right' newline 'p=' num2str(round(pWilPostVsPokeRightBaselineSub,4))])
set(gcf,'position',[400, 100, 400,1000])

nexttile
title("LEFT" + newline + "N =" + num2str(length(right_antenna)))
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [preLeft-preLeft; PokeLeftZBaselineSub; PostPokeBaselineSubLeft];
h = boxplot2(plot_y',groups,'whisker', 0);
set(h.out, 'Marker', 'none')
xticks([1,2,3]);
xticklabels(groupNames);
hold on
%Use scatter to plot filled circles of each mean
s = scatter(groups, plot_y);
set(s, 'Marker', {'none'});
line(groups, plot_y)
ylabel('Angular Velocity (°/s)')
% ylim([-1 0.4])
text(1,0.3, ['Pre vs Left' newline 'p=' num2str(round(pWilPreVsPokeLeftBaselineSub,4))])
text(3,0.3, ['Post vs Left' newline 'p=' num2str(round(pWilPostVsPokeLeftBaselineSub,4))])


PokeRightZBaselineSub = [];
PokeLeftZBaselineSub = [];
PrePokeBaselineSub = [];
PostPokeBaselineSub = [];

% Subtract baseline from pre, poke and post
PokeRightZBaselineSub = meanPokeZRight- PrePostPooled;
PokeLeftZBaselineSub = meanPokeZLeft - PrePostPooled;
PrePokeBaselineSub = prePooled - PrePostPooled;
PostPokeBaselineSub = postPooled - PrePostPooled;

[pWilPreVsPokeLeftBaselineSub] = signrank(PrePokeBaselineSub, PokeLeftZBaselineSub);
[pWilPreVsPokeRightBaselineSub] = signrank(PrePokeBaselineSub, PokeRightZBaselineSub);
[pWilPostVsPokeLeftBaselineSub] = signrank(PostPokeBaselineSub, PokeLeftZBaselineSub);
[pWilPostVsPokeRightBaselineSub] = signrank(PostPokeBaselineSub, PokeRightZBaselineSub);

figure
title(" Angular Velocity" + newline + "Baseline Pre+Post subtracted")
groupNames = ["Pre"; "Left"; "Right"; "Post"];
groups = [1; 2; 3; 4;];
plot_y = [PrePokeBaselineSub; PokeLeftZBaselineSub; PokeRightZBaselineSub; PostPokeBaselineSub];
h = boxplot2(plot_y',groups);
xticks([1,2,3,4]);
xticklabels([groupNames]);
hold on
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
line([1 2], [0.2 0.2])
line([1 2], [0.3 0.3])
line([3 4], [0.2 0.2])
line([3 4], [0.3 0.3])
% ylim([-0.9 0.5])

text(1,0.23, ['Pre vs Left' newline 'p=' num2str(round(pWilPreVsPokeLeftBaselineSub,3))])
text(3,0.23, ['Post vs Left' newline 'p=' num2str(round(pWilPostVsPokeLeftBaselineSub,3))])
text(1,0.33, ['Pre vs Right' newline 'p=' num2str(round(pWilPreVsPokeRightBaselineSub,3))])
text(3,0.33, ['Post vs Right' newline 'p=' num2str(round(pWilPostVsPokeRightBaselineSub,3))])
ylabel('Angular Velocity (s/°)')

set(gcf,'position',[400, 100, 400,1000])
%% Compare left and right Spike Rate, Forward and Angular Velocity

figure

% Colormap for scatter plots
color_group =  gray(9);

nexttile
groups = [1; 2];
groupNames = ["Left"; "Right"];
plot_y = [PokeLeftSpikeBaselineSub; PokeRightSpikeBaselineSub];
boxplot2(plot_y', groups)
xticks([1,2]);
xticklabels(groupNames);
hold on
scatter(groups, plot_y, [], color_group, 'filled')
line(groups, plot_y)
title("Spike Rate")
ylim([-5 15])
ylabel("Spike Rate (Hz)")
hold on
line([1 2], [13 13])
% Statistic
[~, pval] = ttest(PokeLeftSpikeBaselineSub, PokeRightSpikeBaselineSub);
text(1.25,14, ['p=' num2str(pval)])

nexttile
plot_y = [PokeLeftXBaselineSub; PokeRightXBaselineSub];
boxplot2(plot_y', groups)
hold on 
xticks([1,2]);
xticklabels(groupNames);
scatter(groups, plot_y, [], color_group, 'filled')
line(groups, plot_y)
title("Forward Velocity")
ylim([-1 0.2])
ylabel("Forward Velocity (mm/s)")
hold on
line([1 2], [0.1 0.1])
% Statistic
[~, pval] = ttest(PokeLeftXBaselineSub, PokeRightXBaselineSub);
text(1.25,0.15, ['p=' num2str(pval)])

nexttile
plot_y = [PokeLeftZBaselineSub; PokeRightZBaselineSub];
boxplot2(plot_y', groups)
hold on 
xticks([1,2]);
xticklabels(groupNames);
scatter(groups, plot_y, [], color_group, 'filled')
line(groups, plot_y)
title("Angular Velocity")
ylim([-35 35])
ylabel("Angular Velocity (°/s)")
hold on
line([1 2], [32 32])
% Statistic
[~, pval] = ttest(PokeLeftZBaselineSub, PokeRightZBaselineSub);
text(1.25,34, ['p=' num2str(pval)])

set(gcf,'position',[400, 100, 400,1000])

%% Save statistics % median

T = table(pWilPreVsPokeSpikesRight, pWilPostVsPokeSpikesRight,...
    pWilPreVsPokeSpikesLeft, pWilPostVsPokeSpikesLeft,...
    pWilPreVsPokeXRight, pWilPostVsPokeXRight,...
    pWilPreVsPokeXLeft, pWilPostVsPokeXLeft,...
    pWilPreVsPokeZRight, pWilPostVsPokeZRight,...
    pWilPreVsPokeZLeft, pWilPostVsPokeZLeft,...
    'VariableNames', { 'Pre vs Poke Spikes -- Right', 'Poke vs Post Spikes -- Right',...
    'Pre vs Poke Spikes -- Left', 'Poke vs Post Spikes -- Left',...
    'Pre vs Poke X -- Right', 'Poke vs Post X -- Right',...
    'Pre vs Poke X -- Left', 'Poke vs Post X -- Left',...
    'Pre vs Poke Z -- Right', 'Poke vs Post Z -- Right',...
    'Pre vs Poke Z -- Left', 'Poke vs Post Z -- Left'} );
% Write data to text file
writetable(T, 'Stats.xls')

%Data collection for median table
medianFiringRate = [medianFiringRateRightPre; medianFiringRateRightPoke; medianFiringRateRightPost;
    medianFiringRateLeftPre; medianFiringRateLeftPoke; medianFiringRateLeftPost;
    medianFiringRatePooledPre; medianFiringRatePooledPoke; medianFiringRatePooledPost];
FiringRate = {'Median Firing Rate Right Pre'; 'Median Firing Rate Right Poke'; 'Median Firing Rate Right Post';
    'Median Firing Rate Left Pre'; 'Median Firing Rate Left Poke'; 'Median Firing Rate Left Post';
    'Median Firing Rate Pooled Pre'; 'Median Firing Rate Pooled Poke'; 'Median Firing Rate Pooled Post'};

medianXVelocity = [medianXVelocityRightPre; medianXVelocityRightPoke; medianXVelocityRightPost;
    medianXVelocityLeftPre; medianXVelocityLeftPoke; medianXVelocityLeftPost;
    medianXVelocityPooledPre; medianXVelocityPooledPoke; medianXVelocityPooledPost];
XVelocity = {'Median XVelocity Right Pre'; 'Median XVelocity Right Poke'; 'Median XVelocity Right Post';
    'Median XVelocity Left Pre'; 'Median XVelocity Left Poke'; 'Median XVelocity Left Post';
    'Median XVelocity Pooled Pre'; 'Median XVelocity Pooled Poke'; 'Median XVelocity Pooled Post'};

medianZVelocity = [medianZVelocityRightPre; medianZVelocityRightPoke; medianZVelocityRightPost;
    medianZVelocityLeftPre; medianZVelocityLeftPoke; medianZVelocityLeftPost;
    medianZVelocityPooledPre; medianZVelocityPooledPoke; medianZVelocityPooledPost];
ZVelocity = {'Median ZVelocity Right Pre'; 'Median ZVelocity Right Poke'; 'Median ZVelocity Right Post';
    'Median ZVelocity Left Pre'; 'Median ZVelocity Left Poke'; 'Median ZVelocity Left Post';
    'Median ZVelocity Pooled Pre'; 'Median ZVelocity Pooled Poke'; 'Median ZVelocity Pooled Post'};

T = table(FiringRate, medianFiringRate, XVelocity, medianXVelocity, ZVelocity, medianZVelocity)

% Write data to text file
writetable(T, 'Median.xls')
%% Analyse 1 s window after stim and compare firing rate with X and Z respectively --- LEFT

for i = 1 : length(left_antenna)

    for k = 1 : numPokes

        %Window for analyzing firing rate, x and z velocity, 1 sec
        window_on = left_antenna(i).trigger_on(k);
        window_off = left_antenna(i).trigger_on(k)+1*sampling_rate;

        %Add all data to an array
        StimAnalysisSpikesLeft{i,k} = left_antenna(i).spikes(window_on:window_off);
        StimAnalysisXvelocLeft{i,k} = left_antenna(i).xveloc_in_mm(window_on:window_off);
        StimAnalysisZvelocLeft{i,k} = left_antenna(i).zveloc_in_mm(window_on:window_off);

    end

end

%Sum up the spikes in the 1s analysis window
StimAnalysisSumSpikesLeft = cellfun(@(x) sum(x, 'all', 'omitnan'), StimAnalysisSpikesLeft);
%Either take the mean or the minimum (maximal backward) x velocity
StimAnalysisMeanXvelocLeft = cellfun(@(x) mean(x, 'all', 'omitnan'), StimAnalysisXvelocLeft);
StimAnalysisMinXvelocLeft = cellfun(@(x) min(x), StimAnalysisXvelocLeft);
%Either take the mean or the max(left antenna) z velocity
StimAnalysisMeanZvelocLeft = cellfun(@(x) mean(x, 'all', 'omitnan'), StimAnalysisZvelocLeft);
StimAnalysisMaxZvelocLeft = cellfun(@(x) max(x), StimAnalysisXvelocLeft);


%Create a colormap with 10 different colors
% colors = [
%     0.8500 0.3250 0.0980;   % Reddish-Orange
%     0.9290 0.6940 0.1250;   % Yellow
%     0.4940 0.1840 0.5560;   % Purple
%     0.4660 0.6740 0.1880;   % Green
%     0.3010 0.7450 0.9330;   % Light Blue
%     0.6350 0.0780 0.1840;   % Dark Red
%     0.0000 0.4470 0.7410;   % Blue
%     0.8500 0.3250 0.0980;   % Reddish-Orange (again)
%     0.4940 0.1840 0.5560;   % Purple (again)
%     0.9290 0.6940 0.1250;   % Yellow (again)
% ];
%
% colormap(colors);
% cmap = colormap;


%Scatter plot of Firing rate vs X velocity
figure
for n = 1 : numCells
    %Create a new tile for every cell
    nexttile
    title('Spikes vs. X Veloc --- LEFT')
    for i = 1 : numPokes
        %Plot Firing rate vs Mean X Vel
        scatter(StimAnalysisSumSpikesLeft(n,i), StimAnalysisMeanXvelocLeft(n,i),50 ,...
            'filled', 'jitter','on', 'jitterAmount',0.25)
        hold on
        %Plot Firing rate vs Min X Vel
        scatter(StimAnalysisSumSpikesLeft(n,i), StimAnalysisMinXvelocLeft(n,i),50 ,...
            'jitter','on', 'jitterAmount',0.25)

        xlabel('MDN Firing Rate (Hz)')
        ylabel('X Velocity mm/s')
    end

end

%Scatter plot of Firing rate vs Z velocity
figure
for n = 1 : length(left_antenna)
    %Create a new tile for every cell
    nexttile
    title('Spikes vs. Z Veloc --- LEFT')
    for i = 1 : numPokes
        %Plot Firing rate vs Mean Z Vel
        %Add jitter to discriminate between data points with same firing r.
        scatter(StimAnalysisSumSpikesLeft(n,i), StimAnalysisMeanZvelocLeft(n,i),50, ...
            'filled', 'jitter','on', 'jitterAmount',0.25)
        hold on
        %Plot Firing rate vs Max Z Vel
        scatter(StimAnalysisSumSpikesLeft(n,i), StimAnalysisMaxZvelocLeft(n,i),50 ,...
            'jitter', 'on', 'jitterAmount',0.25)
        xlabel('MDN Firing Rate (Hz)')
        ylabel('Z Velocity mm/s')
    end
end


%% Analyse 1 s window after stim and compare firing rate with X and Z respectively --- RIGHT

for i = 1 : numCells

    for k = 1 : numPokes

        window_on = right_antenna(i).trigger_on(k);
        window_off = right_antenna(i).trigger_on(k)+1*sampling_rate;

        StimAnalysisSpikesRight{i,k} = right_antenna(i).spikes(window_on:window_off);
        StimAnalysisXvelocRight{i,k} = right_antenna(i).xveloc_in_mm(window_on:window_off);
        StimAnalysisZvelocRight{i,k} = right_antenna(i).zveloc_in_mm(window_on:window_off);

    end

end
%
% figure
% plot(StimAnalysisSpikesRight{1,1}, 'o')



StimAnalysisSumSpikesRight = cellfun(@(x) sum(x, 'all', 'omitnan'), StimAnalysisSpikesRight);
StimAnalysisMeanXvelocRight = cellfun(@(x) mean(x, 'all', 'omitnan'), StimAnalysisXvelocRight);
StimAnalysisMinXvelocRight = cellfun(@(x) min(x), StimAnalysisXvelocRight);
StimAnalysisMeanZvelocRight = cellfun(@(x) mean(x, 'all', 'omitnan'), StimAnalysisZvelocRight);
StimAnalysisMinZvelocRight = cellfun(@(x) min(x), StimAnalysisXvelocRight);


colors = [
    0.8500 0.3250 0.0980;   % Reddish-Orange
    0.9290 0.6940 0.1250;   % Yellow
    0.4940 0.1840 0.5560;   % Purple
    0.4660 0.6740 0.1880;   % Green
    0.3010 0.7450 0.9330;   % Light Blue
    0.6350 0.0780 0.1840;   % Dark Red
    0.0000 0.4470 0.7410;   % Blue
    0.8500 0.3250 0.0980;   % Reddish-Orange (again)
    0.4940 0.1840 0.5560;   % Purple (again)
    0.9290 0.6940 0.1250;   % Yellow (again)
    ];

colormap(colors);
cmap = colormap;

figure;
for n = 1 : numCells
    %nexttile
    title('name', 'Spikes vs. X Veloc --- RIGHT')
    %figure('name', 'Spikes vs. X Veloc --- RIGHT')
    for i = 1 : numPokes

        scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMeanXvelocRight(n,i),50 ,...
            'filled', 'jitter','on', 'jitterAmount',0.25)
        hold on

        scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMinXvelocRight(n,i),50 ,...
            'jitter','on', 'jitterAmount',0.25)

        xlabel('MDN Firing Rate (Hz)')
        ylabel('X Velocity mm/s')
    end

end

figure;
for n = 1 : numCells
    %nexttile
    title('name', 'Spikes vs. Z Veloc --- RIGHT')
    %figure('name', 'Spikes vs. Z Veloc --- RIGHT')
    for i = 1 : numPokes

        scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMeanZvelocRight(n,i),50, ...
            'filled', 'jitter','on', 'jitterAmount',0.25)
        hold on

        scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMinZvelocRight(n,i),50 ,...
            'jitter', 'on', 'jitterAmount',0.25)
        xlabel('MDN Firing Rate (Hz)')
        ylabel('Z Velocity mm/s')
    end
end

%% Find min X and Z velocity (= max backward locomotion) and get spike max 100 ms before ---- RIGHT
clearvars minXvelocIdxRight StimAnalysisMinXvelocRight beforeStimSpikes...
    beforeStimSpikesX beforeStimSpikesZ

%Find the minimal X velocity in the previous generated 1s analysis window
%and also save the index
[StimAnalysisMinXvelocRight, minXvelocIdxRight] = cellfun(@(x) min(x),...
    StimAnalysisXvelocRight);
%Same for Z velocity
[StimAnalysisMinZvelocRight, minZvelocIdxRight] = cellfun(@(x) min(x),...
    StimAnalysisZvelocRight);

%Find all spikes in a determined window previous to the maximal X or Z
%velocity

%Size of window for spike analysis in ms before max x/z velocity
preMaxInMs = 100;
preMaxFactor = (1000/preMaxInMs);
preMaxWindow = sampling_rate/(1000/preMaxInMs);
%Calculated to data points

for i = 1 : numCells
    for k = 1 : numPokes
        %Anaysis window starts x ms before max x/z velocity
        startX = minXvelocIdxRight(i,k)-preMaxWindow;
        stopX = minXvelocIdxRight(i,k);
        %Break loop if analysis window lies outside of 1s window
        if startX < 0
            break
        end

        %Loop through spiking date and save spikes within the analysis window
        beforeStimSpikesX{i,k} = StimAnalysisSpikesRight{i,k} (startX:stopX,1);

    end
end

%Same as previous loop but for max Z velocities
for i = 1 : numCells
    for k = 1 : numPokes

        startZ = minZvelocIdxRight(i,k)-preMaxWindow;
        stopZ = minZvelocIdxRight(i,k);

        if startZ < 0
            break
        end
        beforeStimSpikesZ{i,k} = StimAnalysisSpikesRight{i,k} (startZ:stopZ,1);

    end
end

% Sum up spikes and calculate frequency
beforeStimSpikesXSum = cellfun(@(x) sum(x, 'all', 'omitnan'), beforeStimSpikesX);
beforeStimSpikesZSum = cellfun(@(x) sum(x, 'all', 'omitnan'), beforeStimSpikesZ);
beforeStimSpikeXHz = beforeStimSpikesXSum*preMaxFactor;
beforeStimSpikeZHz = beforeStimSpikesZSum*preMaxFactor;

% Replace emtpy cells with nans
% beforeStimSpikesX(cellfun('isempty', beforeStimSpikesX)) = {NaN};



figure('name', 'Spikes vs. X Veloc --- RIGHT short time window')
t = tiledlayout("flow");
title(t, preMaxInMs, 'ms')
for n = 1 : numCells
    nexttile
    for i = 1 : numPokes

        scatter(beforeStimSpikeXHz(n,i), StimAnalysisMinXvelocRight(n,i),50 ,...
            'filled', 'jitter','on', 'jitterAmount',0.25)
        hold on
        %
        %         scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMeanXvelocRight(n,i),50 ,...
        %             'jitter','on', 'jitterAmount',0.25)
        %         hold on
        %
        %         scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMinXvelocRight(n,i),50 ,...
        %             's', 'jitter','on', 'jitterAmount',0.25)
        %
        xlabel('MDN Firing Rate (Hz)')
        ylabel('X Velocity mm/s')
    end
    %title(n,'Interpreter','none')
end


figure('name', 'Spikes vs. X Veloc --- RIGHT short time window')

scatter(beforeStimSpikeXHz, StimAnalysisMinXvelocRight,50 ,...
    'filled', 'jitter','on', 'jitterAmount',0.25)

figure('name', 'Spikes vs. Z Veloc --- RIGHT short time window')
t = tiledlayout("flow");
title(t, preMaxInMs, 'ms')
for n = 1 : numCells
    nexttile
    for i = 1 : numPokes

        scatter(beforeStimSpikeZHz(n,i), StimAnalysisMinZvelocRight(n,i),50 ,...
            'filled', 'jitter','on', 'jitterAmount',0.25)
        hold on
        %
        %         scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMeanXvelocRight(n,i),50 ,...
        %             'jitter','on', 'jitterAmount',0.25)
        %         hold on
        %
        %         scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMinXvelocRight(n,i),50 ,...
        %             's', 'jitter','on', 'jitterAmount',0.25)
        %
        xlabel('MDN Firing Rate (Hz)')
        ylabel('Z Velocity mm/s')
    end
    title(n,'Interpreter','none')
end
%% Find min X velocity (= max backward locomotion) and get spike max 100 ms before ---- LEFT
clearvars minXvelocIdxRight StimAnalysisMinXvelocRight beforeStimSpikeHz...
    beforeStimSpikesX beforeStimSpikesZ beforeStimSpikeZHz beforeStimSpikeXHz

[StimAnalysisMinXvelocLeft, minXvelocIdxLeft] = cellfun(@(x) min(x),...
    StimAnalysisXvelocLeft);
[StimAnalysisMaxZvelocLeft, maxZvelocIdxLeft] = cellfun(@(x) max(x),...
    StimAnalysisXvelocLeft);

StimAnalysisSpikesLeft_temp = cellfun(@(c) fillmissing(c,'constant',0),...
    StimAnalysisSpikesLeft,'UniformOutput',false)


%Find all spikes in a determined window previous to the maximal X or Z
%velocity

%Size of window for spike analysis in ms before max x/z velocity
preMaxInMs = 100;
preMaxFactor = (1000/preMaxInMs);
preMaxWindow = sampling_rate/(1000/preMaxInMs);
%Calculated to data points

for i = 1 : numCells
    for k = 1 : numPokes
        %Anaysis window starts x ms before max x/z velocity
        startX = minXvelocIdxLeft(i,k)-preMaxWindow
        stopX = minXvelocIdxLeft(i,k)
        %Break loop if analysis window lies outside of 1s window
        if startX < 0
            break
        end

        %Loop through spiking date and save spikes within the analysis window
        beforeStimSpikesX{i,k} = StimAnalysisSpikesLeft{i,k} (startX:stopX,1);

    end
end

%Same as previous loop but for max Z velocities
for i = 1 : numCells
    for k = 2 : numPokes

        startZ = maxZvelocIdxLeft(i,k)-preMaxWindow;
        stopZ = maxZvelocIdxLeft(i,k);

        if startZ < 0
            break
        end
        beforeStimSpikesZ{i,k} = StimAnalysisSpikesLeft{i,k} (startZ:stopZ,1)

    end
end

beforeStimSpikesXSum = cellfun(@(x) sum(x, 'all', 'omitnan'), beforeStimSpikesX);
beforeStimSpikesZSum = cellfun(@(x) sum(x, 'all', 'omitnan'), beforeStimSpikesZ);
beforeStimSpikeXHz = beforeStimSpikesXSum*10;
beforeStimSpikeZHz = beforeStimSpikesZSum*10;


figure('name', 'Spikes vs. X Veloc --- Left')
t = tiledlayout("flow");
title(t, preMaxInMs, 'ms')
for n = 1 : numCells
    nexttile
    for i = 1 : numPokes

        scatter(beforeStimSpikeXHz(n,i), StimAnalysisMinXvelocLeft(n,i),50 ,...
            'filled', 'jitter','on', 'jitterAmount',0.25)
        hold on
        %
        %         scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMeanXvelocRight(n,i),50 ,...
        %             'jitter','on', 'jitterAmount',0.25)
        %         hold on
        %
        %         scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMinXvelocRight(n,i),50 ,...
        %             's', 'jitter','on', 'jitterAmount',0.25)
        %
        xlabel('MDN Firing Rate (Hz)')
        ylabel('X Velocity mm/s')
    end
    title(n,'Interpreter','none')
end

figure('name', 'Spikes vs. Z Veloc --- Left')
t = tiledlayout("flow");
title(t, preMaxInMs, 'ms')
for n = 1 : numCells
    nexttile
    for i = 1 : numPokes

        scatter(beforeStimSpikeZHz(n,i), StimAnalysisMaxZvelocLeft(n,i),50 ,...
            'filled', 'jitter','on', 'jitterAmount',0.25)
        hold on
        %
        %         scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMeanXvelocRight(n,i),50 ,...
        %             'jitter','on', 'jitterAmount',0.25)
        %         hold on
        %
        %         scatter(StimAnalysisSumSpikesRight(n,i), StimAnalysisMinXvelocRight(n,i),50 ,...
        %             's', 'jitter','on', 'jitterAmount',0.25)
        %
        xlabel('MDN Firing Rate (Hz)')
        ylabel('Z Velocity mm/s')
    end
    title(n,'Interpreter','none')
end


%% Plot all spike events

figure
n = 0;
for i = 1 : length(right_antenna)
    for k = 1 : numPokes
        windowOnRight = right_antenna(i).trigger_on(k);
        windowOffRight = right_antenna(i).trigger_off(k);
        windowOnLeft = left_antenna(i).trigger_on(k);
        windowOffLeft = left_antenna(i).trigger_off(k);
        x = [1:length(right_antenna(i).spikes(windowOnRight-20000:windowOffRight+20000))]/sampling_rate;
        plot(x, right_antenna(i).spikes(windowOnRight-20000:windowOffRight+20000)+n,'m.')
        hold on
        x = [1:length(left_antenna(i).spikes(windowOnLeft-20000:windowOffLeft+20000))]/sampling_rate;
        plot(x, left_antenna(i).spikes(windowOnLeft-20000:windowOffLeft+20000)+n,'k.')
        n = n + 3;
    end
    n = n+5;
end

rectangle('Position',[1 0 2 350], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')

%print(['spiketrains_all' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
%print(['spiketrains_all' '.png'], '-dpng','-r300', '-vector')

%% Cross correlation between spikes and x velocity during poking
clearvars meanTrialspikesXvelocRight meanTrialspikesXvelocLeft lags spikesXvelocAll meanspikesXvelocAll

for k = 1 : length(right_antenna)
    clearvars tempRight tempLeft
    %figure('Name','X Correlation: spikes X xvelocity');

    for i = 1 : numPokes
        %Set window for cross correlation
        windowOnRight = floor(right_antenna(k).trigger_on(i)/binsize);
        windowOffRight = floor(right_antenna(k).trigger_off(i)/binsize);
        windowOnLeft = floor(left_antenna(k).trigger_on(i)/binsize);
        windowOffLeft = floor(left_antenna(k).trigger_off(i)/binsize);

        [tempRight(i,:),lags] = xcorr(zscore(right_antenna(k).spikesbinned(windowOnRight:windowOffRight)),...
            zscore(right_antenna(k).xveloc_in_mm_binned(windowOnRight:windowOffRight)), 10 ,'normalized');

        %plot(lags/(1/binsize_factor),spikesZvelocRight)


        [tempLeft(i,:),lags] = xcorr(zscore(left_antenna(k).spikesbinned(windowOnLeft:windowOffLeft)),...
            zscore(left_antenna(k).xveloc_in_mm_binned(windowOnLeft:windowOffLeft)), 10 ,'normalized');
        %plot(lags/(1/binsize_factor),spikesZvelocLeft)


    end
    meanTrialspikesXvelocRight(k,:) = mean(tempRight,'omitnan');
    meanTrialspikesXvelocLeft(k,:) = mean(tempLeft,'omitnan');

end
meanOverallspikesXvelocRight = mean(meanTrialspikesXvelocRight);
meanOverallspikesXvelocLeft = mean(meanTrialspikesXvelocLeft);
meanOverallspikesXvelocBoth = [meanOverallspikesXvelocRight + meanOverallspikesXvelocLeft]/2;

figure
tiledlayout(3,1)
nexttile
plot(lags/(1/binsize_factor), meanTrialspikesXvelocRight, 'Color', [0.7 0.7 0.7])
hold on
plot(lags/(1/binsize_factor), meanOverallspikesXvelocRight, 'k', 'LineWidth', 3)
title('RIGHT Antenna', 'Cross correlation of Firing Rate vs X Velocity')
xlabel("Lag (s)")
ylabel("Correlation Coefficient")
ylim([-0.5 0.5])

nexttile
plot(lags/(1/binsize_factor), meanTrialspikesXvelocLeft, 'Color', [0.7 0.7 0.7])
hold on
plot(lags/(1/binsize_factor), meanOverallspikesXvelocLeft, 'k', 'LineWidth', 3)
title('LEFT Antenna', 'Cross correlation of Firing Rate vs X Velocity')
xlabel("Lag (s)")
ylabel("Correlation Coefficient")
ylim([-0.5 0.5])

nexttile
plot(lags/(1/binsize_factor),meanOverallspikesXvelocRight)
hold on
plot(lags/(1/binsize_factor),meanOverallspikesXvelocLeft)
hold on
plot(lags/(1/binsize_factor),meanOverallspikesXvelocBoth)

xlabel("Lag (s)")
ylabel("Correlation Coefficient")
title("Cross correlation of Firing Rate vs X Velocity")
legend('Right', 'Left', 'Both','Location','northwest')
set(gcf,'position',[400, 100, 400,1000])
% 
% print(['crosscorr_Spikes_X' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['crosscorr_Spikes_X' '.png'], '-dpng','-r300', '-vector')
%% Cross correlation between spikes and Z velocity during poking
clearvars spikesXveloc lags spikesXvelocAll meanTrialspikesZvelocRight ...
    meanTrialspikesZvelocLeft

for k = 1 : length(right_antenna)
    clearvars tempRight tempLeft
    %figure('Name','X Correlation: spikes Z xvelocity');

    for i = 1 : numPokes

        windowOnRight = floor(right_antenna(k).trigger_on(i)/binsize);
        windowOffRight = floor(right_antenna(k).trigger_off(i)/binsize);
        windowOnLeft = floor(left_antenna(k).trigger_on(i)/binsize);
        windowOffLeft = floor(left_antenna(k).trigger_off(i)/binsize);

        [tempRight(i,:),lags] = xcorr(zscore(right_antenna(k).spikesbinned(windowOnRight:windowOffRight)),...
            zscore(right_antenna(k).zvelocbinned(windowOnRight:windowOffRight))*(-1), 10,'normalized');

        %plot(lags/(1/binsize_factor),spikesZvelocRight)

        [tempLeft(i,:),lags] = xcorr(zscore(left_antenna(k).spikesbinned(windowOnLeft:windowOffLeft)),...
            zscore(left_antenna(k).zvelocbinned(windowOnLeft:windowOffLeft)), 10 ,'normalized');
        %plot(lags/(1/binsize_factor),spikesZvelocLeft)

    end

    meanTrialspikesZvelocRight(k,:) = mean(tempRight,'omitnan');
    meanTrialspikesZvelocLeft(k,:) = mean(tempLeft,'omitnan');

end
meanOverallspikesZvelocRight = mean(meanTrialspikesZvelocRight);
meanOverallspikesZvelocLeft = mean(meanTrialspikesZvelocLeft);
%%% RIGHT SIDE Z VELOCITY IS MIRRORED (*-1)%%%%
meanOverallspikesZvelocBoth = [meanOverallspikesZvelocRight + meanOverallspikesZvelocLeft]/2;

figure
tiledlayout(3,1)
nexttile
plot(lags/(1/binsize_factor), meanTrialspikesZvelocRight, 'Color', [0.7 0.7 0.7])
hold on
plot(lags/(1/binsize_factor), meanOverallspikesZvelocRight, 'k', 'LineWidth', 3)
title('RIGHT Antenna', 'Cross correlation of Firing Rate vs Z Velocity')
xlabel("Lag (s)")
ylabel("Correlation Coefficient")
ylim([-0.5 0.5])

nexttile
plot(lags/(1/binsize_factor), meanTrialspikesZvelocLeft, 'Color', [0.7 0.7 0.7])
hold on
plot(lags/(1/binsize_factor), meanOverallspikesZvelocLeft, 'k', 'LineWidth', 3)
title('LEFT Antenna', 'Cross correlation of Firing Rate vs Z Velocity')
xlabel("Lag (s)")
ylabel("Correlation Coefficient")
ylim([-0.5 0.5])

nexttile
plot(lags/(1/binsize_factor),meanOverallspikesZvelocRight)
hold on
plot(lags/(1/binsize_factor),meanOverallspikesZvelocLeft)
hold on
plot(lags/(1/binsize_factor),meanOverallspikesZvelocBoth)

xlabel("Lag (s)")
ylabel("Correlation Coefficient")
title("Cross correlation of Firing Rate vs Z Velocity")
legend('Right', 'Left', 'Both','Location','northwest')
set(gcf,'position',[400, 100, 400,1000])

% print(['crosscorr_Spikes_Z' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['crosscorr_Spikes_Z' '.png'], '-dpng','-r300', '-vector')
%% Cross correlation between X and Z velocity during poking
clearvars spikesXveloc lags spikesXvelocAll meanspikesXvelocRight meanspikesXvelocLeft

for k = 1 : length(right_antenna)
    clearvars tempRight tempLeft
    %figure('Name','X Correlation: spikes Z xvelocity');

    for i = 1 : numPokes

        windowOnRight = floor(right_antenna(k).trigger_on(i)/binsize);
        windowOffRight = floor(right_antenna(k).trigger_off(i)/binsize);
        windowOnLeft = floor(left_antenna(k).trigger_on(i)/binsize);
        windowOffLeft = floor(left_antenna(k).trigger_off(i)/binsize);

        [tempRight(i,:),lags] = xcorr(zscore(right_antenna(k).xvelocbinned(windowOnRight:windowOffRight)),...
            zscore(right_antenna(k).zvelocbinned(windowOnRight:windowOffRight))*(-1) ,10,'normalized');

        %plot(lags/(1/binsize_factor),spikesZvelocRight)


        [tempLeft(i,:),lags] = xcorr(zscore(left_antenna(k).xvelocbinned(windowOnLeft:windowOffLeft)),...
            zscore(left_antenna(k).zvelocbinned(windowOnLeft:windowOffLeft)), 10, 'normalized');
        %plot(lags/(1/binsize_factor),spikesZvelocLeft)


    end
    meanTrialXvelocZvelocRight(k,:) = mean(tempRight,'omitnan');
    meanTrialXvelocZvelocLeft(k,:) = mean(tempLeft,'omitnan');
end

%%%% MIRROR RIGHT Z?%%%%%%
meanOverallXvelocZvelocRight = mean(meanTrialXvelocZvelocRight);
meanOverallXvelocZvelocLeft = mean(meanTrialXvelocZvelocLeft);
meanOverallXvelocZvelocBoth = [meanOverallXvelocZvelocRight + meanOverallXvelocZvelocLeft]/2;

figure
tiledlayout(3,1)

nexttile
plot(lags/(1/binsize_factor), meanTrialXvelocZvelocRight, 'Color', [0.7 0.7 0.7])
hold on
plot(lags/(1/binsize_factor), meanOverallXvelocZvelocRight, 'k', 'LineWidth', 3)
title('RIGHT Antenna', 'Cross correlation of X Velocity vs Z Velocity')
xlabel("Lag (s)")
ylabel("Correlation Coefficient")
ylim([-0.5 0.5])

nexttile
plot(lags/(1/binsize_factor), meanTrialXvelocZvelocLeft, 'Color', [0.7 0.7 0.7])
hold on
plot(lags/(1/binsize_factor), meanOverallXvelocZvelocLeft, 'k', 'LineWidth', 3)
title('LEFT Antenna', 'Cross correlation of X Velocity vs Z Velocity')
xlabel("Lag (s)")
ylabel("Correlation Coefficient")
ylim([-0.5 0.5])

nexttile
plot(lags/(1/binsize_factor), meanOverallXvelocZvelocRight)
hold on
plot(lags/(1/binsize_factor), meanOverallXvelocZvelocLeft)
hold on
plot(lags/(1/binsize_factor), meanOverallXvelocZvelocBoth)

xlabel('Lag (s)')
ylabel('Correlation Coefficient')
title('Cross correlation of X Velocity vs Z Velocity')
legend('Right', 'Left', 'Both','Location','northwest')

set(gcf,'position',[400, 100, 400,1000])

% print(['crosscorr_X_Z' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['crosscorr_X_Z' '.png'], '-dpng','-r300', '-vector')


%% Binning for data 05122023 excluding trials

for n = 1 : length(datasets)

    data = datasets{n};

    clearvars binstart binstarts binsize spikesbinned binned_spikerate_right

    binstarts = nan;
    binsize = binsize_factor*sampling_rate;
    binned_spikerate= cell(1,length(data));

    for k = 1 : length(data)
        temp_spikes = [];
        temp_xveloc = [];
        temp_zveloc_degree = [];

        for m = 1 : numPokes
            window_on = data(k).trigger_on(m)-1*sampling_rate; %1s before poking
            window_off = data(k).trigger_off(m)+1*sampling_rate; %1s afte poking

            binstarts = window_on:binsize:window_off;
            for bin = 1 : length(binstarts)-1
                temp_spikes(m,bin) = sum(data(k).spikes...
                    (binstarts(bin):binstarts(bin+1)),'omitnan');
                temp_xveloc(m,bin) = sum(data(k).xveloc_in_mm...
                    (binstarts(bin):binstarts(bin+1)),'omitnan');
                temp_zveloc_degree(m,bin) = sum(data(k).zveloc_in_degree_per_s...
                    (binstarts(bin):binstarts(bin+1)),'omitnan');
            end
            if to_be_removed_ID{n}(k,m) == 1
                temp_spikes(m,:) = NaN;
                temp_xveloc(m,:) = NaN;
                temp_zveloc_degree(m,:) = NaN;

            else
            end


        end
        temp_xveloc(any(isnan(temp_xveloc), 2), :) = [];
        temp_spikes(any(isnan(temp_spikes), 2), :) = [];
        temp_zveloc_degree(any(isnan(temp_zveloc_degree), 2), :) = [];
        % right_antenna(k).binnedSpikerate = cellfun(@(x)x/binsize_factor, temp_spikes);
        data(k).binnedSpikerate = temp_spikes/binsize_factor;
        data(k).binnedSpikerateMean = mean(data(k).binnedSpikerate);
        data(k).binnedXveloc = temp_xveloc/binsize;
        data(k).binnedXvelocMean = mean(data(k).binnedXveloc);
        data(k).binnedZvelocDegree = temp_zveloc_degree/binsize;
        data(k).binnedZvelocDegreeMean = mean(data(k).binnedZvelocDegree);

    end

    datasets{n} = data;
end

right_antenna_filtered = datasets{1};
left_antenna_filtered = datasets{2};


%% Plot example cell 

fly_number = 5;

data_R = right_antenna(fly_number);
data_L = left_antenna(fly_number);

figure
tiledlayout('flow')
nexttile
k = 0;
for i = 1 : numPokes    
    window_on_R = data_R.trigger_on(i)-sampling_rate;
    window_off_R = data_R.trigger_off(i)+sampling_rate;
    window_on_L = data_L.trigger_on(i)-sampling_rate;
    window_off_L = data_L.trigger_off(i)+sampling_rate;
    plot(data_R.VM(window_on_R:window_off_R)+k, 'k')
    hold on
    plot(data_L.VM(window_on_L:window_off_L)+(k+10), 'm')
    k = k + 10;
end
nexttile
plot(xbin, meanSpikerateCellRight{fly_number}, 'k', xbin, meanSpikerateCellLeft{fly_number}, 'm');
nexttile
plot(xbin, meanXvelocCellRight{fly_number}, 'k', xbin, meanXvelocCellLeft{fly_number}, 'm');
nexttile
plot(xbin, meanZvelocCellRight{fly_number}, 'k',xbin, meanZvelocCellLeft{fly_number}, 'm');
% nexttile
% plot(xbin, meanVM_medCellRight{fly_number}, xbin, meanVM_medCellLeft{fly_number});
legend('Right', 'Left')

% print(['example_cell' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['example_cell' '.png'], '-dpng','-r300', '-vector')


%% Compare trials left vs right side

%Change in spike frequency, x and z velocity before vs during poke
for i = 1 : numCells
deltaSpikesLeft{i} = spikesDuringPokeLeft{i} - spikesPrePokeLeft{i};
deltaSpikesRight{i}  = spikesDuringPokeRight{i} - spikesPrePokeRight{i};
deltaXVelocLeft{i}  = xVelocityDuringPokeLeft{i} - xVelocityPrePokeLeft{i};
deltaXVelocRight{i} = xVelocityDuringPokeRight{i} - xVelocityPrePokeRight{i};
deltaZVelocLeft{i}  = zVelocityDuringPokeLeft{i} - zVelocityPrePokeLeft{i};
deltaZVelocRight{i}  = zVelocityDuringPokeRight{i} - zVelocityPrePokeRight{i};

%statistics
pValSpikeRate(i) = signrank(deltaSpikesLeft{i}, deltaSpikesRight{i});
pValAngularVel(i) = signrank(deltaZVelocLeft{i}, deltaZVelocRight{i});
pValForwardVel(i) = signrank(deltaXVelocLeft{i}, deltaXVelocRight{i});
end

% is there tuning? No. of MDNs withsignificant differences in left right
% response to antennal poking in spike rate
pValSpikeRate(pValSpikeRate < 0.05) = 0;
pValSpikeRate(pValSpikeRate > 0.05) = 1;

pValAngularVel(pValAngularVel < 0.05) = 0;
pValAngularVel(pValAngularVel > 0.05) = 1;

pValForwardVel(pValForwardVel < 0.05) = 0;
pValForwardVel(pValForwardVel > 0.05) = 1;

isNotSignificantSpikeRate = sum(pValSpikeRate);
isNotSignificantAngularVel =sum(pValAngularVel);
isNotSignificantForwardVel =sum(pValForwardVel);

y1 = [isNotSignificantSpikeRate; numCells-isNotSignificantSpikeRate];
y2 = [isNotSignificantAngularVel; numCells-isNotSignificantAngularVel];
y3 = [isNotSignificantForwardVel; numCells-isNotSignificantForwardVel];
y = [isNotSignificantSpikeRate numCells-isNotSignificantSpikeRate;...
    isNotSignificantForwardVel numCells-isNotSignificantForwardVel;...
    isNotSignificantAngularVel numCells-isNotSignificantAngularVel];
x = categorical({'delta Spike Rate', 'delta Forward Velocity', 'delta Angular Velocity'});
x = reordercats(x,{'delta Spike Rate', 'delta Forward Velocity', 'delta Angular Velocity'});
figure
nexttile
bar(x, y, 'stacked')
xtickangle(45)
text(1,5, [num2str(isNotSignificantSpikeRate) '/' num2str(numCells)])
text(2,5, [num2str(isNotSignificantForwardVel) '/' num2str(numCells)])
text(3,5, [num2str(isNotSignificantAngularVel) '/' num2str(numCells)])
legend("Not Significant, p > 0.05", "Significant, p < 0.05")
title("Difference in tuning between Left and Right Antenna")


k = 0;
figure
tiledlayout(3,1);

for i = 1 : numCells
%Tile for spikes     
t1 = nexttile(1);
y = deltaSpikesLeft{i};
scatter(k,y,'o')
hold on 
y2 = deltaSpikesRight{i};
scatter(k+1,y2,'o')
for n = 1 : length(deltaSpikesRight{i})
    line([k k+1], [y(n) y2(n)])
end
yline(0)
ylabel('delta Spike Frequency (Hz)')
title('Change in Spike Frequency during Poking in Left vs Right Antenna')

%Tile for X Velocity     
t2 = nexttile(2);
y = deltaXVelocLeft{i};
scatter(k,y,'o')
hold on 
y2 = deltaXVelocRight{i};
scatter(k+1,y2,'o')
for n = 1 : length(deltaXVelocRight{i})
    line([k k+1], [y(n) y2(n)])
end
yline(0)
ylabel('delta Forward Velocity (mm/s)')
title('Change in Forward Velocity during Poking in Left vs Right Antenna')

%Tile for Z Velocity
t3 = nexttile(3);
y = deltaZVelocLeft{i};
scatter(k,y,'o')
hold on 
y2 = deltaZVelocRight{i};
scatter(k+1,y2,'o')
for n = 1 : length(deltaZVelocRight{i})
    line([k k+1], [y(n) y2(n)])
end
ylabel('delta Angular Velocity (°/s)')
k = k+2;
end
yline(0)
linkaxes([t1 t2 t3],'x')
xlim([-1 20])
title('Change in Angular Velocity during Poking in Left vs Right Antenna')

%% Total walking acivity 
clearvars totalDistancePerFly

datasets = {right_antenna, left_antenna};
datasets_names = {'right_antenna', 'left_antenna'};

for n = 1 : length(datasets)

    data = datasets{n};

    for i = 1 : length(data)
        data(i).motion = [data(i).VelocX, data(i).VelocY, data(i).VelocZ];
        currentMotion  = data(i).motion(1:data_reduction_for_plotting:end, :);
        allTrajectories = ball2trajectory(currentMotion);
        T = allTrajectories;
        x = T(1:end,2);
        y = T(1:end,3);
        z = zeros(size(x));
        d = hypot(diff(x), diff(y));                            % Distance Of Each Segment
        totalDistancePerFly(i,:) = sum(d);
    end

    datasets{n} = data;
end

right_antenna= datasets{1};
left_antenna = datasets{2};

%% Plot poking events in trajectory plots

clearvars trigger_on trigger_off

for k = 1% : length(right_antenna)

    for n = 1 : numPokes
        trigger_on(n,k) = right_antenna(k).trigger_on(n);
        trigger_off(n,k) = right_antenna(k).trigger_off(n)-20000;
    end

    currentMotion  = right_antenna(k).motion(1:data_reduction_for_plotting:end, :);

    allTrajectories = ball2trajectory(currentMotion);


    figure('Name','Trajecory');
        title(right_antenna(k).ID,'Interpreter','none') %% Interpreter 'none' required for preseting underscores (_) properly

    hold on
    xline(0,'Color',[0.8,0.8,0.8]);
    yline(0,'Color',[0.8,0.8,0.8]);
    caxis('manual')
    caxis([0 6])
    T = allTrajectories;
    x = T(1:end,2);
    y = T(1:end,3);
    z = zeros(size(x));

    scatter(x(1), y(1), [3000], '.')

    % scatter(x(1:end-2), y(1:end-2), [], 'm',  'filled', 'o');
    for m = 1 : length(trigger_on)
     
        xPlot = x(floor(trigger_on(m,k)/400):floor(trigger_off(m,k)/400));
        yPlot = y(floor(trigger_on(m,k)/400):floor(trigger_off(m,k)/400));
        xPlot = xPlot-xPlot(1);
        yPlot = yPlot-yPlot(1);
        
        c = right_antenna(k).spikesbinned_50Hz(floor(trigger_on(m,k)/400):floor(trigger_off(m,k)/400));

    s = scatter(xPlot, yPlot, [], 'c', 'filled', 'o');
    line(xPlot,yPlot, 'LineWidth', 3)       
    s.Marker = 'none';
    hold on
    end

    colormap(flipud(cool))
    colorbar;
    % caxis([0 1])


    plot(0,0,...
        'Marker','o', 'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    %Total distance walked
    d = hypot(diff(x), diff(y));                            % Distance Of Each Segment
    d_tot = sum(d);

    xlabel('X dimension [mm]');
    ylabel('Y dimension [mm]');
    axis equal;

end

%% %% Binning for trajcetory firing rate plot

clearvars binstart binstarts binsize

binsize_factor_50Hz = 0.02;
%binsize_factor_50Hz = 0.1;      % 100 ms bins

binstarts = nan;
binsize = binsize_factor_50Hz*sampling_rate;

for k = 1 : length(right_antenna)
    
    binstarts = 1:binsize:length(right_antenna(k).spikes);
    spikesbinned_50Hz = nan;
    for bin = 1 : length(binstarts)-1
        spikesbinned_50Hz(:,bin) = sum(right_antenna(k).spikes(binstarts(bin):binstarts(bin+1)),'omitnan');
    end
    spikesbinned_50Hz = spikesbinned_50Hz/binsize_factor;
    spikesnorm_binned = spikesbinned_50Hz/max(spikesbinned_50Hz);
    right_antenna(k).spikesbinned_50Hz = spikesbinned_50Hz;
    right_antenna(k).spikesnorm_binned = spikesnorm_binned';
       
end
