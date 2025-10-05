%% Load data

load('\prePoker_behavior.mat')

filepath='\behavior_only';

smoothing_factor = 2;
sampling_rate = 20000;
data_reduction_for_plotting = 400;
binsize_factor = 0.2; %200 ms
numPokes = 10;
%t_ramp = 0.5; % 0.5s ramp time
t_poke = 2; % 2s poke hold
interstim_interval = 10; %10s between pokes

binsize = 4000;
bin = 20;

%% trajectory

for k = 1 : length(analysis)
    analysis(k).motion = [analysis(k).VelocX, analysis(k).VelocY, analysis(k).VelocZ];
end
figure('Name', 'Trajectory')
%tiledlayout(round(length(analysis)/4), 4)
%tiledlayout(1,2)

for k = 1 : length(analysis)

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
    yleft = T(1:end,3);
    z = zeros(size(x));

    scatter(x(1), yleft(1), [3000], '.')

    %c = analysis(i).spikesnorm_binned;

    scatter(x(1:end-2), yleft(1:end-2), [], 'c',  'filled', 'o');


    %scatter(x(1:end-2), y(1:end-2), [], 'filled', 'o');

    colormap(flipud(cool))
    colorbar;
    caxis([0 1])


    plot(0,0,...
        'Marker','o', 'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    % Total distance walked
    d = hypot(diff(x), diff(yleft));                            % Distance Of Each Segment
    d_tot = sum(d);

    xlabel('X dimension [mm]');
    ylabel('Y dimension [mm]');
    axis equal;

    analysis(k).mean_velocity = d_tot/((length(currentMotion)*400)/sampling_rate);
end

mean_velocity = d_tot/((length(currentMotion)*400)/sampling_rate); %%mean velocity in mm/s



%% X velocity binning during stimulation window
for k = 1 : length(analysis)
    analysis(k).xveloc_in_mm = analysis(k).VelocX(:,1)*8.79;
end
%%%%% Right Antenna
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

%% Z velocity preliminary

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

    for m = 1 : numPokes
        window_on = analysis(k).trigger_on(m)-1*sampling_rate; %1s before poking
        window_off = analysis(k).trigger_off(m)+1*sampling_rate; %1s after poking

        binstarts = window_on:binsize:window_off;
        for bin = 1 : length(binstarts)-1
            temp_zveloc{m,bin} = sum(analysis(k).zveloc_in_mm(binstarts(bin):binstarts(bin+1)),'omitnan');
        end
    end
    analysis(k).binnedZveloc = cellfun(@(x)x/binsize, temp_zveloc);
    analysis(k).binnedZvelocMean = mean(analysis(k).binnedZveloc);
end

%% Analyze poker effect
clearvars right_antenna_control left_antenna_control right_antenna_kir...
    left_antenna_kir

for i = 1 : length(analysis)
    if strcmp(analysis(i).genotype, 'control')

        if strcmp(analysis(i).antenna, 'right')
            right_antenna_control(i)= analysis(i);

        else
            left_antenna_control(i) = analysis(i);

        end
    else
        if strcmp(analysis(i).antenna, 'right')
            right_antenna_kir(i)= analysis(i);

        else
            left_antenna_kir(i) = analysis(i);

        end
    end
end
%
right_antenna_control = right_antenna_control(~cellfun(@isempty,{right_antenna_control.VelocX}));
left_antenna_control = left_antenna_control(~cellfun(@isempty,{left_antenna_control.VelocX}));

right_antenna_kir = right_antenna_kir(~cellfun(@isempty,{right_antenna_kir.VelocX}));
left_antenna_kir = left_antenna_kir(~cellfun(@isempty,{left_antenna_kir.VelocX}));

% 
save('postPoker_right_control_all','right_antenna_control','-v7.3');
save('postPoker_left_control_all','left_antenna_control','-v7.3');

save('postPoker_right_kir_all','right_antenna_kir','-v7.3');
save('postPoker_left_kir_all','left_antenna_kir','-v7.3');

save('postAnalysis', 'analysis', '-v7.3');


%% Quality control
% clearvars right_antenna_control_min_velocity left_antenna_control_min_velocity
%
% right_antenna_control_min_velocity = right_antenna_control(~cellfun(@(x)x<0.5,{right_antenna_control.mean_velocity}));
% left_antenna_control_min_velocity  = left_antenna_control(~cellfun(@(x)x<0.5,{left_antenna_control.mean_velocity}));


%% Plotting mean and binned spikerates during stimulation with 1 s before and 1 s after ---- CONTROL
clearvars binned_spikerate_right_mean binned_spikerate_left_mean binned_xveloc_right_mean...
    binned_xveloc_left_mean binned_zveloc_right_mean binned_zveloc_left_mean xbin


meanXvelocCellControlRight = arrayfun(@(s) s.binnedXvelocMean, right_antenna_control, 'UniformOutput', false);
meanXvelocOverallControlRight = mean(cat(1, meanXvelocCellControlRight{:}));

meanXvelocCellControlLeft = arrayfun(@(s) s.binnedXvelocMean, left_antenna_control, 'UniformOutput', false);
meanXvelocOverallControlLeft = mean(cat(1, meanXvelocCellControlLeft{:}));

meanZvelocCellControlRight = arrayfun(@(s) s.binnedZvelocMean, right_antenna_control, 'UniformOutput', false);
meanZvelocOverallControlRight = mean(cat(1, meanZvelocCellControlRight{:}));

meanZvelocCellControlLeft = arrayfun(@(s) s.binnedZvelocMean, left_antenna_control, 'UniformOutput', false);
meanZvelocOverallControlLeft = mean(cat(1, meanZvelocCellControlLeft{:}));

%xbin = [1:bin]/(sampling_rate/binsize);
xbin = ((binsize/2)/sampling_rate):binsize/sampling_rate:(bin*(binsize/sampling_rate));

% Create a colormap
cmap = hsv(numel(right_antenna_control));

%Plot right antenna firing rate and x velocity
figure('name', 'antenna poking 2s poke 0.2 bins')
t = tiledlayout(2,2);
title(t,'CONTROL')

nexttile
for i = 1 : length(right_antenna_control)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanXvelocCellControlLeft(i)), 'Color', color)
    hold on
end
plot(xbin, meanXvelocOverallControlLeft, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('X velocity LEFT')
ylim([-2 1])
ylabel('X Velocity (mm)');

nexttile
for i = 1 : length(right_antenna_control)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanXvelocCellControlRight(i)), 'Color', color)
    hold on
end
plot(xbin, meanXvelocOverallControlRight, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('X velocity RIGHT')
ylim([-2 1])
xlabel('X Velocity');
nexttile

for i = 1 : length(right_antenna_control)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanZvelocCellControlLeft(i)), 'Color', color)
    hold on
end
plot(xbin, meanZvelocOverallControlLeft, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('Z velocity LEFT')
ylim([-2.5 2.5])
ylabel('Z Velocity (mm)');
xlabel('Time (s)')

nexttile
for i = 1 : length(right_antenna_control)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanZvelocCellControlRight(i)), 'Color', color)
    hold on
end
plot(xbin, meanZvelocOverallControlRight, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('Z velocity RIGHT')
ylim([-2.5 2.5])
xlabel('Time (s)')

% Save the figure as a high-quality PNG image
%print(gcf,'poke1.png','-dpng','-r300');

figure('name', 'Overview')
t = tiledlayout('flow');
title(t,'CONTROL');
nexttile

% calculate error of x velocity for right and left antenna
stdev_left = std(cat(1, meanXvelocCellControlLeft{:}));
stderror_left = stdev_left/sqrt(length(cat(1, meanXvelocCellControlLeft{:})));
stdev_right = std(cat(1, meanXvelocCellControlRight{:}));
stderror_right = stdev_right/sqrt(length(cat(1, meanXvelocCellControlRight{:})));

% %Add SEM to y upper and limit and subtract from lower limit
% y_upper_limit_X_Control_Left = meanXvelocOverallControlLeft + stderror_left;
% y_lower_limit_X_Control_Left = meanXvelocOverallControlLeft - stderror_left;
% y_upper_limit_X_Control_Right = meanXvelocOverallControlRight + stderror_right;
% y_lower_limit_X_Control_Right = meanXvelocOverallControlRight - stderror_right;

%Add SD to y upper and limit and subtract from lower limit
y_upper_limit_X_Control_Left = meanXvelocOverallControlLeft + stdev_left;
y_lower_limit_X_Control_Left = meanXvelocOverallControlLeft - stdev_left;
y_upper_limit_X_Control_Right = meanXvelocOverallControlRight + stdev_right;
y_lower_limit_X_Control_Right = meanXvelocOverallControlRight - stdev_right;

%plot Left antenna x velocity with error as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_X_Control_Left fliplr(y_lower_limit_X_Control_Left)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallControlLeft, 'LineWidth',3)

%plot Left antenna x velocity with error as shaded area
hold on
fill([xbin fliplr(xbin)], [y_upper_limit_X_Control_Right fliplr(y_lower_limit_X_Control_Right)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanXvelocOverallControlRight, 'LineWidth',3)

title('Control with SEM')
legend( '','left antenna','','right antenna','Location','northwest')
title('X velocity')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'w')
ylim([-1 0.5])
ylabel('X Velocity (mm))');
legend( '','left antenna','','right antenna','Location','northwest')

nexttile

% calculate error of Z velocity for right and left antenna

stdev_left = std(cat(1, meanZvelocCellControlLeft{:}));
stderror_left = stdev_left/sqrt(length(cat(1, meanZvelocCellControlLeft{:})));
stdev_right = std(cat(1, meanZvelocCellControlRight{:}));
stderror_right = stdev_right/sqrt(length(cat(1, meanZvelocCellControlRight{:})));

%
% %Add SEM to y upper and limit and subtract from lower limit
% y_upper_limit_Z_Control_Left = meanZvelocOverallControlLeft + stderror_left;
% y_lower_limit_Z_Control_Left = meanZvelocOverallControlLeft - stderror_left;
% y_upper_limit_Z_Control_Right = meanZvelocOverallControlRight + stderror_right;
% y_lower_limit_Z_Control_Right = meanZvelocOverallControlRight - stderror_right;

%Add error to y upper and limit and subtract from lower limit
y_upper_limit_Z_Control_Left = meanZvelocOverallControlLeft + stdev_left;
y_lower_limit_Z_Control_Left = meanZvelocOverallControlLeft - stdev_left;
y_upper_limit_Z_Control_Right = meanZvelocOverallControlRight + stdev_right;
y_lower_limit_Z_Control_Right = meanZvelocOverallControlRight - stdev_right;
%plot Left antenna x velocity with error as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Control_Left fliplr(y_lower_limit_Z_Control_Left)],...
    [1 0 0], 'FaceAlpha', 0.3,'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallControlLeft, 'LineWidth',3)

%plot Left antenna x velocity with error as shaded area
hold on
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Control_Right fliplr(y_lower_limit_Z_Control_Right)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanZvelocOverallControlRight, 'LineWidth',3)
title('Z velocity')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylim([-1 2])
ylabel('Z Velocity (mm))');
xlabel('Time (s)');

% Save the figure as a high-quality PNG image
%print(gcf,'poke2.png','-dpng','-r300');

%Plot all cells with firing rate, x and z velocities next to each other
figure('name', 'Cell by cell comparison -- 2s poke 0.2 bins')

t = tiledlayout(2,length(right_antenna_control),"TileSpacing","compact");
title(t,'All cells --- CONTROL')

for i = 1 : length(right_antenna_control)
    nexttile;
    plot(xbin, meanXvelocCellControlRight{i})
    hold on
    plot(xbin, meanXvelocCellControlLeft{i})
    set(gca, 'XTickLabel', [])
    ylim([-1.5 0.5])
end

for i = 1 : length(right_antenna_control)
    nexttile;
    plot(xbin, meanZvelocCellControlRight{i})
    hold on
    plot(xbin, meanZvelocCellControlLeft{i})
    ylim([-3 3])
end

xlabel(t,'Time (s)')

% print(['control_overview' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['control_overview' '.png'], '-dpng','-r300', '-vector')

%% Plotting mean and binned spikerates during stimulation with 1 s before and 1 s after ---- SILENCED
clearvars binned_spikerate_right_mean binned_spikerate_left_mean binned_xveloc_right_mean...
    binned_xveloc_left_mean binned_zveloc_right_mean binned_zveloc_left_mean xbin


meanXvelocCellKirRight = arrayfun(@(s) s.binnedXvelocMean, right_antenna_kir, 'UniformOutput', false);
meanXvelocOverallKirRight = mean(cat(1, meanXvelocCellKirRight{:}));

meanXvelocCellKirLeft = arrayfun(@(s) s.binnedXvelocMean, left_antenna_kir, 'UniformOutput', false);
meanXvelocOverallKirLeft = mean(cat(1, meanXvelocCellKirLeft{:}));

meanZvelocCellKirRight = arrayfun(@(s) s.binnedZvelocMean, right_antenna_kir, 'UniformOutput', false);
meanZvelocOverallKirRight = mean(cat(1, meanZvelocCellKirRight{:}));

meanZvelocCellKirLeft = arrayfun(@(s) s.binnedZvelocMean, left_antenna_kir, 'UniformOutput', false);
meanZvelocOverallKirLeft = mean(cat(1, meanZvelocCellKirLeft{:}));

%xbin = [1:bin]/(sampling_rate/binsize);
xbin = ((binsize/2)/sampling_rate):binsize/sampling_rate:(bin*(binsize/sampling_rate));

% Create a colormap
cmap = hsv(numel(right_antenna_kir));

%Plot right antenna firing rate and x velocity
figure('name', 'Right antenna poking 2s poke 0.2 bins')
t = tiledlayout(2,2)
title(t,'SILENCED')

nexttile
for i = 1 : length(left_antenna_kir)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanXvelocCellKirLeft(i)), 'Color', color)
    hold on
end
plot(xbin, meanXvelocOverallKirLeft, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('X velocity LEFT')
ylim([-2 1])
xlabel('X Velocity');

nexttile
for i = 1 : length(right_antenna_kir)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanZvelocCellKirRight(i)), 'Color', color)
    hold on
end
plot(xbin, meanZvelocOverallKirRight, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('Z velocity RIGHT')
ylim([-2.5 2.5])
ylabel('Z Velocity (mm)');
xlabel('Time (s)')

nexttile
for i = 1 : length(right_antenna_kir)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanZvelocCellKirLeft(i)), 'Color', color)
    hold on
end
plot(xbin, meanZvelocOverallKirLeft, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('Z velocity LEFT')
ylim([-2.5 2.5])
xlabel('Time (s)')
nexttile

for i = 1 : length(right_antenna_kir)
    color = cmap(i,:);
    plot(xbin, cell2mat(meanXvelocCellKirRight(i)), 'Color', color)
    hold on
end
plot(xbin, meanXvelocOverallKirRight, 'k','Linewidth', 3)
%rectangle('Position',[2 -1 0.5 3], 'FaceColor', [1 0 1 0.15], 'EdgeColor', 'w')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0 0 1 0.15], 'EdgeColor', 'w')
title('X velocity RIGHT')
ylim([-2 1])
ylabel('X Velocity (mm)');
% Save the figure as a high-quality PNG image
%print(gcf,'poke1.png','-dpng','-r300');
figure('name', 'Overview')
t = tiledlayout('flow');
title(t,'SILENCED');
nexttile

% calculate error of x velocity for right and left antenna
stdev_left = std(cat(1, meanXvelocCellKirLeft{:}));
stderror_left = stdev_left/sqrt(length(cat(1, meanXvelocCellKirLeft{:})));
stdev_right = std(cat(1, meanXvelocCellKirRight{:}));
stderror_right = stdev_right/sqrt(length(cat(1, meanXvelocCellKirRight{:})));

% %Add error to y upper and limit and subtract from lower limit
% y_upper_limit_X_Kir_Left = meanXvelocOverallKirLeft + stderror_left;
% y_lower_limit_X_Kir_Left = meanXvelocOverallKirLeft - stderror_left;
% y_upper_limit_X_Kir_Right = meanXvelocOverallKirRight + stderror_right;
% y_lower_limit_X_Kir_Right = meanXvelocOverallKirRight - stderror_right;

%Add error to y upper and limit and subtract from lower limit
y_upper_limit_X_Kir_Left = meanXvelocOverallKirLeft + stdev_left;
y_lower_limit_X_Kir_Left = meanXvelocOverallKirLeft - stdev_left;
y_upper_limit_X_Kir_Right = meanXvelocOverallKirRight + stdev_right;
y_lower_limit_X_Kir_Right = meanXvelocOverallKirRight - stdev_right;
%plot Left antenna x velocity with error as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_X_Kir_Left fliplr(y_lower_limit_X_Kir_Left)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallKirLeft, 'LineWidth',3)

%plot Left antenna x velocity with error as shaded area
hold on
fill([xbin fliplr(xbin)], [y_upper_limit_X_Kir_Right fliplr(y_lower_limit_X_Kir_Right)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanXvelocOverallKirRight, 'LineWidth',3)

legend( '','left antenna','','right antenna','Location','northwest')
title('X velocity')
rectangle('Position',[1 -2 2 3], 'FaceColor', [0.7 0.7 0.7 0.2], 'EdgeColor', 'w')
ylim([-1 0.5])
ylabel('X Velocity (mm))');
legend( '','left antenna','','right antenna','Location','northwest')

clearvars stdev_left stdev_right stderror_left stderror_right
% calculate error of Z velocity for right and left antenna
stdev_left = std(cat(1, meanZvelocCellKirLeft{:}));
stderror_left = stdev_left/sqrt(length(cat(1, meanZvelocCellKirLeft{:})));
stdev_right = std(cat(1, meanZvelocCellKirRight{:}));
stderror_right = stdev_right/sqrt(length(cat(1, meanZvelocCellKirRight{:})));

% %Add error to y upper and limit and subtract from lower limit
% y_upper_limit_Z_Kir_Left = meanZvelocOverallKirLeft + stderror_left;
% y_lower_limit_Z_Kir_Left = meanZvelocOverallKirLeft - stderror_left;
% y_upper_limit_Z_Kir_Right = meanZvelocOverallKirRight + stderror_right;
% y_lower_limit_Z_Kir_Right = meanZvelocOverallKirRight - stderror_right;

%Add error to y upper and limit and subtract from lower limit
y_upper_limit_Z_Kir_Left = meanZvelocOverallKirLeft + stdev_left;
y_lower_limit_Z_Kir_Left = meanZvelocOverallKirLeft - stdev_left;
y_upper_limit_Z_Kir_Right = meanZvelocOverallKirRight + stdev_right;
y_lower_limit_Z_Kir_Right = meanZvelocOverallKirRight - stdev_right;

nexttile
%plot Left antenna x velocity with error as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Kir_Left fliplr(y_lower_limit_Z_Kir_Left)],...
    [1 0 0], 'FaceAlpha', 0.3,'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallKirLeft, 'LineWidth',3)

%plot Left antenna x velocity with error as shaded area
hold on
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Kir_Right fliplr(y_lower_limit_Z_Kir_Right)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanZvelocOverallKirRight, 'LineWidth',3)
title('Z velocity')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylim([-1 2])
ylabel('Z Velocity (mm))');
xlabel('Time (s)');

% Save the figure as a high-quality PNG image
%print(gcf,'poke2.png','-dpng','-r300');

%Plot all cells with firing rate, x and z velocities next to each other
figure('name', 'Cell by cell comparison -- 2s poke 0.2 bins')

t = tiledlayout(2,length(right_antenna_kir),"TileSpacing","compact");
title(t,'All cells --- SILENCED')


for i = 1 : length(right_antenna_kir)
    nexttile;
    plot(xbin, meanXvelocCellKirRight{i})
    hold on
    plot(xbin, meanXvelocCellKirLeft{i})
    set(gca, 'XTickLabel', [])
    ylim([-1.5 0.5])
end

for i = 1 : length(right_antenna_kir)
    nexttile;
    plot(xbin, meanZvelocCellKirRight{i})
    hold on
    plot(xbin, meanZvelocCellKirLeft{i})
    ylim([-3 3])
end

xlabel(t,'Time (s)')

% print(['silenced_overview' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['silenced_overview' '.png'], '-dpng','-r300', '-vector')
%% Compare CONTROL vs SILENCED

figure('name', 'CONTROL vs SILENCED')

tiledlayout(2,2)

nexttile;
% fill([xbin fliplr(xbin)], [y_upper_limit_X_Control_Left fliplr(y_lower_limit_X_Control_Left)],...
% [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallControlLeft, 'LineWidth',3)
% fill([xbin fliplr(xbin)], [y_upper_limit_X_Kir_Left fliplr(y_lower_limit_X_Kir_Left)],...
% [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallKirLeft, 'LineWidth',3)

title('LEFT')
ylabel('X Velocity (mm))');
yline(0,'--');
% legend('', 'control +- SD', '', 'silenced +- SD','Location','southeast','FontSize',12)
legend('control', 'silenced','Location','southeast','FontSize',12)
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylim([-1 0.5])

nexttile;
% fill([xbin fliplr(xbin)], [y_upper_limit_X_Control_Right fliplr(y_lower_limit_X_Control_Right)],...
% [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallControlRight, 'LineWidth',3)
% fill([xbin fliplr(xbin)], [y_upper_limit_X_Kir_Right fliplr(y_lower_limit_X_Kir_Right)],...
% [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallKirRight, 'LineWidth',3)
title('RIGHT')
ylabel('X Velocity (mm))');
yline(0,'--');
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylim([-1 0.5])

nexttile
% fill([xbin fliplr(xbin)], [y_upper_limit_Z_Control_Left fliplr(y_lower_limit_Z_Control_Left)],...
% [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallControlLeft, 'LineWidth',3)
% fill([xbin fliplr(xbin)], [y_upper_limit_Z_Kir_Left fliplr(y_lower_limit_Z_Kir_Left)],...
% [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallKirLeft, 'LineWidth',3)

ylabel('Z Velocity (mm))');
yline(0,'--');
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylim([-1 1])

nexttile
% fill([xbin fliplr(xbin)], [y_upper_limit_Z_Control_Right fliplr(y_lower_limit_Z_Control_Right)],...
% [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallControlRight, 'LineWidth',3)
% fill([xbin fliplr(xbin)], [y_upper_limit_Z_Kir_Right fliplr(y_lower_limit_Z_Kir_Right)],...
% [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallKirRight, 'LineWidth',3)
ylabel('Z Velocity (mm))');
yline(0,'--');
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylim([-1 1])

% print(['control_vs_silence_overview' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['control_vs_silence_overview' '.png'], '-dpng','-r300', '-vector')

%% Analyse 1 s window after stim and compare X and Z velocity ---- CONTROL

%analysis window size in seconds
window_size = 0.25;

for i = 1 : length(left_antenna_control)

    for k = 1 : numPokes

        %Window for analyzing firing rate, x and z velocity, 1 sec
        window_on = left_antenna_control(i).trigger_on(k);
        window_off = left_antenna_control(i).trigger_on(k)+window_size*sampling_rate;

        %Add all data to an array
        XvelocLeft{i,k} = left_antenna_control(i).xveloc_in_mm(window_on:window_off);
        ZvelocLeft{i,k} = left_antenna_control(i).zveloc_in_mm(window_on:window_off);

    end

end

for i = 1 : length(right_antenna_control)

    for k = 1 : numPokes

        %Window for analyzing firing rate, x and z velocity, 1 sec
        window_on = right_antenna_control(i).trigger_on(k);
        window_off = right_antenna_control(i).trigger_on(k)+window_size*sampling_rate;

        %Add all data to an array
        XvelocRight{i,k} = right_antenna_control(i).xveloc_in_mm(window_on:window_off);
        ZvelocRight{i,k} = right_antenna_control(i).zveloc_in_mm(window_on:window_off);

    end

end

%Either take the mean or the minimum (maximal backward) x velocity
MeanXvelocLeft = cellfun(@(x) mean(x, 'all', 'omitnan'), XvelocLeft);
MedianXvelocLeft = cellfun(@(x) median(x, 'all', 'omitnan'), XvelocLeft);
MinXvelocLeft = cellfun(@(x) min(x), XvelocLeft);
%Either take the mean or the max(left antenna) z velocity
MeanZvelocLeft = cellfun(@(x) mean(x, 'all', 'omitnan'), ZvelocLeft);
MedianZvelocLeft = cellfun(@(x) median(x, 'all', 'omitnan'), ZvelocLeft);
MaxZvelocLeft = cellfun(@(x) max(x), XvelocLeft);

%Either take the mean or the minimum (maximal backward) x velocity
MeanXvelocRight= cellfun(@(x) mean(x, 'all', 'omitnan'), XvelocRight);
MedianXvelocRight = cellfun(@(x) median(x, 'all', 'omitnan'), XvelocRight);
MinXvelocRight = cellfun(@(x) min(x), XvelocRight);
%Either take the mean or the min(right antenna) z velocity
MeanZvelocRight = cellfun(@(x) mean(x, 'all', 'omitnan'), ZvelocRight);
MedianZvelocRight = cellfun(@(x) median(x, 'all', 'omitnan'), ZvelocRight);
MinZvelocRight = cellfun(@(x) min(x), XvelocRight);


figure
title(['CONTROL \newline MEAN x velocity for all flies ', num2str(window_size),  's window'])

groups = 1:length(right_antenna_control);


for i = 1 : length(right_antenna_control)

    %Plot Min X velocity of right and left antenna in one group per fly number
    plot_y = [MeanXvelocLeft(i,:)', MeanXvelocRight(i,:)'];

    % %save median of right and left x velocity per fly
    % med(i,:) = median(plot_y);
    % med_absolute = abs(med);
    % difference_left_right = abs(med_absolute(:,1) - med_absolute(:,2));
    %
    % %if difference in median between righ and left is larger than X dont plot
    %  if difference_left_right(i) > 0.5
    %      continue
    %  else
    %  end

    hbox{k} = boxplot2(plot_y, groups(i),'whisker', 0);

    %Plot individual values for each poking trial on top of the boxplots
    hold on
    scatter(groups(i)-0.15, plot_y(:,1), 250, '.', 'k')
    hold on
    scatter(groups(i)+0.15, plot_y(:,2), 250, '.', 'm')
    hold on
    plot(groups(i)-0.15, mean(plot_y(:,1)), 'kx', 'MarkerSize', 10)
    hold on
    plot(groups(i)+0.15, mean(plot_y(:,2)), 'mx', 'MarkerSize', 10)


    %set boxplots parameters,  no outliers, box colors etc
    set(hbox{k}.box(1, 1), 'Linewidth', 1, 'Color','k')
    set(hbox{k}.box(1, 2), 'Linewidth', 1, 'Color','m')
    set(hbox{k}.out, 'Marker', 'none')
    set(hbox{k}.med, 'Linewidth', 2 )
    set(hbox{k}.uadj, 'LineStyle', 'none')
    set(hbox{k}.ladj, 'LineStyle', 'none')
end

legend( 'left antenna', 'right antenna', 'Location', 'northwest');
xticks(1:length(right_antenna_control));
xlabel('#Fly')
ylabel('X Velocity (mm/s)')


figure
title(['CONTROL \newline absolute max z velocity for all flies ' , num2str(window_size),  's window'])

groups = 1:length(right_antenna_control);


for i = 1 : length(right_antenna_control)

    %Plot Min X velocity of right and left antenna in one group per fly number
    plot_y = [MaxZvelocLeft(i,:)', abs(MinZvelocRight(i,:)')];

    hbox{k} = boxplot2(plot_y, groups(i),'whisker', 0)

    %Plot individual values for each poking trial on top of the boxplots
    hold on
    scatter(groups(i)-0.15, plot_y(:,1), 250, '.', 'k')
    hold on
    scatter(groups(i)+0.15, plot_y(:,2), 250, '.', 'm')

    %set boxplots parameters,  no outliers, box colors etc
    set(hbox{k}.box(1, 1), 'Linewidth', 1, 'Color','k')
    set(hbox{k}.box(1, 2), 'Linewidth', 1, 'Color','m')
    set(hbox{k}.out, 'Marker', 'none')
    set(hbox{k}.med, 'Linewidth', 2 )
    set(hbox{k}.uadj, 'LineStyle', 'none')
    set(hbox{k}.ladj, 'LineStyle', 'none')

end

xticks(1:length(right_antenna_control));
xlabel('#Fly')
ylabel('absolute Z Velocity (mm/s)')

%% Analyse 1 s window after stim and compare X and Z velocity ---- SILENCED

%analysis window size in seconds
window_size = 1;

for i = 1 : length(left_antenna_kir)

    for k = 1 : numPokes

        %Window for analyzing firing rate, x and z velocity, 1 sec
        window_on = left_antenna_kir(i).trigger_on(k);
        window_off = left_antenna_kir(i).trigger_on(k)+window_size*sampling_rate;

        %Add all data to an array
        XvelocSilencedLeft{i,k} = left_antenna_kir(i).xveloc_in_mm(window_on:window_off);
        ZvelocSilencedLeft{i,k} = left_antenna_kir(i).zveloc_in_mm(window_on:window_off);

    end

end

for i = 1 : length(right_antenna_kir)

    for k = 1 : numPokes

        %Window for analyzing firing rate, x and z velocity, 1 sec
        window_on = right_antenna_kir(i).trigger_on(k);
        window_off = right_antenna_kir(i).trigger_on(k)+window_size*sampling_rate;

        %Add all data to an array
        XvelocSilencedRight{i,k} = right_antenna_kir(i).xveloc_in_mm(window_on:window_off);
        ZvelocSilencedRight{i,k} = right_antenna_kir(i).zveloc_in_mm(window_on:window_off);

    end

end

%Either take the mean or the minimum (maximal backward) x velocity
MeanXvelocSilencedLeft = cellfun(@(x) mean(x, 'all', 'omitnan'), XvelocSilencedLeft);
MedianXvelocSilencedLeft = cellfun(@(x) median(x, 'all', 'omitnan'), XvelocSilencedLeft);
MinXvelocSilencedLeft = cellfun(@(x) min(x), XvelocSilencedLeft);
%Either take the mean or the max(left antenna) z velocity
MeanZvelocSilencedLeft = cellfun(@(x) mean(x, 'all', 'omitnan'), ZvelocSilencedLeft);
MedianZvelocSilencedLeft = cellfun(@(x) median(x, 'all', 'omitnan'), ZvelocSilencedLeft);
MaxZvelocSilencedLeft = cellfun(@(x) max(x), ZvelocSilencedLeft);

%Either take the mean or the minimum (maximal backward) x velocity
MeanXvelocSilencedRight= cellfun(@(x) mean(x, 'all', 'omitnan'), XvelocSilencedRight);
MedianXvelocSilencedRight = cellfun(@(x) median(x, 'all', 'omitnan'), XvelocSilencedRight);
MinXvelocSilencedRight = cellfun(@(x) min(x), XvelocSilencedRight);
%Either take the mean or the min(right antenna) z velocity
MeanZvelocSilencedRight = cellfun(@(x) mean(x, 'all', 'omitnan'), ZvelocSilencedRight);
MedianZvelocSilencedRight = cellfun(@(x) median(x, 'all', 'omitnan'), ZvelocSilencedRight);
MinZvelocSilencedRight = cellfun(@(x) min(x), ZvelocSilencedRight);


figure
title(['SILENCED \newline min x velocity for all flies ', num2str(window_size),  's window'])

groups = 1:length(right_antenna_kir);

for i = 1 : length(right_antenna_kir)

    %Plot Min X velocity of right and left antenna in one group per fly number
    plot_y = [MinXvelocSilencedLeft(i,:)', MinXvelocSilencedRight(i,:)'];

    % %save median of right and left x velocity per fly
    % med(i,:) = median(plot_y);
    % med_absolute = abs(med);
    % difference_left_right = abs(med_absolute(:,1) - med_absolute(:,2));
    %
    % %if difference in median between righ and left is larger than X dont plot
    %  if difference_left_right(i) > 0.5
    %      continue
    %  else
    %  end

    hbox{k} = boxplot2(plot_y, groups(i),'whisker', 0);

    %Plot individual values for each poking trial on top of the boxplots
    hold on
    scatter(groups(i)-0.15, plot_y(:,1), 250, '.', 'k')
    hold on
    scatter(groups(i)+0.15, plot_y(:,2), 250, '.', 'm')

    %set boxplots parameters,  no outliers, box colors etc
    set(hbox{k}.box(1, 1), 'Linewidth', 1, 'Color','k')
    set(hbox{k}.box(1, 2), 'Linewidth', 1, 'Color','m')
    set(hbox{k}.out, 'Marker', 'none')
    set(hbox{k}.med, 'Linewidth', 2 )
    set(hbox{k}.uadj, 'LineStyle', 'none')
    set(hbox{k}.ladj, 'LineStyle', 'none')
end

legend( 'left antenna', 'right antenna', 'Location', 'northwest');
xticks(1:length(right_antenna_kir));
xlabel('#Fly')
ylabel('X Velocity (mm/s)')


figure
title(['SILENCED \newline absolute max z velocity for all flies ' , num2str(window_size),  's window'])

groups = 1:length(right_antenna_kir);


for i = 1 : length(right_antenna_kir)

    %Plot Min X velocity of right and left antenna in one group per fly number
    plot_y = [MaxZvelocSilencedLeft(i,:)', abs(MinZvelocSilencedRight(i,:)')];

    hbox{k} = boxplot2(plot_y, groups(i),'whisker', 0)

    %Plot individual values for each poking trial on top of the boxplots
    hold on
    scatter(groups(i)-0.15, plot_y(:,1), 250, '.', 'k')
    hold on
    scatter(groups(i)+0.15, plot_y(:,2), 250, '.', 'm')

    %set boxplots parameters,  no outliers, box colors etc
    set(hbox{k}.box(1, 1), 'Linewidth', 1, 'Color','k')
    set(hbox{k}.box(1, 2), 'Linewidth', 1, 'Color','m')
    set(hbox{k}.out, 'Marker', 'none')
    set(hbox{k}.med, 'Linewidth', 2 )
    set(hbox{k}.uadj, 'LineStyle', 'none')
    set(hbox{k}.ladj, 'LineStyle', 'none')

end

xticks(1:length(right_antenna_kir));
xlabel('#Fly')
ylabel('absolute Z Velocity (mm/s)')


%% Calculate mean X velocitx and create boxplots for CONTROL

%Size of poking analysis window in s
analysis_window_size = 0.5;


for i = 1 : length(right_antenna_control)
    for k = 1 : numPokes
        windowOnRight = right_antenna_control(i).trigger_on(k);
        windowOffRight = right_antenna_control(i).trigger_off(k);
        windowOnLeft = left_antenna_control(i).trigger_on(k);
        windowOffLeft = left_antenna_control(i).trigger_off(k);
        %Sum spikes 1 s pre poke
        xVelocityPrePokeControlRight{i}(k) = mean(right_antenna_control(i).xveloc_in_mm...
            (windowOnRight-20000:windowOnRight),'omitnan');
        xVelocityPrePokeControlLeft{i}(k) = mean(left_antenna_control(i).xveloc_in_mm...
            (windowOnLeft-20000:windowOnLeft),'omitnan');
        %Sum spikes 1st second of poke
        xVelocityDuringPokeControlRight{i}(k) = mean(right_antenna_control(i).xveloc_in_mm...
            (windowOnRight:windowOnRight+analysis_window_size*sampling_rate),'omitnan');
        xVelocityDuringPokeControlLeft{i}(k) = mean(left_antenna_control(i).xveloc_in_mm...
            (windowOnLeft:windowOnLeft+analysis_window_size*sampling_rate),'omitnan');
        %Sum spikes after poke
        xVelocityPostPokeControlRight{i}(k)  = mean(right_antenna_control(i).xveloc_in_mm...
            (windowOffRight:windowOffRight+20000),'omitnan');
        xVelocityPostPokeControlLeft{i}(k)  = mean(left_antenna_control(i).xveloc_in_mm...
            (windowOffLeft:windowOffLeft+20000),'omitnan');

    end
end

%Calculate means of X velocity
meanPreXRightControl = cellfun(@mean, xVelocityPrePokeControlRight);
meanPokeXRightControl = cellfun(@mean, xVelocityDuringPokeControlRight);
meanPostXRightControl = cellfun(@mean, xVelocityPostPokeControlRight);

meanPreXLeftControl = cellfun(@mean, xVelocityPrePokeControlLeft);
meanPokeXLeftControl = cellfun(@mean, xVelocityDuringPokeControlLeft);
meanPostXLeftControl = cellfun(@mean, xVelocityPostPokeControlLeft);

%Statistics
[pWilPreVsPokeRightControl] = signrank(meanPreXRightControl, meanPokeXRightControl);
[pWilPostVsPokeRightControl] = signrank(meanPostXRightControl, meanPokeXRightControl);

[pWilPreVsPokeLeftControl] = signrank(meanPreXLeftControl, meanPokeXLeftControl);
[pWilPostVsPokeLeftControl] = signrank(meanPostXLeftControl, meanPokeXLeftControl);

%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean firing rates pre, during 1st  500 ms and post stimulus -- RIGHT')
t = tiledlayout('flow');
title(t, 'CONTROL')
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna_control)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXRightControl; meanPokeXRightControl; meanPostXRightControl];
h = boxplot2(plot_y',groups);
xticks([1,2,3]);
xticklabels(groupNames);
hold on
plot_y = [meanPreXRightControl; meanPokeXRightControl; meanPostXRightControl];
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,1, ['p=' num2str(pWilPreVsPokeRightControl)])
text(3,1, ['p=' num2str(pWilPostVsPokeRightControl)])
ylim([-1 1])
ylabel('X Velocity (mm/s)')

medianXvelocControlRight = h.med(2).YData(1);

nexttile
title("LEFT" + newline + "N =" + num2str(length(left_antenna_control)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXLeftControl; meanPokeXLeftControl; meanPostXLeftControl];
h = boxplot2(plot_y',groups);
xticks([1,2,3]);
xticklabels(groupNames);
hold on
plot_y = [meanPreXLeftControl; meanPokeXLeftControl; meanPostXLeftControl];
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,1, ['p=' num2str(pWilPreVsPokeLeftControl)])
text(3,1, ['p=' num2str(pWilPostVsPokeLeftControl)])
ylim([-1 1])
ylabel('X Velocity (mm/s)')

medianXvelocControlLeft = h.med(2).YData(1);

% print(['X_control_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['X_control_boxplot' '.png'], '-dpng','-r300', '-vector')


%% Calculate mean Z velocity and create boxplots for CONTROL

%Size of poking analysis window in s
analysis_window_size = 0.5;


for i = 1 : length(right_antenna_control)
    for k = 1 : numPokes
        windowOnRight = right_antenna_control(i).trigger_on(k);
        windowOffRight = right_antenna_control(i).trigger_off(k);
        windowOnLeft = left_antenna_control(i).trigger_on(k);
        windowOffLeft = left_antenna_control(i).trigger_off(k);
        %Sum spikes 1 s pre poke
        zVelocityPrePokeControlRight{i}(k) = mean(right_antenna_control(i).zveloc_in_degree_per_s...
            (windowOnRight-20000:windowOnRight),'omitnan');
        zVelocityPrePokeControlLeft{i}(k) = mean(left_antenna_control(i).zveloc_in_degree_per_s...
            (windowOnLeft-20000:windowOnLeft),'omitnan');
        %Sum spikes 1st second of poke
        zVelocityDuringPokeControlRight{i}(k) = mean(right_antenna_control(i).zveloc_in_degree_per_s...
            (windowOnRight:windowOnRight+analysis_window_size*sampling_rate),'omitnan');
        zVelocityDuringPokeControlLeft{i}(k) = mean(left_antenna_control(i).zveloc_in_degree_per_s...
            (windowOnLeft:windowOnLeft+analysis_window_size*sampling_rate),'omitnan');
        %Sum spikes after poke
        zVelocityPostPokeControlRight{i}(k)  = mean(right_antenna_control(i).zveloc_in_degree_per_s...
            (windowOffRight:windowOffRight+20000),'omitnan');
        zVelocityPostPokeControlLeft{i}(k)  = mean(left_antenna_control(i).zveloc_in_degree_per_s...
            (windowOffLeft:windowOffLeft+20000),'omitnan');

    end
end

%Calculate means of X velocity
meanPreZRightControl = cellfun(@mean, zVelocityPrePokeControlRight);
meanPokeZRightControl = cellfun(@mean, zVelocityDuringPokeControlRight);
meanPostZRightControl = cellfun(@mean, zVelocityPostPokeControlRight);

meanPreZLeftControl = cellfun(@mean, zVelocityPrePokeControlLeft);
meanPokeZLeftControl = cellfun(@mean, zVelocityDuringPokeControlLeft);
meanPostZLeftControl = cellfun(@mean, zVelocityPostPokeControlLeft);

%Statistics
[pWilPreVsPokeZRightControl] = signrank(meanPreZRightControl, meanPokeZRightControl)
[pWilPostVsPokeZRightControl] = signrank(meanPostZRightControl, meanPokeZRightControl)

[pWilPreVsPokeZLeftControl] = signrank(meanPreZLeftControl, meanPokeZLeftControl)
[pWilPostVsPokeZLeftControl] = signrank(meanPostZLeftControl, meanPokeZLeftControl)

%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean firing rates pre, during 1st  500 ms and post stimulus -- RIGHT')
t = tiledlayout('flow');
title(t, 'CONTROL')
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna_control)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreZRightControl; meanPokeZRightControl; meanPostZRightControl];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
plot_y = [meanPreZRightControl; meanPokeZRightControl; meanPostZRightControl];
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,1, ['p=' num2str(pWilPreVsPokeZRightControl)])
text(3,1, ['p=' num2str(pWilPostVsPokeZRightControl)])
ylim([-1 1])
ylabel('Z Velocity (mm/s)')


nexttile
title("LEFT" + newline + "N =" + num2str(length(left_antenna_control)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreZLeftControl; meanPokeZLeftControl; meanPostZLeftControl];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
plot_y = [meanPreZLeftControl; meanPokeZLeftControl; meanPostZLeftControl];
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,1, ['p=' num2str(pWilPreVsPokeZLeftControl)])
text(3,1, ['p=' num2str(pWilPostVsPokeZLeftControl)])
ylim([-1 1])
ylabel('Z Velocity (mm/s)')

% print(['Z_control_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['Z_control_boxplot' '.png'], '-dpng','-r300', '-vector')

%% Calculate mean X velocitx and create boxplots for SILENCED

clearvars xVelocityPrePokeRight xVelocityPrePokeLeft xVelocityDuringPokeRight...
    xVelocityDuringPokeLeft xVelocityDuringPokeLeft xVelocityPostPokeLeft...
    xVelocityPostPokeRight

%Size of poking analysis window in s
analysis_window_size = 0.5;

for i = 1 : length(right_antenna_kir)
    for k = 1 : numPokes
        windowOnRight = right_antenna_kir(i).trigger_on(k);
        windowOffRight = right_antenna_kir(i).trigger_off(k);
        windowOnLeft = left_antenna_kir(i).trigger_on(k);
        windowOffLeft = left_antenna_kir(i).trigger_off(k);
        %Sum spikes 1 s pre poke
        xVelocityPrePokeKirRight{i}(k) = mean(right_antenna_kir(i).xveloc_in_mm...
            (windowOnRight-20000:windowOnRight),'omitnan');
        xVelocityPrePokeKirLeft{i}(k) = mean(left_antenna_kir(i).xveloc_in_mm...
            (windowOnLeft-20000:windowOnLeft),'omitnan');
        %Sum spikes 1st second of poke
        xVelocityDuringPokeKirRight{i}(k) = mean(right_antenna_kir(i).xveloc_in_mm...
            (windowOnRight:windowOnRight+analysis_window_size*sampling_rate),'omitnan');
        xVelocityDuringPokeKirLeft{i}(k) = mean(left_antenna_kir(i).xveloc_in_mm...
            (windowOnLeft:windowOnLeft+analysis_window_size*sampling_rate),'omitnan');
        %Sum spikes after poke
        xVelocityPostPokeKirRight{i}(k)  = mean(right_antenna_kir(i).xveloc_in_mm...
            (windowOffRight:windowOffRight+20000),'omitnan');
        xVelocityPostPokeKirLeft{i}(k)  = mean(left_antenna_kir(i).xveloc_in_mm...
            (windowOffLeft:windowOffLeft+20000),'omitnan');
    end
end

%Calculate means of X velocity
meanPreXRightKir = cellfun(@mean, xVelocityPrePokeKirRight);
meanPokeXRightKir = cellfun(@mean, xVelocityDuringPokeKirRight);
meanPostXRightKir = cellfun(@mean, xVelocityPostPokeKirRight);

meanPreXLeftKir = cellfun(@mean, xVelocityPrePokeKirLeft);
meanPokeXLeftKir = cellfun(@mean, xVelocityDuringPokeKirLeft);
meanPostXLeftKir = cellfun(@mean, xVelocityPostPokeKirLeft);

%Statistics
[pWilPreVsPokeXRightKir] = signrank(meanPreXRightKir, meanPokeXRightKir);
[pWilPostVsPokeXRightKir] = signrank(meanPostXRightKir, meanPokeXRightKir);

[pWilPreVsPokeXLeftKir] = signrank(meanPreXLeftKir, meanPokeXLeftKir);
[pWilPostVsPokeXLeftKir] = signrank(meanPostXLeftKir, meanPokeXLeftKir);


%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean firing rates pre, during 1st  500 ms and post stimulus')
t = tiledlayout('flow');
title(t, 'SILENCED')
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna_kir)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXRightKir; meanPokeXRightKir; meanPostXRightKir];
h = boxplot2(plot_y',groups);
xticks([1,2,3]);
xticklabels(groupNames);
hold on
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,1, ['p=' num2str(pWilPreVsPokeXRightKir)])
text(3,1, ['p=' num2str(pWilPostVsPokeXRightKir)])
ylim([-1 1])
ylabel('X Velocity (mm/s)')

medianXvelocKirRight = h.med(2).YData(1);


nexttile
title("LEFT" + newline + "N =" + num2str(length(left_antenna_kir)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXLeftKir; meanPokeXLeftKir; meanPostXLeftKir];
h = boxplot2(plot_y',groups);
xticks([1,2,3]);
xticklabels(groupNames);
hold on
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
ylim([-1 1])
ylabel('X Velocity (mm/s)')
text(1,1, ['p=' num2str(pWilPreVsPokeXLeftKir)])
text(3,1, ['p=' num2str(pWilPostVsPokeXLeftKir)])

medianXvelocKirRight = h.med(2).YData(1);

% print(['X_silenced_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['X_silenced_boxplot' '.png'], '-dpng','-r300', '-vector')

%% Calculate mean Z velocity and create boxplots for SILENCED

clearvars xVelocityPrePokeRight xVelocityPrePokeLeft xVelocityDuringPokeRight...
    xVelocityDuringPokeLeft xVelocityDuringPokeLeft xVelocityPostPokeLeft...
    xVelocityPostPokeRight

%Size of poking analysis window in s
analysis_window_size = 0.5;

for i = 1 : length(right_antenna_kir)
    for k = 1 : numPokes
        windowOnRight = right_antenna_kir(i).trigger_on(k);
        windowOffRight = right_antenna_kir(i).trigger_off(k);
        windowOnLeft = left_antenna_kir(i).trigger_on(k);
        windowOffLeft = left_antenna_kir(i).trigger_off(k);
        %Sum spikes 1 s pre poke
        zVelocityPrePokeKirRight{i}(k) = mean(right_antenna_kir(i).zveloc_in_degree_per_s...
            (windowOnRight-20000:windowOnRight),'omitnan');
        zVelocityPrePokeKirLeft{i}(k) = mean(left_antenna_kir(i).zveloc_in_degree_per_s...
            (windowOnLeft-20000:windowOnLeft),'omitnan');
        %Sum spikes 1st second of poke
        zVelocityDuringPokeKirRight{i}(k) = mean(right_antenna_kir(i).zveloc_in_degree_per_s...
            (windowOnRight:windowOnRight+analysis_window_size*sampling_rate),'omitnan');
        zVelocityDuringPokeKirLeft{i}(k) = mean(left_antenna_kir(i).zveloc_in_degree_per_s...
            (windowOnLeft:windowOnLeft+analysis_window_size*sampling_rate),'omitnan');
        %Sum spikes after poke
        zVelocityPostPokeKirRight{i}(k)  = mean(right_antenna_kir(i).zveloc_in_degree_per_s...
            (windowOffRight:windowOffRight+20000),'omitnan');
        zVelocityPostPokeKirLeft{i}(k)  = mean(left_antenna_kir(i).zveloc_in_degree_per_s...
            (windowOffLeft:windowOffLeft+20000),'omitnan');
    end
end

%Calculate means of X velocity
meanPreZRightKir = cellfun(@mean, zVelocityPrePokeKirRight);
meanPokeZRightKir = cellfun(@mean, zVelocityDuringPokeKirRight);
meanPostZRightKir = cellfun(@mean, zVelocityPostPokeKirRight);

meanPreZLeftKir = cellfun(@mean, zVelocityPrePokeKirLeft);
meanPokeZLeftKir = cellfun(@mean, zVelocityDuringPokeKirLeft);
meanPostZLeftKir = cellfun(@mean, zVelocityPostPokeKirLeft);

%Statistics
[pWilPreVsPokeZRightKir] = signrank(meanPreZRightKir, meanPokeZRightKir);
[pWilPostVsPokeZRightKir] = signrank(meanPostZRightKir, meanPokeZRightKir);

[pWilPreVsPokeZLeftKir] = signrank(meanPreZLeftKir, meanPokeZLeftKir);
[pWilPostVsPokeZLeftKir] = signrank(meanPostZLeftKir, meanPokeZLeftKir);


%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean firing rates pre, during 1st  500 ms and post stimulus')
t = tiledlayout('flow');
title(t, 'SILENCED')
nexttile
title("RIGHT" + newline + "N =" + num2str(length(right_antenna_kir)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreZRightKir; meanPokeZRightKir; meanPostZRightKir];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,1, ['p=' num2str(pWilPreVsPokeZRightKir)])
text(3,1, ['p=' num2str(pWilPostVsPokeZRightKir)])
ylim([-1 1])
ylabel('Z Velocity (°/s)')


nexttile
title("LEFT" + newline + "N =" + num2str(length(left_antenna_kir)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreZLeftKir; meanPokeZLeftKir; meanPostZLeftKir];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
ylim([-1 1])
ylabel('Z Velocity (°/s)')
text(1,1, ['p=' num2str(pWilPreVsPokeZLeftKir)])
text(3,1, ['p=' num2str(pWilPostVsPokeZLeftKir)])

% print(['Z_silenced_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['Z_silenced_boxplot' '.png'], '-dpng','-r300', '-vector')


%% Statistics test X Velocity SILENCED vs CONTROL

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanPokeXRightControlPlot, meanPokeXRightKirPlot ] = padarrays(meanPokeXRightControl, meanPokeXRightKir);
[meanPokeXLeftControlPlot, meanPokeXLeftKirPlot ] = padarrays(meanPokeXLeftControl, meanPokeXLeftKir);

[pWilControlvsSilencedPokeXRight] = ranksum(meanPokeXRightControl, meanPokeXRightKir)
[pWilControlvsSilencedPokeXLeft] = ranksum(meanPokeXLeftControl, meanPokeXLeftKir)

figure
t = tiledlayout('flow');
nexttile;
title('Poke Control vs Silenced Right')
groupNames = ["Control "+'N = ' + num2str(length(right_antenna_control));...
    "Silenced "+'N = ' + num2str(length(right_antenna_kir))];
groups = [1; 2];
plot_y = [meanPokeXRightControlPlot; meanPokeXRightKirPlot];
boxplot2(plot_y', groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
ylabel('X Velocity (mm/s)')
ylim([-1.5 1.5])
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeXRight)])


nexttile
title('Poke Control vs Silenced Left')
groupNames = ["Control "+'N = ' + num2str(length(left_antenna_control));...
    "Silenced "+'N = ' + num2str(length(left_antenna_kir))];
groups = [1; 2];
plot_y = [meanPokeXLeftControlPlot; meanPokeXLeftKirPlot];
boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
ylabel('X Velocity (mm/s)')
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeXLeft)])
ylim([-1.5 1.5])

title(t, "NEW Dataset")

% print(['X_control_vs_silenced' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['X_control_vs_silenced' '.png'], '-dpng','-r300', '-vector')

% Statistics Control Left vs Right

[pWilControlLeftvsRightPokeX] = signrank(meanPokeXRightControl, meanPokeXLeftControl)

figure
t = tiledlayout('flow');
nexttile;
title('Poke Control Left vs Right')
groupNames = ["Left "+'N = ' + num2str(length(left_antenna_control));...
    "Right "+'N = ' + num2str(length(right_antenna_control))];
groups = [1; 2];
plot_y = [meanPokeXLeftControl; meanPokeXRightControl];
h = boxplot2(plot_y', groups);
xticks([1,2]);
xticklabels(groupNames);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
ylabel('X Velocity (mm/s)')
ylim([-1.5 1.5])
yline(0,'--');

text(1,1, ['p=' num2str(pWilControlLeftvsRightPokeX)])
medianXvelocLeftVsRight = h.med(1).YData(1);
medianXvelocLeftVsRight = h.med(2).YData(1);

% print(['X_control_left_vs_right' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['X_control_left_vs_right' '.png'], '-dpng','-r300', '-vector')
%% Statistics test Z Velocity SILENCED vs CONTROL

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanPokeZRightControlPlot, meanPokeZRightKirPlot ] = padarrays(meanPokeZRightControl, meanPokeZRightKir);
[meanPokeZLeftControlPlot, meanPokeZLeftKirPlot ] = padarrays(meanPokeZLeftControl, meanPokeZLeftKir);

[pWilControlvsSilencedPokeZRight] = ranksum(meanPokeZRightControl, meanPokeZRightKir);
[pWilControlvsSilencedPokeZLeft] = ranksum(meanPokeZLeftControl, meanPokeZLeftKir);

figure
t = tiledlayout('flow');
nexttile;
title('Poke Control vs Silenced Right')
groupNames = ["Control "+'N = ' + num2str(length(right_antenna_control));...
    "Silenced "+'N = ' + num2str(length(right_antenna_kir))];
groups = [1; 2];
plot_y = [meanPokeZRightControlPlot; meanPokeZRightKirPlot];
boxplot2(plot_y', groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
ylabel('Z Velocity (mm/s)')
ylim([-1.5 1.5])
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeZRight)])


nexttile
title('Poke Control vs Silenced Left')
groupNames = ["Control "+'N = ' + num2str(length(left_antenna_control));...
    "Silenced "+'N = ' + num2str(length(left_antenna_kir))];
groups = [1; 2];
plot_y = [meanPokeZLeftControlPlot; meanPokeZLeftKirPlot];
boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
ylabel('Z Velocity (mm/s)')
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeZLeft)])
ylim([-1.5 1.5])

title(t, "NEW Dataset")

% print(['Z_control_vs_silenced' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['Z_control_vs_silenced' '.png'], '-dpng','-r300', '-vector')

% Statistics Control Left vs Right

[pWilControlLeftvsRightPokeZ] = signrank(meanPokeZRightControl*(-1), meanPokeZLeftControl);

figure
t = tiledlayout('flow');
nexttile;
title('Poke Control Left vs Right')
groupNames = ["Left "+'N = ' + num2str(length(left_antenna_control));...
    "Right "+'N = ' + num2str(length(right_antenna_control))];
groups = [1; 2];
plot_y = [meanPokeZLeftControl; meanPokeZRightControl*(-1)];
h = boxplot2(plot_y', groups);
xticks([1,2]);
xticklabels(groupNames);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
ylabel('Z Velocity (°/s)')
% ylim([-1.5 1.5])
yline(0,'--');

text(1,1, ['p=' num2str(pWilControlLeftvsRightPokeZ)])
medianZvelocLeftVsLeft = h.med(1).YData(1);
medianZvelocLeftVsRight = h.med(2).YData(1);

 % print(['Z_control_left_vs_right' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
 % print(['Z_control_left_vs_right' '.png'], '-dpng','-r300', '-vector')
%% POOL RIGHT and LEFT antenna X and compare CONTROL vs SILENCED

%Pool means 1st 500ms of poke for right and left antenna
meanXControlPooled = [meanPokeXRightControl meanPokeXLeftControl];
meanXKirPooled = [meanPokeXRightKir meanPokeXLeftKir];

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanControlPooledPlot, meanKirPooledPlot ] = padarrays(meanXControlPooled, meanXKirPooled);

%Ranksum test for control vs silenced pooled 1st 500ms of poke
[pWilControlvsSilencedPoke] = ranksum(meanXControlPooled, meanXKirPooled);

%Plot means in a boxplot, control vs silenced
figure
t = tiledlayout('flow');
title(t, 'Poke Control vs Silenced POOLED')
nexttile
groupNames = ["Control " + "N = " + num2str(length(meanXControlPooled));...
    "Silenced " + "N = "+ num2str(length(meanXKirPooled))];
groups = [1; 2];
plot_y = [meanControlPooledPlot; meanKirPooledPlot];
boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled', 'jitter','on', 'jitterAmount', 0.1)
ylabel('X Velocity (mm/s)')
ylim([-1.5 1.5])
text(1.5,1, ['p=' num2str(pWilControlvsSilencedPoke)])

%Pool x vel curve for right and left  plot control vs silenced
meanXCurveControlPooled = [meanXvelocCellControlRight meanXvelocCellControlLeft];
meanXCurveOverallControl = mean(cat(1, meanXCurveControlPooled{:}));

meanXCurveKirPooled = [meanXvelocCellKirRight meanXvelocCellKirLeft];
meanXCurveOverallKir = mean(cat(1, meanXCurveKirPooled{:}));

nexttile
plot(xbin,meanXCurveOverallControl, 'LineWidth', 3)
hold on
plot(xbin,meanXCurveOverallKir, 'LineWidth', 3)
rectangle('Position',[1 -0.5 2 1], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylabel('X Velocity (mm/s)')
xlabel('Time (s)')
legend('control', 'silenced','Location','northwest')


% print(['X_pooled_control_vs_silenced_line' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['X_pooled_control_vs_silenced_line' '.png'], '-dpng','-r300', '-vector')

%% POOL RIGHT and LEFT antenna Z and compare CONTROL vs SILENCED

%Pool means 1st 500ms of poke for right and left antenna
meanZControlPooled = [meanPokeZRightControl meanPokeZLeftControl];
meanZKirPooled = [meanPokeZRightKir meanPokeZLeftKir];

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanZControlPooledPlot, meanZKirPooledPlot ] = padarrays(meanZControlPooled, meanZKirPooled);

%Ranksum test for control vs silenced pooled 1st 500ms of poke
[pWilControlvsSilencedZPoke] = ranksum(meanZControlPooled, meanZKirPooled);

%Plot means in a boxplot, control vs silenced
figure
t = tiledlayout('flow');
title(t, 'Poke Control vs Silenced POOLED')
nexttile
groupNames = ["Control " + "N = " + num2str(length(meanZControlPooled));...
    "Silenced " + "N = "+ num2str(length(meanZKirPooled))];
groups = [1; 2];
plot_y = [meanZControlPooledPlot; meanZKirPooledPlot];
boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled', 'jitter','on', 'jitterAmount', 0.1)
ylabel('Z Velocity (mm/s)')
ylim([-1.5 1.5])
text(1.5,1, ['p=' num2str(pWilControlvsSilencedZPoke)])

% print(['Z_pooled_control_vs_silenced_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')

%Pool x vel curve for right and left  plot control vs silenced
meanZCurveControlPooled = [meanZvelocCellControlRight meanZvelocCellControlLeft];
meanZCurveOverallControl = mean(cat(1, meanZCurveControlPooled{:}));

meanZCurveKirPooled = [meanZvelocCellKirRight meanZvelocCellKirLeft];
meanZCurveOverallKir = mean(cat(1, meanZCurveKirPooled{:}));

nexttile
plot(xbin,meanZCurveOverallControl, 'LineWidth', 3)
hold on
plot(xbin,meanZCurveOverallKir, 'LineWidth', 3)
rectangle('Position',[1 -0.5 2 1], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylabel('Z Velocity (mm/s)')
xlabel('Time (s)')
legend('control', 'silenced','Location','northwest')


% print(['Z_pooled_control_vs_silenced_line' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% print(['Z_pooled_control_vs_silenced_line' '.png'], '-dpng','-r300', '-vector')

%% Sort for no walking before poking

%Set x  velocity threshold and window size
windowBeforeInSec = 1;
walkingThreshold = 10;

for i = 1 : length(right_antenna_control)
    for k = 1 : numPokes
        windowOffRight = right_antenna_control(i).trigger_on(k);
        windowOnRight = right_antenna_control(i).trigger_on(k)-...
            (sampling_rate*windowBeforeInSec);
        XVelocitybefore{i,k} = right_antenna_control(i).xveloc_in_mm...
            (windowOnRight:windowOffRight);
    end
end

meanXVelocitybefore = cellfun(@mean, XVelocitybefore);
medianXVelocitybefore = cellfun(@mean, XVelocitybefore);

%Find index of trials exceeding the treshold
logicalMask = [meanXVelocitybefore > walkingThreshold];

%Create duplicate of data for filtering
xVelocityDuringPokeControlRightFiltered = xVelocityDuringPokeControlRight;

%Go through the logical mask and remove values in the X velocity where
%logcialMask == 1 (threshold is exceeded)
for i = 1 : size(logicalMask, 1)
    for k = 1 : size(logicalMask, 2)
        if logicalMask(i,k) == 1
            xVelocityDuringPokeControlRightFiltered{i}(k) = nan;
        else
        end
    end
end

meanPokeRightControlFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    xVelocityDuringPokeControlRightFiltered);

clearvars medianXVelocitybefore XVelocitybefore meanXVelocitybefore...
    logicalMask

for i = 1 : length(left_antenna_control)
    for k = 1 : numPokes
        windowOffLeft = left_antenna_control(i).trigger_on(k);
        windowOnLeft = left_antenna_control(i).trigger_on(k)-...
            (sampling_rate*windowBeforeInSec);
        XVelocitybefore{i,k} = left_antenna_control(i).xveloc_in_mm...
            (windowOnLeft:windowOffLeft);
    end
end

meanXVelocitybefore = cellfun(@mean, XVelocitybefore);
medianXVelocitybefore = cellfun(@mean, XVelocitybefore);

%Find index of trials exceeding the treshold
logicalMask = [meanXVelocitybefore > walkingThreshold];

%Create duplicate of data for filtering
xVelocityDuringPokeControlLeftFiltered = xVelocityDuringPokeControlLeft;

%Go through the logical mask and remove values in the X velocity where
%logcialMask == 1 (threshold is exceeded)
for i = 1 : size(logicalMask, 1)
    for k = 1 : size(logicalMask, 2)
        if logicalMask(i,k) == 1
            xVelocityDuringPokeControlLeftFiltered{i}(k) = nan;
        else
        end
    end
end

meanPokeLeftControlFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    xVelocityDuringPokeControlLeftFiltered);


%% Statistics test SILENCED vs CONTROL with x velocity threshold fitered data

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanPokeRightControlPlot, meanPokeRightKirPlot ] = padarrays(meanPokeRightControlFiltered, meanPokeXRightKir);
[meanPokeLeftControlPlot, meanPokeLeftKirPlot ] = padarrays(meanPokeLeftControlFiltered, meanPokeXLeftKir);



[pWilControlvsSilencedPokeRight] = ranksum(meanPokeRightControlFiltered, meanPokeXRightKir);
[pWilControlvsSilencedPokeLeft] = ranksum(meanPokeLeftControlFiltered, meanPokeXLeftKir);

figure
t = tiledlayout('flow');
nexttile;
title('Poke Control vs Silenced Right')
groupNames = ["Control "+'N = ' + num2str(length(right_antenna_control));...
    "Silenced "+'N = ' + num2str(length(right_antenna_kir))];
groups = [1; 2];
plot_y = [meanPokeRightControlPlot; meanPokeRightKirPlot];
boxplot2(plot_y', groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
ylabel('X Velocity (mm/s)')
ylim([-1.5 1.5])
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeRight)])


nexttile
title('Poke Control vs Silenced Left')
groupNames = ["Control "+'N = ' + num2str(length(left_antenna_control));...
    "Silenced "+'N = ' + num2str(length(left_antenna_kir))];
groups = [1; 2];
plot_y = [meanPokeLeftControlPlot; meanPokeLeftKirPlot];
boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
ylabel('X Velocity (mm/s)')
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeLeft)])
ylim([-1.5 1.5])

title(t, "NEW Dataset")

% print(['control_vs_silenced' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
