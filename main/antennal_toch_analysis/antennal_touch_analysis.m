
%% Load preprocessed data and set analysis parameter

load('data/preprocessed.mat')

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

for k = 1 : length(analysis)
    analysis(k).xveloc_in_mm = analysis(k).VelocX(:,1)*8.79;
    analysis(k).zveloc_in_degree_per_s = analysis(k).VelocZ(:,1)*158.9;
end
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
    y = T(1:end,3);
    z = zeros(size(x));

    scatter(x(1), y(1), [3000], '.')

    %c = analysis(i).spikesnorm_binned;

    scatter(x(1:end-2), y(1:end-2), [], 'c',  'filled', 'o');

    %scatter(x(1:end-2), y(1:end-2), [], 'filled', 'o');

    colormap(flipud(cool))
    colorbar;
    caxis([0 1])


    plot(0,0,...
        'Marker','o', 'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    % Total distance walked
    d = hypot(diff(x), diff(y));                            % Distance Of Each Segment
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

clearvars binstart binstarts binsize binned_xveloc_right

binstarts = nan;
binsize = binsize_factor*sampling_rate;

%binned_xveloc_right = cell(1,length(analysis));

for k = 1 : length(analysis)

    temp_xveloc = [];

    for m = 1 : numPokes
        window_on = analysis(k).trigger_on(m)-1*sampling_rate; %2s before poking
        window_off = analysis(k).trigger_off(m)+1*sampling_rate; %2s after poking

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

    end
    analysis(k).binnedZveloc = cellfun(@(x)x/binsize, temp_zveloc);
    analysis(k).binnedZvelocMean = mean(analysis(k).binnedZveloc);
    analysis(k).binnedZvelocDegree = cellfun(@(x)x/binsize, temp_zveloc_degree);
    analysis(k).binnedZvelocDegreeMean = mean(analysis(k).binnedZvelocDegree);
end

%% Y velocity preliminary

for k = 1 : length(analysis)
    analysis(k).yveloc_in_mm = analysis(k).VelocY(:,1)*8.79;
end

%%%%% Right Antenna
clearvars binstart binstarts binsize

binstarts = nan;
binsize = binsize_factor*sampling_rate;

%binned_zveloc_right = cell(1,length(analysis));

for k = 1 : length(analysis)

    temp_zveloc = [];

    for m = 1 : numPokes
        window_on = analysis(k).trigger_on(m)-1*sampling_rate; %2s before poking
        window_off = analysis(k).trigger_off(m)+1*sampling_rate; %2s after poking

        binstarts = window_on:binsize:window_off;
        for bin = 1 : length(binstarts)-1
            temp_yveloc{m,bin} = sum(analysis(k).yveloc_in_mm(binstarts(bin):binstarts(bin+1)),'omitnan');
        end
    end
    analysis(k).binnedYveloc = cellfun(@(x)x/binsize, temp_yveloc);
    analysis(k).binnedYvelocMean = mean(analysis(k).binnedYveloc);
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
save('postPoker_right_control_pooled','right_antenna_control','-v7.3');
save('postPoker_left_control_pooled','left_antenna_control','-v7.3');

save('postPoker_right_kir_pooled','right_antenna_kir','-v7.3');
save('postPoker_left_kir_pooled','left_antenna_kir','-v7.3');

save('postAnalysis_pooled', 'analysis', '-v7.3');

%% Quality control

clearvars right_antenna_control_min_velocity left_antenna_control_min_velocity

% Exlucde trials with low walking activity 
right_antenna_control_min_velocity = right_antenna_control(~cellfun(@(x)x<0.5,{right_antenna_control.mean_velocity}));
left_antenna_control_min_velocity  = left_antenna_control(~cellfun(@(x)x<0.5,{left_antenna_control.mean_velocity}));


%% Plotting Forward and Angular Velocitys during stimulation with 2 s before and 2 s after ---- CONTROL
clearvars binned_spikerate_right_mean binned_spikerate_left_mean binned_xveloc_right_mean...
    binned_xveloc_left_mean binned_zveloc_right_mean binned_zveloc_left_mean xbin


meanXvelocCellControlRight = arrayfun(@(s) s.binnedXvelocMean, right_antenna_control, 'UniformOutput', false);
meanXvelocOverallControlRight = mean(cat(1, meanXvelocCellControlRight{:}));

meanXvelocCellControlLeft = arrayfun(@(s) s.binnedXvelocMean, left_antenna_control, 'UniformOutput', false);
meanXvelocOverallControlLeft = mean(cat(1, meanXvelocCellControlLeft{:}));

meanXvelocCellControlPooled = [meanXvelocCellControlRight meanXvelocCellControlLeft];
meanXvelocOverallControlPooled = mean(cat(1, meanXvelocCellControlPooled{:}));

meanZvelocCellControlRight = arrayfun(@(s) s.binnedZvelocDegreeMean, right_antenna_control, 'UniformOutput', false);
meanZvelocOverallControlRight = mean(cat(1, meanZvelocCellControlRight{:}));

meanZvelocCellControlLeft = arrayfun(@(s) s.binnedZvelocDegreeMean, left_antenna_control, 'UniformOutput', false);
meanZvelocOverallControlLeft = mean(cat(1, meanZvelocCellControlLeft{:}));

% Mirror right antenna *(-1)
meanZvelocCellControlRightMirrored = cellfun(@(x) x*(-1),...
    meanZvelocCellControlRight, 'UniformOutput', false);

meanZvelocCellControlPooled = [meanZvelocCellControlRightMirrored meanZvelocCellControlLeft];
meanZvelocOverallControlPooled = mean(cat(1, meanZvelocCellControlPooled{:}));


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
% ylim([-2.5 2.5])
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
% ylim([-2.5 2.5])
xlabel('Time (s)')

% %save the figure as a high-quality PNG image
%%%print(gcf,'poke1.png','-dpng','-r300');

figure('name', 'Overview')
t = tiledlayout('flow');
title(t,'CONTROL');
nexttile

% calculate SEM of x velocity for right and left antenna
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

%plot Left antenna x velocity with SEM as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_X_Control_Left fliplr(y_lower_limit_X_Control_Left)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallControlLeft, 'LineWidth',3)

%plot Left antenna x velocity with SEM as shaded area
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

% calculate SEM of Z velocity for right and left antenna

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

%Add SD to y upper and limit and subtract from lower limit
y_upper_limit_Z_Control_Left = meanZvelocOverallControlLeft + stdev_left;
y_lower_limit_Z_Control_Left = meanZvelocOverallControlLeft - stdev_left;
y_upper_limit_Z_Control_Right = meanZvelocOverallControlRight + stdev_right;
y_lower_limit_Z_Control_Right = meanZvelocOverallControlRight - stdev_right;
%plot Left antenna x velocity with SEM as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Control_Left fliplr(y_lower_limit_Z_Control_Left)],...
    [1 0 0], 'FaceAlpha', 0.3,'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallControlLeft, 'LineWidth',3)

%plot Left antenna x velocity with SEM as shaded area
hold on
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Control_Right fliplr(y_lower_limit_Z_Control_Right)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanZvelocOverallControlRight, 'LineWidth',3)
title('Z velocity')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
% ylim([-1 2])
ylabel('Z Velocity (mm))');
xlabel('Time (s)');

% %save the figure as a high-quality PNG image
%%%print(gcf,'poke2.png','-dpng','-r300');

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

% %%print(['control_overview' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['control_overview' '.png'], '-dpng','-r300', '-vector')

%% Plotting Forward and Angular Velocitys during stimulation with 2 s before and 2 s after ---- SILENCED
clearvars binned_spikerate_right_mean binned_spikerate_left_mean binned_xveloc_right_mean...
    binned_xveloc_left_mean binned_zveloc_right_mean binned_zveloc_left_mean xbin


meanXvelocCellKirRight = arrayfun(@(s) s.binnedXvelocMean, right_antenna_kir, 'UniformOutput', false);
meanXvelocOverallKirRight = mean(cat(1, meanXvelocCellKirRight{:}));

meanXvelocCellKirLeft = arrayfun(@(s) s.binnedXvelocMean, left_antenna_kir, 'UniformOutput', false);
meanXvelocOverallKirLeft = mean(cat(1, meanXvelocCellKirLeft{:}));

meanXvelocCellKirPooled = [meanXvelocCellKirRight meanXvelocCellKirLeft];
meanXvelocOverallKirPooled = mean(cat(1, meanXvelocCellKirPooled{:}));


meanZvelocCellKirRight = arrayfun(@(s) s.binnedZvelocDegreeMean, right_antenna_kir, 'UniformOutput', false);
meanZvelocOverallKirRight = mean(cat(1, meanZvelocCellKirRight{:}));

meanZvelocCellKirLeft = arrayfun(@(s) s.binnedZvelocDegreeMean, left_antenna_kir, 'UniformOutput', false);
meanZvelocOverallKirLeft = mean(cat(1, meanZvelocCellKirLeft{:}));

% Mirror right antenna *(-1)
meanZvelocCellKirRightMirrored = cellfun(@(x) x*(-1),...
    meanZvelocCellKirRight, 'UniformOutput', false);

meanZvelocCellKirPooled = [meanZvelocCellKirRightMirrored meanZvelocCellKirLeft];
meanZvelocOverallKirPooled = mean(cat(1, meanZvelocCellKirPooled{:}));



%xbin = [1:bin]/(sampling_rate/binsize);
xbin = ((binsize/2)/sampling_rate):binsize/sampling_rate:(bin*(binsize/sampling_rate));

% Create a colormap
cmap = hsv(numel(left_antenna_kir));

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
% %save the figure as a high-quality PNG image
%%%print(gcf,'poke1.png','-dpng','-r300');
figure('name', 'Overview')
t = tiledlayout('flow');
title(t,'SILENCED');
nexttile

% calculate SEM of x velocity for right and left antenna
stdev_left = std(cat(1, meanXvelocCellKirLeft{:}));
stderror_left = stdev_left/sqrt(length(cat(1, meanXvelocCellKirLeft{:})));
stdev_right = std(cat(1, meanXvelocCellKirRight{:}));
stderror_right = stdev_right/sqrt(length(cat(1, meanXvelocCellKirRight{:})));

% %Add SEM to y upper and limit and subtract from lower limit
% y_upper_limit_X_Kir_Left = meanXvelocOverallKirLeft + stderror_left;
% y_lower_limit_X_Kir_Left = meanXvelocOverallKirLeft - stderror_left;
% y_upper_limit_X_Kir_Right = meanXvelocOverallKirRight + stderror_right;
% y_lower_limit_X_Kir_Right = meanXvelocOverallKirRight - stderror_right;

%Add SD to y upper and limit and subtract from lower limit
y_upper_limit_X_Kir_Left = meanXvelocOverallKirLeft + stdev_left;
y_lower_limit_X_Kir_Left = meanXvelocOverallKirLeft - stdev_left;
y_upper_limit_X_Kir_Right = meanXvelocOverallKirRight + stdev_right;
y_lower_limit_X_Kir_Right = meanXvelocOverallKirRight - stdev_right;
%plot Left antenna x velocity with SEM as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_X_Kir_Left fliplr(y_lower_limit_X_Kir_Left)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin, meanXvelocOverallKirLeft, 'LineWidth',3)

%plot Left antenna x velocity with SEM as shaded area
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
% calculate SEM of Z velocity for right and left antenna
stdev_left = std(cat(1, meanZvelocCellKirLeft{:}));
stderror_left = stdev_left/sqrt(length(cat(1, meanZvelocCellKirLeft{:})));
stdev_right = std(cat(1, meanZvelocCellKirRight{:}));
stderror_right = stdev_right/sqrt(length(cat(1, meanZvelocCellKirRight{:})));

% %Add SEM to y upper and limit and subtract from lower limit
% y_upper_limit_Z_Kir_Left = meanZvelocOverallKirLeft + stderror_left;
% y_lower_limit_Z_Kir_Left = meanZvelocOverallKirLeft - stderror_left;
% y_upper_limit_Z_Kir_Right = meanZvelocOverallKirRight + stderror_right;
% y_lower_limit_Z_Kir_Right = meanZvelocOverallKirRight - stderror_right;

%Add SD to y upper and limit and subtract from lower limit
y_upper_limit_Z_Kir_Left = meanZvelocOverallKirLeft + stdev_left;
y_lower_limit_Z_Kir_Left = meanZvelocOverallKirLeft - stdev_left;
y_upper_limit_Z_Kir_Right = meanZvelocOverallKirRight + stdev_right;
y_lower_limit_Z_Kir_Right = meanZvelocOverallKirRight - stdev_right;

nexttile
%plot Left antenna x velocity with SEM as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Kir_Left fliplr(y_lower_limit_Z_Kir_Left)],...
    [1 0 0], 'FaceAlpha', 0.3,'linestyle', 'none')
hold all
plot(xbin, meanZvelocOverallKirLeft, 'LineWidth',3)

%plot Left antenna x velocity with SEM as shaded area
hold on
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Kir_Right fliplr(y_lower_limit_Z_Kir_Right)],...
    [0 1 0], 'FaceAlpha', 0.2,'linestyle', 'none')
plot(xbin, meanZvelocOverallKirRight, 'LineWidth',3)
title('Z velocity')
rectangle('Position',[1 -2.5 2 5], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylim([-1 2])
ylabel('Z Velocity (mm))');
xlabel('Time (s)');

% %save the figure as a high-quality PNG image
%%%print(gcf,'poke2.png','-dpng','-r300');

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

% %%print(['silenced_overview' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['silenced_overview' '.png'], '-dpng','-r300', '-vector')
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
% ylim([-1 1])

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
% ylim([-1 1])

% %%print(['control_vs_silence_overview' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['control_vs_silence_overview' '.png'], '-dpng','-r300', '-vector')

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

    % %%save median of right and left x velocity per fly
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

    % %%save median of right and left x velocity per fly
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

dataset_rigth = right_antenna_control;
dataset_left = left_antenna_control;

for i = 1 : length(dataset_rigth)
    for k = 1 : numPokes


        % if to_be_removed_ID{1}(i,k) == 1
        %     xVelocityPrePokeControlRight{i}(k) = NaN;
        %     xVelocityDuringPokeControlRight{i}(k) = NaN;
        %     xVelocityPostPokeControlRight{i}(k) = NaN;
        %     continue
        % els

        windowOnRight = dataset_rigth(i).trigger_on(k);
        windowOffRight = dataset_rigth(i).trigger_off(k);

        %Sum spikes 1 s pre poke
        xVelocityPrePokeControlRight{i}(k) = mean(dataset_rigth(i).xveloc_in_mm...
            (windowOnRight-20000:windowOnRight),'omitnan');
        %Sum spikes 1st second of poke
        xVelocityDuringPokeControlRight{i}(k) = mean(dataset_rigth(i).xveloc_in_mm...
            (windowOnRight:windowOnRight+analysis_window_size*sampling_rate),'omitnan');
        %Sum spikes after poke
        xVelocityPostPokeControlRight{i}(k)  = mean(dataset_rigth(i).xveloc_in_mm...
            (windowOffRight:windowOffRight+20000),'omitnan');

    end
end

for i = 1 : length(dataset_rigth)
    for k = 1 : numPokes


        % if to_be_removed_ID{1}(i,k) == 1
        %     xVelocityPrePokeControlRight{i}(k) = NaN;
        %     xVelocityDuringPokeControlRight{i}(k) = NaN;
        %     xVelocityPostPokeControlRight{i}(k) = NaN;
        %     continue
        % els

        windowOnRight = dataset_rigth(i).trigger_on(k);
        windowOffRight = dataset_rigth(i).trigger_off(k);
        windowOnLeft = dataset_left(i).trigger_on(k);
        windowOffLeft = dataset_left(i).trigger_off(k);
        %Sum spikes 1 s pre poke
        xVelocityPrePokeControlRight{i}(k) = mean(dataset_rigth(i).xveloc_in_mm...
            (windowOnRight-20000:windowOnRight),'omitnan');
        xVelocityPrePokeControlLeft{i}(k) = mean(dataset_left(i).xveloc_in_mm...
            (windowOnLeft-20000:windowOnLeft),'omitnan');
        %Sum spikes 1st second of poke
        xVelocityDuringPokeControlRight{i}(k) = mean(dataset_rigth(i).xveloc_in_mm...
            (windowOnRight:windowOnRight+analysis_window_size*sampling_rate),'omitnan');
        xVelocityDuringPokeControlLeft{i}(k) = mean(dataset_left(i).xveloc_in_mm...
            (windowOnLeft:windowOnLeft+analysis_window_size*sampling_rate),'omitnan');
        %Sum spikes after poke
        xVelocityPostPokeControlRight{i}(k)  = mean(dataset_rigth(i).xveloc_in_mm...
            (windowOffRight:windowOffRight+20000),'omitnan');
        xVelocityPostPokeControlLeft{i}(k)  = mean(dataset_left(i).xveloc_in_mm...
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
[pWilPreVsPokeRightControl] = signrank(meanPreXRightControl, meanPokeXRightControl)
[pWilPostVsPokeRightControl] = signrank(meanPostXRightControl, meanPokeXRightControl)

[pWilPreVsPokeLeftControl] = signrank(meanPreXLeftControl, meanPokeXLeftControl)
[pWilPostVsPokeLeftControl] = signrank(meanPostXLeftControl, meanPokeXLeftControl)

%Plot means of firing rate, x & z velocity 1 sec before poke, 1st sec of
%poke, 2nd sec of poke and 1 sec post poke
figure('name', 'Mean firing rates pre, during 1st  500 ms and post stimulus -- RIGHT')
t = tiledlayout('flow');
title(t, 'CONTROL')
nexttile
title("RIGHT" + newline + "N =" + num2str(length(dataset_rigth)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXRightControl; meanPokeXRightControl; meanPostXRightControl];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
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

nexttile
title("LEFT" + newline + "N =" + num2str(length(dataset_left)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXLeftControl; meanPokeXLeftControl; meanPostXLeftControl];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
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

% %%print(['X_control_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['X_control_boxplot' '.png'], '-dpng','-r300', '-vector')


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

% %%print(['Z_control_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['Z_control_boxplot' '.png'], '-dpng','-r300', '-vector')

%% Calculate mean X velocitx and create boxplots for SILENCED

clearvars xVelocityPrePokeRight xVelocityPrePokeLeft xVelocityDuringPokeRight...
    xVelocityDuringPokeLeft xVelocityDuringPokeLeft xVelocityPostPokeLeft...
    xVelocityPostPokeRight

%Size of poking analysis window in s
analysis_window_size = 0.5;

dataset_right = right_antenna_kir;
dataset_left = left_antenna_kir;

for i = 1 : length(dataset_right)
    for k = 1 : numPokes



        windowOnRight = dataset_right(i).trigger_on(k);
        windowOffRight = dataset_right(i).trigger_off(k);
        windowOnLeft = dataset_left(i).trigger_on(k);
        windowOffLeft = dataset_left(i).trigger_off(k);
        %Sum spikes 1 s pre poke
        xVelocityPrePokeKirRight{i}(k) = mean(dataset_right(i).xveloc_in_mm...
            (windowOnRight-20000:windowOnRight),'omitnan');
        xVelocityPrePokeKirLeft{i}(k) = mean(dataset_left(i).xveloc_in_mm...
            (windowOnLeft-20000:windowOnLeft),'omitnan');
        %Sum spikes 1st second of poke
        xVelocityDuringPokeKirRight{i}(k) = mean(dataset_right(i).xveloc_in_mm...
            (windowOnRight:windowOnRight+analysis_window_size*sampling_rate),'omitnan');
        xVelocityDuringPokeKirLeft{i}(k) = mean(dataset_left(i).xveloc_in_mm...
            (windowOnLeft:windowOnLeft+analysis_window_size*sampling_rate),'omitnan');
        %Sum spikes after poke
        xVelocityPostPokeKirRight{i}(k)  = mean(dataset_right(i).xveloc_in_mm...
            (windowOffRight:windowOffRight+20000),'omitnan');
        xVelocityPostPokeKirLeft{i}(k)  = mean(dataset_left(i).xveloc_in_mm...
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
title("RIGHT" + newline + "N =" + num2str(length(dataset_right)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXRightKir; meanPokeXRightKir; meanPostXRightKir];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
text(1,1, ['p=' num2str(pWilPreVsPokeXRightKir)])
text(3,1, ['p=' num2str(pWilPostVsPokeXRightKir)])
ylim([-1 1])
ylabel('X Velocity (mm/s)')


nexttile
title("LEFT" + newline + "N =" + num2str(length(dataset_left)))
%Labels for groups
groupNames = ["Pre"; "Poke"; "Post"];
groups = [1; 2; 3];
plot_y = [meanPreXLeftKir; meanPokeXLeftKir; meanPostXLeftKir];
boxplot2(plot_y',groups)
xticks([1,2,3]);
xticklabels([groupNames]);
hold on
groups = [1,2,3];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled')
line(groups, plot_y)
ylim([-1 1])
ylabel('X Velocity (mm/s)')
text(1,1, ['p=' num2str(pWilPreVsPokeXLeftKir)])
text(3,1, ['p=' num2str(pWilPostVsPokeXLeftKir)])

% %%print(['X_silenced_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['X_silenced_boxplot' '.png'], '-dpng','-r300', '-vector')

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
ylabel('Z Velocity (mm/s)')


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
ylabel('Z Velocity (mm/s)')
text(1,1, ['p=' num2str(pWilPreVsPokeZLeftKir)])
text(3,1, ['p=' num2str(pWilPostVsPokeZLeftKir)])

% %%print(['Z_silenced_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['Z_silenced_boxplot' '.png'], '-dpng','-r300', '-vector')


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

% title(t, "NEW Dataset")

% %%print(['X_control_vs_silenced' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['X_control_vs_silenced' '.png'], '-dpng','-r300', '-vector')

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

% title(t, "NEW Dataset")

% %%print(['Z_control_vs_silenced' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['Z_control_vs_silenced' '.png'], '-dpng','-r300', '-vector')

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
b = boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled', 'jitter','on', 'jitterAmount', 0.1)
ylabel('X Velocity (mm/s)')
ylim([-1.5 1.5])
text(1.5,1, ['p=' num2str(pWilControlvsSilencedPoke)])

% save medians and IQR
medians = zeros(size(b.med));  % preallocate

for i = 1:numel(b.med)
    ydata = get(b.med(i), 'YData');  % YData has 2 identical values
    medians(i) = ydata(1);           % or mean(ydata), doesn't matter
end

numBoxes = numel(b.box);
iqr_ranges = zeros(numBoxes, 2);  % each row: [Q1, Q3]

for i = 1:numBoxes
    ydata = get(b.box(i), 'YData');
    Q1 = min(ydata);  % bottom of box
    Q3 = max(ydata);  % top of box
    iqr_ranges(i, :) = [Q1, Q3];
end

% Generate box labels
labels = arrayfun(@(i) sprintf('Box_%d', i), (1:numBoxes)', 'UniformOutput', false);

% Create table
T = table(labels, medians', iqr_ranges(:,1), iqr_ranges(:,2), ...
    'VariableNames', {'Label', 'Median', 'Q1', 'Q3'});

% Write to text file
% writetable(T, 'boxplot_stats.txt', 'Delimiter', '\t');
% writetable(T, 'fwd_veloc_stats.xlsx');



%Pool x vel curve for right and left  plot control vs silenced
meanXCurveControlPooled = [meanXvelocCellControlRight meanXvelocCellControlLeft];
meanXCurveOverallControl = mean(cat(1, meanXCurveControlPooled{:}));

meanXCurveKirPooled = [meanXvelocCellKirRight meanXvelocCellKirLeft];
meanXCurveOverallKir = mean(cat(1, meanXCurveKirPooled{:}));


temp_matrix = [];
temp_matrix = cell2mat(meanXvelocCellControlPooled'); 
mean_Xveloc_Control = mean(temp_matrix, 1);         
std_Xveloc_Control = std(temp_matrix, 0, 1);        
n = size(temp_matrix, 1);                 
sem_Xveloc_Control = std_Xveloc_Control / sqrt(n);  

temp_matrix = [];
temp_matrix = cell2mat(meanXvelocCellKirPooled'); 
mean_Xveloc_Kir = mean(temp_matrix, 1);         
std_Xveloc_Kir = std(temp_matrix, 0, 1);        
n = size(temp_matrix, 1);                 
sem_Xveloc_Kir = std_Xveloc_Kir / sqrt(n);  

%Add SD to y upper and limit and subtract from lower limit
y_upper_limit_X_Control = meanXvelocOverallControlPooled + sem_Xveloc_Control;
y_lower_limit_X_Control = meanXvelocOverallControlPooled - sem_Xveloc_Control;

%Add SD to y upper and limit and subtract from lower limit
y_upper_limit_X_Kir = meanXvelocOverallKirPooled + sem_Xveloc_Kir;
y_lower_limit_X_Kir = meanXvelocOverallKirPooled - sem_Xveloc_Kir;

nexttile
%plot x velocity with SEM as shaded area
fill([xbin fliplr(xbin)], [y_upper_limit_X_Control fliplr(y_lower_limit_X_Control)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin,meanXCurveOverallControl, 'LineWidth', 3)
hold on
fill([xbin fliplr(xbin)], [y_upper_limit_X_Kir fliplr(y_lower_limit_X_Kir)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin,meanXCurveOverallKir, 'LineWidth', 3)
rectangle('Position',[1 -0.5 2 1], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylabel('X Velocity (mm/s)')
xlabel('Time (s)')
legend('control', 'silenced','Location','northwest')


%print(['X_pooled_control_vs_silenced_line' '.eps'], '-depsc2', '-tiff', '-r300')
%print(['X_pooled_control_vs_silenced_line' '.png'], '-dpng','-r300')

%exportgraphics(gcf, 'X_pooled_control_vs_silenced_line.eps', 'ContentType', 'vector', 'Resolution', 300)

%% POOL RIGHT and LEFT antenna Z and compare CONTROL vs SILENCED

%Pool means 1st 500ms of poke for right and left antenna
%-----------------------------------------------------
%---------RIGHT ANTENNA IS MIRRORED -> *(-1)----------
%-----------------------------------------------------

meanZControlPooled = [meanPokeZRightControl*(-1) meanPokeZLeftControl];
meanZKirPooled = [meanPokeZRightKir*(-1) meanPokeZLeftKir];

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
b = boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled', 'jitter','on', 'jitterAmount', 0.1)
ylabel('Z Velocity (mm/s)')
% ylim([-1.5 1.5])
text(1.5,1, ['p=' num2str(pWilControlvsSilencedZPoke)])

% save medians and IQR
medians = zeros(size(b.med));  % preallocate

for i = 1:numel(b.med)
    ydata = get(b.med(i), 'YData');  % YData has 2 identical values
    medians(i) = ydata(1);           % or mean(ydata), doesn't matter
end

numBoxes = numel(b.box);
iqr_ranges = zeros(numBoxes, 2);  % each row: [Q1, Q3]

for i = 1:numBoxes
    ydata = get(b.box(i), 'YData');
    Q1 = min(ydata);  % bottom of box
    Q3 = max(ydata);  % top of box
    iqr_ranges(i, :) = [Q1, Q3];
end

% Generate box labels
labels = arrayfun(@(i) sprintf('Box_%d', i), (1:numBoxes)', 'UniformOutput', false);

% Create table
T = table(labels, medians', iqr_ranges(:,1), iqr_ranges(:,2), ...
    'VariableNames', {'Label', 'Median', 'Q1', 'Q3'});

% Write to text file
writetable(T, 'angular_veloc_stats.xlsx');


% %%print(['Z_pooled_control_vs_silenced_boxplot' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')

%Pool x vel curve for right and left  plot control vs silenced
meanZvelocCellControlRightMirrored = cellfun(@(x) x*(-1),...
    meanZvelocCellControlRight, 'UniformOutput', false);
meanZCurveControlPooledMirrored = [meanZvelocCellControlRightMirrored meanZvelocCellControlLeft];
meanZCurveOverallControlMirrored = mean(cat(1, meanZCurveControlPooledMirrored{:}));

meanZvelocCellKirRightMirrored = cellfun(@(x) x*(-1),...
    meanZvelocCellKirRight, 'UniformOutput', false);
meanZCurveKirPooledMirrored = [meanZvelocCellKirRightMirrored meanZvelocCellKirLeft];
meanZCurveOverallKirMirrored = mean(cat(1, meanZCurveKirPooledMirrored{:}));

% meanZCurveControlPooled = [meanZvelocCellControlRight meanZvelocCellControlLeft];
% meanZCurveOverallControl = mean(cat(1, meanZCurveControlPooled{:}));
% 
% meanZCurveKirPooled = [meanZvelocCellKirRight meanZvelocCellKirLeft];
% meanZCurveOverallKir = mean(cat(1, meanZCurveKirPooled{:}));


temp_matrix = [];
temp_matrix = cell2mat(meanZvelocCellControlPooled'); 
mean_Zveloc_Control = mean(temp_matrix, 1);         
std_Zveloc_Control= std(temp_matrix, 0, 1);        
n = size(temp_matrix, 1);                 
sem_Zveloc_Control = std_Zveloc_Control / sqrt(n);  

temp_matrix = [];
temp_matrix = cell2mat(meanZvelocCellKirPooled'); 
mean_Zveloc_Kir = mean(temp_matrix, 1);         
std_Zveloc_Kir = std(temp_matrix, 0, 1);        
n = size(temp_matrix, 1);                 
sem_Zveloc_Kir = std_Zveloc_Kir / sqrt(n);  

%Add SD to y upper and limit and subtract from lower limit
y_upper_limit_Z_Control = meanZvelocOverallControlPooled + sem_Zveloc_Control;
y_lower_limit_Z_Control = meanZvelocOverallControlPooled - sem_Zveloc_Control;

%Add SD to y upper and limit and subtract from lower limit
y_upper_limit_Z_Kir = meanZvelocOverallKirPooled + sem_Zveloc_Kir;
y_lower_limit_Z_Kir = meanZvelocOverallKirPooled - sem_Zveloc_Kir;

nexttile
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Control fliplr(y_lower_limit_Z_Control)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin,meanZCurveOverallControlMirrored, 'LineWidth', 3)
hold on
fill([xbin fliplr(xbin)], [y_upper_limit_Z_Kir fliplr(y_lower_limit_Z_Kir)],...
    [1 0 0], 'FaceAlpha', 0.2, 'linestyle', 'none')
hold all
plot(xbin,meanZCurveOverallKirMirrored, 'LineWidth', 3)
rectangle('Position',[1 -0.5 2 1], 'FaceColor', [0.7 0.7 0.7 0.3], 'EdgeColor', 'w')
ylabel('Z Velocity (mm/s)')
xlabel('Time (s)')
legend('control', 'silenced','Location','northwest')


% %%print(['Z_pooled_control_vs_silenced_line' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['Z_pooled_control_vs_silenced_line' '.png'], '-dpng','-r300', '-vector')
%exportgraphics(gcf, 'Z_pooled_control_vs_silenced_line.eps', 'ContentType', 'vector', 'Resolution', 300)

%% Sort for walking / no-walking before poking in CONTROL

walkingThreshold = [10, 0.5, 0.01];

for n = 1 : length(walkingThreshold)
 
    
clearvars medianXVelocitybefore XVelocitybefore meanXVelocitybefore...
    logicalMask meanZVelocitybefore meanYVelocitybefore total_activity...
    XVelocitybefore ZVelocitybefore YVelocitybefore

%Set x  velocity threshold and window size
windowBeforeInSec = 1;
%Walking threshold in mm/s
% walkingThreshold = 0.01;

for i = 1 : length(right_antenna_control)
    for k = 1 : numPokes
        windowOnRight = right_antenna_control(i).trigger_on(k)-...
            (sampling_rate*windowBeforeInSec);
        windowOffRight = right_antenna_control(i).trigger_on(k);
        XVelocitybefore{i,k} = right_antenna_control(i).xveloc_in_mm...
            (windowOnRight:windowOffRight);
        ZVelocitybefore{i,k} = right_antenna_control(i).zveloc_in_mm...
            (windowOnRight:windowOffRight);
%         YVelocitybefore{i,k} = right_antenna_control(i).yveloc_in_mm...
%             (windowOnRight:windowOffRight);

    end
end

meanXVelocitybefore = cellfun(@mean, XVelocitybefore);
meanZVelocitybefore = cellfun(@mean, ZVelocitybefore);
% meanYVelocitybefore = cellfun(@mean, YVelocitybefore);

total_activity = (meanXVelocitybefore);% + abs(meanZVelocitybefore) + abs(meanYVelocitybefore);

%Find index of trials exceeding the treshold
logicalMask = total_activity > walkingThreshold(n);
nRejectedTrialsRightControl = sum(logicalMask, 'all');
nTotalTrialsRightControl = numel(logicalMask);

%Create duplicate of data for filtering
xVelocityDuringPokeControlRightFiltered = xVelocityDuringPokeControlRight;
zVelocityDuringPokeControlRightFiltered = zVelocityDuringPokeControlRight;

%Go through the logical mask and remove values in the X velocity where
%logcialMask == 1 (threshold is exceeded)
for i = 1 : size(logicalMask, 1)
    for k = 1 : size(logicalMask, 2)
        if logicalMask(i,k) == 1
            xVelocityDuringPokeControlRightFiltered{i}(k) = nan;
            zVelocityDuringPokeControlRightFiltered{i}(k) = nan;

        else
        end
    end
end

meanPokeRightControlFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    xVelocityDuringPokeControlRightFiltered);
meanPokeZRightControlFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    zVelocityDuringPokeControlRightFiltered);

clearvars medianXVelocitybefore XVelocitybefore meanXVelocitybefore...
    logicalMask meanZVelocitybefore meanYVelocitybefore total_activity...
    XVelocitybefore ZVelocitybefore YVelocitybefore

for i = 1 : length(left_antenna_control)
    for k = 1 : numPokes
        windowOffLeft = left_antenna_control(i).trigger_on(k);
        windowOnLeft = left_antenna_control(i).trigger_on(k)-...
            (sampling_rate*windowBeforeInSec);
        XVelocitybefore{i,k} = left_antenna_control(i).xveloc_in_mm...
            (windowOnLeft:windowOffLeft);
         ZVelocitybefore{i,k} = left_antenna_control(i).zveloc_in_mm...
            (windowOnRight:windowOffRight);
%         YVelocitybefore{i,k} = left_antenna_control(i).yveloc_in_mm...
%             (windowOnRight:windowOffRight);
    end
end

meanXVelocitybefore = cellfun(@mean, XVelocitybefore);
meanZVelocitybefore = cellfun(@mean, ZVelocitybefore);
% meanYVelocitybefore = cellfun(@mean, YVelocitybefore);

total_activity = (meanXVelocitybefore);% + abs(meanZVelocitybefore) + abs(meanYVelocitybefore);

%Find index of trials exceeding the treshold
logicalMask = total_activity > walkingThreshold(n);
nRejectedTrialsLeftControl = sum(logicalMask, 'all');
nTotalTrialsLeftControl = numel(logicalMask);

%Create duplicate of data for filtering
xVelocityDuringPokeControlLeftFiltered = xVelocityDuringPokeControlLeft;
zVelocityDuringPokeControlLeftFiltered = zVelocityDuringPokeControlLeft;

%Go through the logical mask and remove values in the X velocity where
%logcialMask == 1 (threshold is exceeded)
for i = 1 : size(logicalMask, 1)
    for k = 1 : size(logicalMask, 2)
        if logicalMask(i,k) == 1
            xVelocityDuringPokeControlLeftFiltered{i}(k) = nan;
            zVelocityDuringPokeControlLeftFiltered{i}(k) = nan;
           
        else
        end
    end
end

meanPokeLeftControlFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    xVelocityDuringPokeControlLeftFiltered);
meanPokeZLeftControlFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    zVelocityDuringPokeControlLeftFiltered);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort for no walking before poking in SILENCED


clearvars medianXVelocitybefore XVelocitybefore meanXVelocitybefore...
    logicalMask meanZVelocitybefore meanYVelocitybefore total_activity...
    XVelocitybefore ZVelocitybefore YVelocitybefore

for i = 1 : length(right_antenna_kir)
    for k = 1 : numPokes
        windowOffRight = right_antenna_kir(i).trigger_on(k);
        windowOnRight = right_antenna_kir(i).trigger_on(k)-...
            (sampling_rate*windowBeforeInSec);
        XVelocitybefore{i,k} = right_antenna_kir(i).xveloc_in_mm...
            (windowOnRight:windowOffRight);
                ZVelocitybefore{i,k} = right_antenna_kir(i).zveloc_in_mm...
            (windowOnRight:windowOffRight);
%         YVelocitybefore{i,k} = right_antenna_kir(i).yveloc_in_mm...
%             (windowOnRight:windowOffRight);
    end
end

meanXVelocitybefore = cellfun(@mean, XVelocitybefore);
meanZVelocitybefore = cellfun(@mean, ZVelocitybefore);
% meanYVelocitybefore = cellfun(@mean, YVelocitybefore);

total_activity = (meanXVelocitybefore);% + abs(meanZVelocitybefore) + abs(meanYVelocitybefore);

%Find index of trials exceeding the treshold
logicalMask = total_activity > walkingThreshold(n);
nRejectedTrialsRightKir = sum(logicalMask, 'all');
nTotalTrialsRightKir = numel(logicalMask);

%Create duplicate of data for filtering
xVelocityDuringPokeKirRightFiltered = xVelocityDuringPokeKirRight;
zVelocityDuringPokeKirRightFiltered = zVelocityDuringPokeKirRight;


%Go through the logical mask and remove values in the X velocity where
%logcialMask == 1 (threshold is exceeded)
for i = 1 : size(logicalMask, 1)
    for k = 1 : size(logicalMask, 2)
        if logicalMask(i,k) == 1
            xVelocityDuringPokeKirRightFiltered{i}(k) = nan;
            zVelocityDuringPokeKirRightFiltered{i}(k) = nan;

        else
        end
    end
end

meanPokeRightKirFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    xVelocityDuringPokeKirRightFiltered);
meanPokeZRightKirFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    zVelocityDuringPokeKirRightFiltered);

clearvars medianXVelocitybefore XVelocitybefore meanXVelocitybefore...
    logicalMask meanZVelocitybefore meanYVelocitybefore total_activity...
    XVelocitybefore ZVelocitybefore YVelocitybefore

for i = 1 : length(left_antenna_kir)
    for k = 1 : numPokes
        windowOffLeft = left_antenna_kir(i).trigger_on(k);
        windowOnLeft = left_antenna_kir(i).trigger_on(k)-...
            (sampling_rate*windowBeforeInSec);
        XVelocitybefore{i,k} = left_antenna_kir(i).xveloc_in_mm...
            (windowOnLeft:windowOffLeft);
        ZVelocitybefore{i,k} = left_antenna_kir(i).zveloc_in_mm...
            (windowOnRight:windowOffRight);
%         YVelocitybefore{i,k} = left_antenna_kir(i).yveloc_in_mm...
%             (windowOnRight:windowOffRight);
    end
end

meanXVelocitybefore = cellfun(@mean, XVelocitybefore);
meanZVelocitybefore = cellfun(@mean, ZVelocitybefore);
% meanYVelocitybefore = cellfun(@mean, YVelocitybefore);

total_activity = (meanXVelocitybefore);% + abs(meanZVelocitybefore) + abs(meanYVelocitybefore);

%Find index of trials exceeding the treshold
logicalMask = total_activity > walkingThreshold(n);
nRejectedTrialsLeftKir = sum(logicalMask, 'all');
nTotalTrialsLeftKir = numel(logicalMask);

%Create duplicate of data for filtering
xVelocityDuringPokeKirLeftFiltered = xVelocityDuringPokeKirLeft;
zVelocityDuringPokeKirLeftFiltered = zVelocityDuringPokeKirLeft;

%Go through the logical mask and remove values in the X velocity where
%logcialMask == 1 (threshold is exceeded)
for i = 1 : size(logicalMask, 1)
    for k = 1 : size(logicalMask, 2)
        if logicalMask(i,k) == 1
            xVelocityDuringPokeKirLeftFiltered{i}(k) = nan;
            zVelocityDuringPokeKirLeftFiltered{i}(k) = nan;

        else
        end
    end
end

meanPokeLeftKirFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    xVelocityDuringPokeKirLeftFiltered);
meanPokeZLeftKirFiltered = cellfun(@(x) mean(x, 'all', 'omitnan'),...
    zVelocityDuringPokeKirLeftFiltered);

% Statistics test SILENCED vs CONTROL with x velocity threshold fitered data!

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanPokeRightControlFilteredPlot, meanPokeRightKirFilteredPlot ] = ...
    padarrays(meanPokeRightControlFiltered, meanPokeRightKirFiltered);
[meanPokeLeftControlFilteredPlot, meanPokeLeftKirFilteredPlot ] = ...
    padarrays(meanPokeLeftControlFiltered, meanPokeLeftKirFiltered);

[pWilControlvsSilencedPokeRightFiltered] = ranksum(meanPokeRightControlFiltered,...
    meanPokeRightKirFilteredPlot);
[pWilControlvsSilencedPokeLeftFiltered] = ranksum(meanPokeLeftControlFiltered,...
    meanPokeLeftKirFilteredPlot);

figure
t = tiledlayout('flow');
nexttile;
title("Right - walking before < " + num2str(walkingThreshold(n)) + "mm/s")
groupNames = ["Control "+'n = ' + num2str(nTotalTrialsRightControl-nRejectedTrialsRightControl)...
    + "/" + num2str(nTotalTrialsRightControl); "Silenced "+'n = ' + ...
    num2str(nTotalTrialsRightKir-nRejectedTrialsRightKir) + "/" + num2str(nTotalTrialsRightKir)];
groups = [1; 2];
plot_y = [meanPokeRightControlFilteredPlot; meanPokeRightKirFilteredPlot];

boxplot2(plot_y', groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled', 'jitter','on', 'jitterAmount', 0.1)
ylabel('X Velocity (mm/s)')
ylim([-2 2])
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeRightFiltered)])


nexttile
title("Left - walking before < " + num2str(walkingThreshold(n)) + "mm/s")
groupNames = ["Control "+'n = ' + num2str(nTotalTrialsLeftControl-nRejectedTrialsLeftControl) +...
    "/" + num2str(nTotalTrialsLeftControl); "Silenced "+'n = ' + ...
    num2str(nTotalTrialsLeftKir-nRejectedTrialsLeftKir) + "/" + num2str(nTotalTrialsLeftKir)];
groups = [1; 2];
plot_y = [meanPokeLeftControlFilteredPlot; meanPokeLeftKirFilteredPlot];
boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled', 'jitter','on', 'jitterAmount', 0.1)
ylabel('X Velocity (mm/s)')
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeLeftFiltered)])
ylim([-2 2])

% title(t, "NEW Dataset")

% %%print(['no_walking_filtered' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['no_walking_filtered' '.png'], '-dpng','-r300', '-vector')



%%% POOL RIGHT and LEFT antenna X and compare CONTROL vs SILENCED

%Pool means 1st 500ms of poke for right and left antenna
meanXControlPooled = [meanPokeRightControlFiltered meanPokeLeftControlFiltered];
meanXKirPooled = [meanPokeRightKirFiltered meanPokeLeftKirFiltered];

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanControlPooledPlot, meanKirPooledPlot ] = padarrays(meanXControlPooled, meanXKirPooled);

%Ranksum test for control vs silenced pooled 1st 500ms of poke
[pWilControlvsSilencedPoke] = ranksum(meanXControlPooled, meanXKirPooled);

%Plot means in a boxplot, control vs silenced

nexttile
title("Pooled - walking before < " + num2str(walkingThreshold(n)) + "mm/s")
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
ylim([-2 2])
text(1.5,1, ['p=' num2str(pWilControlvsSilencedPoke)])

set(gcf,'position',[400, 100, 400,1000])
% %%print(['Walking_filtered' '.eps'], '-depsc2', '-tiff', '-r300')
% %%print(['Walking_filtered' '.png'], '-dpng','-r300')

% Statistics test SILENCED vs CONTROL with x velocity threshold fitered data!

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanPokeZRightControlFilteredPlot, meanPokeZRightKirFilteredPlot ] = ...
    padarrays(meanPokeZRightControlFiltered, meanPokeZRightKirFiltered);
[meanPokeZLeftControlFilteredPlot, meanPokeZLeftKirFilteredPlot ] = ...
    padarrays(meanPokeZLeftControlFiltered, meanPokeZLeftKirFiltered);

[pWilControlvsSilencedPokeRightFiltered] = ranksum(meanPokeZRightControlFiltered,...
    meanPokeZRightKirFilteredPlot);
[pWilControlvsSilencedPokeLeftFiltered] = ranksum(meanPokeZLeftControlFiltered,...
    meanPokeZLeftKirFilteredPlot);


filenameEps = "x veloc threshold " + num2str(walkingThreshold(n)) + '.eps';
filenamePng = "x veloc threshold " + num2str(walkingThreshold(n)) + '.png';

%%print([filenameEps], '-depsc2', '-tiff', '-r300')
%%print([filenamePng], '-dpng','-r300')

figure
t = tiledlayout('flow');
nexttile;
title("Right - walking before < " + num2str(walkingThreshold(n)) + "mm/s")
groupNames = ["Control "+'n = ' + num2str(nTotalTrialsRightControl-nRejectedTrialsRightControl)...
    + "/" + num2str(nTotalTrialsRightControl); "Silenced "+'n = ' + ...
    num2str(nTotalTrialsRightKir-nRejectedTrialsRightKir) + "/" + num2str(nTotalTrialsRightKir)];
groups = [1; 2];
plot_y = [meanPokeZRightControlFilteredPlot; meanPokeZRightKirFilteredPlot];
boxplot2(plot_y', groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled', 'jitter','on', 'jitterAmount', 0.1)
ylabel('Z Velocity (mm/s)')
ylim([-7 7])
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeRightFiltered)])


nexttile
title("Left - walking before < " + num2str(walkingThreshold(n)) + "mm/s")
groupNames = ["Control "+'n = ' + num2str(nTotalTrialsLeftControl-nRejectedTrialsLeftControl) +...
    "/" + num2str(nTotalTrialsLeftControl); "Silenced "+'n = ' + ...
    num2str(nTotalTrialsLeftKir-nRejectedTrialsLeftKir) + "/" + num2str(nTotalTrialsLeftKir)];
groups = [1; 2];
plot_y = [meanPokeZLeftControlFilteredPlot; meanPokeZLeftKirFilteredPlot];
boxplot2(plot_y',groups)
xticks([1,2]);
xticklabels([groupNames]);
hold on
groups = [1,2];
%Use scatter to plot filled circles of each mean
scatter(groups, plot_y, 'filled', 'jitter','on', 'jitterAmount', 0.1)
ylabel('Z Velocity (mm/s)')
text(1,1, ['p=' num2str(pWilControlvsSilencedPokeLeftFiltered)])
ylim([-7 7])

% title(t, "NEW Dataset")

% %%print(['no_walking_filtered' '.eps'], '-depsc2', '-tiff', '-r300', '-vector')
% %%print(['no_walking_filtered' '.png'], '-dpng','-r300', '-vector')



%%% POOL RIGHT and LEFT antenna X and compare CONTROL vs SILENCED
meanZControlPooled = [];
meanZKirPooled = [];

%Pool means 1st 500ms of poke for right and left antenna
meanZControlPooled = [meanPokeZRightControlFiltered*(-1) meanPokeZLeftControlFiltered];
meanZKirPooled = [meanPokeZRightKirFiltered*(-1) meanPokeZLeftKirFiltered];

%Fill up arrays of different lengths with NANs to the same size for
%PLOTTING ONLY, function: padarrays
[meanControlPooledPlot, meanKirPooledPlot ] = padarrays(meanZControlPooled, meanZKirPooled);

%Ranksum test for control vs silenced pooled 1st 500ms of poke
[pWilControlvsSilencedPoke] = ranksum(meanZControlPooled, meanZKirPooled);

%Plot means in a boxplot, control vs silenced

nexttile
title("Pooled - walking before < " + num2str(walkingThreshold(n)) + "mm/s")
groupNames = ["Control " + "N = " + num2str(length(meanZControlPooled));...
    "Silenced " + "N = "+ num2str(length(meanZKirPooled))];
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
ylim([-7 7])
text(1.5,1, ['p=' num2str(pWilControlvsSilencedPoke)])

set(gcf,'position',[400, 100, 400,1000])

filenameEps = "z veloc threshold " + num2str(walkingThreshold(n)) + '.eps';
filenamePng = "z veloc threshold " + num2str(walkingThreshold(n)) + '.png';

%%print([filenameEps], '-depsc2', '-tiff', '-r300')
%%print([filenamePng], '-dpng','-r300')

end
%% Compare overall acitvity between CONTROL and SILENCED

ActivityControlRight = ([right_antenna_control.mean_velocity]);
ActivityControlLeft = ([left_antenna_control.mean_velocity]);
ActivityKirRight = ([right_antenna_kir.mean_velocity]);
ActivityKirLeft = ([left_antenna_kir.mean_velocity]);

[ActivityControlRightPlot, ActivityKirRightPlot] = padarrays(ActivityControlRight, ActivityKirRight);
[ActivityControlLeftPlot, ActivityKirLeftPlot] = padarrays(ActivityControlLeft, ActivityKirLeft);

[pWilControlvsSilencedActivityRight] = ranksum(ActivityControlRight,...
    ActivityKirRight);
[pWilControlvsSilencedActivityLeft] = ranksum(ActivityControlLeft,...
    ActivityKirLeft);

figure
title("Mean velocity Control vs Silenced")
plotY = [ActivityControlRightPlot; ActivityKirRightPlot;... 
    ActivityControlLeftPlot; ActivityKirLeftPlot];
groups = [1; 2; 3; 4];
groupNames = ["Control Right"; "Silenced Right"; "Control Left"; "Silenced Left"]
boxplot2(plotY', groups)
xticks([1,2,3,4]);
xticklabels([groupNames]);
ylabel("Mean velocity (mm/s)")
text(1,1, ['p=' num2str(pWilControlvsSilencedActivityRight)])
text(3,1, ['p=' num2str(pWilControlvsSilencedActivityLeft)])

