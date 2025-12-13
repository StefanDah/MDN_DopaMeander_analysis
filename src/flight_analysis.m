%==========================================================================
% File Name: flight_analysis.m
%
% Description:
%   Analyse flight bouts detected during preprocessing
%
% Original Author:
%   Sander Ließem (Neurobiology and Genetics, Julius-Maximilians-University
%   of Würzburg)
%
% Notes:
%   Used and modified with permission.

%% Cleaning up and defining workspace
clear all
close all

% Find current path and load data
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);  % only works if file is saved
filepath = fileparts(scriptPath);

%% Introducing Global Variables, Path and predifined Paramteres

cd(filepath)
addpath(genpath(filepath));

files = dir('20*.mat'); %Search for all files
flies = [1:length(files)];
analysis_file = dir('Analysis_TACHO.mat');

% % % % % sampling_rate = 20000;
% ftrigwin=5;                         % size of pre/post window in seconds
binsize_factor=0.250;               % binsize=binsize_factor*sampling_rate
min_flight_freq = 300;              % when set to 150: is 133Hz at 20.000 Sampling Rate, 300 is 66 Hz
min_flight_duration_in_sec = 4;
artifact_during_flight = 1500;      % loss of signal during one cont. flight bout tollerance: 1000 is 50 ms at 20000 sampling rate
min_interflight_interval = 200000;  % only takes flight bouts into account in which 15s (when set to 300000 at sampling rate 20000) no flight is detected for and after one bout
% % % % % analysis_window_pre_onset = 5*sampling_rate; % in seconds*sampling rate Bin Window arround flight onset
% % % % % analysis_window_pre_offset = 15*sampling_rate; % Bin Window arround flight offset
checkbox_show_all_results_onset = 1;      % Shows onset of all Flies analyzed
checkbox_manual_evaluation = 0;     % If set to zero, only the final result graphs will be presented, no corrections possible
input_time_width_average_window = 8;   %later used for plotting graph, if 4 by binsize of 0.5 => 2 seconds
input_current_amplitude_injection = 4.5; %threshold, used to evaluate whether current has been injected. If True, results will be saved seperatly
minimum_spike_amplitude = 5;       %Unless current is injected, data will not be used in analysis

idx_excluded_flight_bouts = 0;      % Indicator how many trials have been discarded due to interflight interval

allspikes_allflies=[];
intra_binned=[];
flight_binned=[];
flight_onset_trig=[];
idx_fly_counter = 0;% Used later to assign each fly in averaged data
idx_fly_counter_clean = 0;% Used later to assign each fly in averaged data
file_modifier = 'clean';

checkbox_flyID = 1;


%Analysis Paramter
analysis_membrane_potential_prior_onset_in_sec = 1;
analysis_membrane_potential_post_onset_in_sec = 1;
analysis_membrane_potential_prior_offset_in_sec = 1;
analysis_membrane_potential_post_offset_in_sec = 1;

analysis_window_pre_onset_in_sec = 2;
analysis_window_post_offset_in_sec = 2;
analysis_window_post_onset_in_sec = 2;
analysis_window_pre_offset_in_sec = 2;

smoothing_factor_tacho = 5;

dlgTitle    = 'Overwrite Exsiting Analysis Winows of Mat Files (Preprocessing)';
dlgQuestion = 'Overwrite Exsiting Analysis Winows of Mat Files (Preprocessing)? New values must be equal or shorter!';
choice = questdlg(dlgQuestion,dlgTitle,'YES','NO', 'NO');
if length(choice)==2
    overwrite_analysis_windows = 0;
else
    overwrite_analysis_windows = 1;
end

% Figure Paramters
dlgTitle    = 'Plott Figures in Black enviroment?';%NEW 2021_08_13
dlgQuestion = 'Plott Figures in Black enviroment?';
save_choice = questdlg(dlgQuestion,dlgTitle,'YES','NO', 'NO');
if length(save_choice)==3  %clean
    plotting_figures_in_black = 1;
    figure_plotting_color = [0 0 0];
else
    plotting_figures_in_black = 0;
    figure_plotting_color = [1 1 1];
end



%Check if Analysis has already been conducted
if isempty(analysis_file)==0
    
    dlgTitle    = 'Start over?';
    dlgQuestion = 'Analysis has already been made. Start over?';
    choice = questdlg(dlgQuestion,dlgTitle,'YES','NO', 'YES');
    if length(choice)==2
        analysis_already_performed = 1;
        eval(['load ' analysis_file.name]);
    else
        analysis_already_performed = 0;
    end
    sampling_rate = 20000;
    binsize=binsize_factor*sampling_rate;
else
    analysis_already_performed = 0;
end

% Save individual Plots?
dlgTitle    = 'Save upcoming Result Figures?';
dlgQuestion = 'Do you wish to save all upcoming result figures?';
save_choice = questdlg(dlgQuestion,dlgTitle,'YES','NO', 'NO');
if length(save_choice)==3  %save all
    input_save_figures = 1;
else
    input_save_figures = 0;
end

%clrs=generate_colors_insulin_flight(1, 2, 3, 4, 5, 6,7);
%flies=length(files);

% % % % % if analysis_already_performed == 0
%% ************* Start of Analysis on Files *************

for fly=flies%length(files);
    clearvars intra_binned flight_binned binstarts binends flight_onset_trig%Clearing Workspace
    eval(['load ' files(fly).name ]);       %Load Data from Mat file for particular fly
    
    %IF FILE ID HAST NOT BEEN DETERMINED SO FAR:
    if checkbox_flyID ==1
    if strcmp(num2str(flyID), 'FLY_CELL_TRIAL') == 1
         %Read Fly Data from excel sheet (FlyID)
    [FileID, FlyIDENT, CellID, Trial] = fcn_Excelimport("X:\MATLAB\data_flight\paper_2025\Fly_ID.xlsx", "Sheet1", [1, 111]); %last input is the number of flies in Excel sheet
    
    %Shorten These Variables in case empty values need to be excluded
    FileID(find(isnan(Trial)==1)) = [];
    FlyIDENT(find(isnan(Trial)==1)) = [];
    CellID(find(isnan(Trial)==1)) = [];
    Trial(find(isnan(Trial)==1)) = [];
    file = files(fly).name;
        idx_flyID = find(FileID==file(1:15));%was end-4
                current_File = FileID(idx_flyID);
                current_FlyNumber = FlyIDENT(idx_flyID);
                current_CellNumber = CellID(idx_flyID);
                current_Trial = Trial(idx_flyID);
                
                flyID = strcat(['FLY', num2str(current_FlyNumber), '_CELL', num2str(current_CellNumber), '_TRIAL', num2str(current_Trial)]);
%                 disp(['Current FlyID is:' flyID])
    end
    end
    
    disp(['Current Fly is: ' num2str(flyID) ' from MAT-File: ' num2str(current_File)])
    
    %Hardcoded
    binsize=binsize_factor*sampling_rate;
    smoothing_factor_membrane = sampling_rate/5; %Used to smooth intra trace and velocity trace for plotting
    
    if overwrite_analysis_windows == 1
        analysis_window_pre_onset = analysis_window_pre_onset_in_sec*sampling_rate;
        analysis_window_post_offset = analysis_window_post_offset_in_sec*sampling_rate;
        analysis_window_post_onset = analysis_window_post_onset_in_sec*sampling_rate;
        analysis_window_pre_offset = analysis_window_pre_offset_in_sec*sampling_rate;
    end
    
    
    % Giving Variables proper names and cleanup workspace
    flight_onset_trigger = flight_onset;
    flight_offset_trigger = flight_offset;
    intra_smoothed = intra; %The original intra is saved as intra_orig
    clearvars flight_onset flight_offset intra
    
    %In Case no Flight bouts can be used
    if isempty(flight_onset_trigger)==1 || isempty(flight_offset_trigger)==1
        disp('Trial excluded because no Flight Triggers detected')
        continue
    end
    
    
    
    %% Figure Creation for Every Trial
    idx_figure_counter = 1; %Counter used to count amount of figures created (for figure name)
    for i = 1:length(flight_onset_trigger)
        result_figure = figure;
        set(result_figure, 'position', [1, 200, 1900, 800], 'Color', figure_plotting_color);
        subplot(2,1,1)
        hold on
        %plotting intra
        x = (flight_onset_trigger(i)-analysis_window_pre_onset:flight_offset_trigger(i)+analysis_window_post_offset)/sampling_rate;
        plot(x, intra_smoothed(flight_onset_trigger(i)-analysis_window_pre_onset:flight_offset_trigger(i)+analysis_window_post_offset))
        
        %setting axes and enviroment of figur %2021_08_13
        if plotting_figures_in_black == 1
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',12);
            ax_membranepotential_figure = gca;
            set(gca,'color',[0 0 0]);
            set(gcf,'inverthardcopy','off');
            set(ax_membranepotential_figure,'YColor','w', 'XColor','w')
            set(gca,'TickDir','out');
            set(gca, 'Xcolor', 'k');
            set(gcf,'inverthardcopy','off');
            set(ax_membranepotential_figure,'YColor','w', 'XColor','w')
            set(gca,'TickDir','out');
        else
            ax_membranepotential_figure = gca;
            set(gca,'color',[1 1 1]);
            set(gcf,'inverthardcopy','off');
            set(ax_membranepotential_figure,'YColor','k', 'XColor','k')
            set(gca,'TickDir','out');
            set(gca, 'Xcolor', 'k');
            set(gcf,'inverthardcopy','off');
            set(ax_membranepotential_figure,'YColor','k', 'XColor','k')
            set(gca,'TickDir','out');
        end
        ylabel('Vm (mV)');
        
        
        subplot(2,1,2)
        hold on
        %plotting activity
        x = (flight_onset_trigger(i)-analysis_window_pre_onset:flight_offset_trigger(i)+analysis_window_post_offset)/sampling_rate;
        plot(x, tacho(flight_onset_trigger(i)-analysis_window_pre_onset:flight_offset_trigger(i)+analysis_window_post_offset))
        %plotting Flight bout
        x = (flight_onset_trigger(i):flight_offset_trigger(i))/sampling_rate;
        p = patch([(flight_onset_trigger(i))/sampling_rate (flight_onset_trigger(i))/sampling_rate, (flight_offset_trigger(i))/sampling_rate (flight_offset_trigger(i))/sampling_rate], [min(ylim) max(ylim) max(ylim) min(ylim)], [1 0 1], 'EdgeColor', 'none');
        alpha(0.3)
        %plotting Analysis window (gray)
        p = patch([(flight_onset_trigger(i)-analysis_window_pre_onset)/sampling_rate (flight_onset_trigger(i)-analysis_window_pre_onset)/sampling_rate, (flight_onset_trigger(i))/sampling_rate (flight_onset_trigger(i))/sampling_rate], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.8 0.8 0.8], 'EdgeColor', 'none');
        alpha(0.3)
        p = patch([(flight_offset_trigger(i))/sampling_rate (flight_offset_trigger(i))/sampling_rate, (flight_offset_trigger(i)+analysis_window_post_offset)/sampling_rate (flight_offset_trigger(i)+analysis_window_post_offset)/sampling_rate], [min(ylim) max(ylim) max(ylim) min(ylim)], [0.8 0.8 0.8], 'EdgeColor', 'none');
        alpha(0.3)
        x = [flight_onset_trigger(i)-analysis_window_pre_onset flight_onset_trigger(i)-analysis_window_pre_onset]/sampling_rate;
        plot(x, [min(ylim) max(ylim)], '--', 'Color', 'k')
        x = [flight_onset_trigger(i) flight_onset_trigger(i)]/sampling_rate;
        plot(x, [min(ylim) max(ylim)], '--', 'Color', 'k')
        x = [flight_offset_trigger(i)+analysis_window_post_offset flight_offset_trigger(i)+analysis_window_post_offset]/sampling_rate;
        plot(x, [min(ylim) max(ylim)], '--', 'Color', 'k')
        x = [flight_offset_trigger(i) flight_offset_trigger(i)]/sampling_rate;
        plot(x, [min(ylim) max(ylim)], '--', 'Color', 'k')
        %Plotting Analysis window wihtin bout
        x = [flight_onset_trigger(i)+analysis_window_post_onset flight_onset_trigger(i)+analysis_window_post_onset]/sampling_rate;
        plot(x, [min(ylim) max(ylim)], '--', 'Color', 'k')
        x = [flight_offset_trigger(i)-analysis_window_pre_offset flight_offset_trigger(i)-analysis_window_pre_offset]/sampling_rate;
        plot(x, [min(ylim) max(ylim)], '--', 'Color', 'k')
        
        ylabel('Fly Activity'); %2021_08_13
        
        %setting axes and enviroment of figur %2021_08_13
        if plotting_figures_in_black == 1
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',12);
            ax_membranepotential_figure = gca;
            set(gca,'color',[0 0 0]);
            set(gcf,'inverthardcopy','off');
            set(ax_membranepotential_figure,'YColor','w', 'XColor','w')
            set(gca,'TickDir','out');
            set(gca, 'Xcolor', 'k');
            set(gcf,'inverthardcopy','off');
            set(ax_membranepotential_figure,'YColor','w', 'XColor','w')
            set(gca,'TickDir','out');
        else
            ax_membranepotential_figure = gca;
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',12);
            set(gca,'color',[1 1 1]);
            set(gcf,'inverthardcopy','off');
            set(ax_membranepotential_figure,'YColor','k', 'XColor','k')
            set(gca,'TickDir','out');
            set(gca, 'Xcolor', 'k');
            set(gcf,'inverthardcopy','off');
            set(ax_membranepotential_figure,'YColor','k', 'XColor','k')
            set(gca,'TickDir','out');
        end
        
        
        
        if input_save_figures == 1
            filename = [files(fly).name(1:end-4) '_clean_'  num2str(idx_figure_counter) '.tif'];
            saveas(gcf, filename)
            filename = [files(fly).name(1:end-4) '_clean_'  num2str(idx_figure_counter) '.eps'];
            print(filename, '-depsc2', '-tiff', '-r300', '-painters')
            pause(1)
            close gcf
        end
        close gcf
        idx_figure_counter = idx_figure_counter + 1;
    end
    
    
    
    %% ************* MAIN RESULT STRUCT: Writing intra data for flight on and offset windows into a cell, together with bin information and amount of spikes *************
    %only if flight triggers were detected
    if isempty(flight_onset_trigger)==0
        
        
        
        
        binsize=binsize_factor*sampling_rate;
        % Variable were everythin is saved in: "Analysis_Data_Flight. Next step,
        % determine wheter this variable already contained data from previouse
        % flys
        if exist('Analysis_Data_Flight') == 1
            bout_counter = length(Analysis_Data_Flight(:,1));
        else
            bout_counter = 0;
        end
        
        
        %In case flight_offset time is to long: Intra trace is
        %filled up with the last value of intra
        if length(flight_onset_trigger) ~= 0 && flight_offset_trigger(end)+analysis_window_pre_offset > length(intra_smoothed)
            flight_onset_trigger(end) = [];
            flight_offset_trigger(end)= [];
        end
        
        %In case flight_offset time is to long: Intra trace is
        %filled up with the last value of intra
        if length(flight_onset_trigger) ~= 0 && flight_onset_trigger(1)-analysis_window_pre_onset < 0;
            flight_onset_trigger(1) = [];
            flight_offset_trigger(1)= [];
        end
        
        for i = bout_counter + 1:bout_counter + length(flight_onset_trigger)
            
            %first general info about the fly, cell and bout
            Analysis_Data_Flight{i,1} = flyID;
            Analysis_Data_Flight{i,2} = files(fly).name(1:end-4);
            Analysis_Data_Flight{i,3} = current;
            
            %Intra Trace shortly before onset of walking
            if flight_onset_trigger(i-bout_counter)-analysis_window_pre_onset<0 %In case fly is flying in the beginning
                analysis_window_pre_onset_back = analysis_window_pre_onset;
                analysis_window_pre_onset = 0;
                disp(['One flying Bout at ' num2str(flight_onset_trigger(i-bout_counter)) 'was excluded because it was at the beginning of the walking'])
            else analysis_window_pre_onset_back = analysis_window_pre_onset;
            end
            Analysis_Data_Flight{i,4} = intra_smoothed(flight_onset_trigger(i-bout_counter)-analysis_window_pre_onset:flight_onset_trigger(i-bout_counter)+analysis_window_post_onset-1);
            
            %NEW: Spikes within this time window
            Analysis_Data_Flight{i,12} = spikes(flight_onset_trigger(i-bout_counter)-analysis_window_pre_onset+binsize:flight_onset_trigger(i-bout_counter)+analysis_window_post_onset+binsize); %-binsize to have the right amount of bins
            %**
            
            %binning starts
            Analysis_Data_Flight{i,5} = flight_onset_trigger(i-bout_counter)-analysis_window_pre_onset-binsize:binsize:flight_onset_trigger(i-bout_counter)+analysis_window_post_onset;
            
            %binning ends
            Analysis_Data_Flight{i,6} = Analysis_Data_Flight{i,5}+binsize;
            
            %spikes in bin
            for z=1:length(Analysis_Data_Flight{i,5})
                %calculates number of spikes wihtin each bin
                Analysis_Data_Flight{i,7}(z) = nansum(spikes(Analysis_Data_Flight{i,5}(z):Analysis_Data_Flight{i,6}(z)));
            end
            
            
            %Intra Trace shortly before and after offset of walking
            Analysis_Data_Flight{i,8} = intra_smoothed(flight_offset_trigger(i-bout_counter)-analysis_window_pre_onset:flight_offset_trigger(i-bout_counter)+analysis_window_post_onset-1);
            
            %Saves the onset (9) and offset (10) trigger for later
            Analysis_Data_Flight{i,9} = flight_onset_trigger(i-bout_counter);
            
            Analysis_Data_Flight{i,10} = flight_offset_trigger(i-bout_counter);
            
            %Exact Spikes times in each flying Bouts (11)
            Analysis_Data_Flight{i,11} = find(spikes(flight_onset_trigger(i-bout_counter)-analysis_window_pre_onset:flight_onset_trigger(i-bout_counter)+analysis_window_post_onset-1)==1);
            
            
            %*************** NEW BINNING ARROUND ONSET:
            %New Try of binning: until exact onset (1):
            Analysis_Data_Flight{i,13} = flight_onset_trigger(i-bout_counter)-analysis_window_pre_onset:binsize:flight_onset_trigger(i-bout_counter);
            
            %New Try of binning: from onset to end of time window (2, not offset):
            Analysis_Data_Flight{i,14} = flight_onset_trigger(i-bout_counter):binsize:flight_onset_trigger(i-bout_counter)+analysis_window_post_onset;
            
            %NEW: Spikes within these time windows
            %** (1) BINNNED
            for z=1:length(Analysis_Data_Flight{i,13})-1
                %calculates number of spikes wihtin each bin
                Analysis_Data_Flight{i,15}(z) = nansum(spikes(Analysis_Data_Flight{i,13}(z):Analysis_Data_Flight{i,13}(z+1)));
            end
            %** (1) UNBINNNED = SPIKEVECTOR
            Analysis_Data_Flight{i,16} = spikes(Analysis_Data_Flight{i,13}(1):Analysis_Data_Flight{i,13}(end)-1);
            Analysis_Data_Flight{i,16}(isnan(Analysis_Data_Flight{i,16}))=0;
            Analysis_Data_Flight{i,16} = Analysis_Data_Flight{i,16}';
            
            %** (2)
            for z=1:length(Analysis_Data_Flight{i,14})-1
                %calculates number of spikes wihtin each bin
                Analysis_Data_Flight{i,17}(z) = nansum(spikes(Analysis_Data_Flight{i,14}(z):Analysis_Data_Flight{i,14}(z+1)));
            end
            %** (2) UNBINNNED = SPIKEVECTOR
            Analysis_Data_Flight{i,18} = spikes(Analysis_Data_Flight{i,14}(1):Analysis_Data_Flight{i,14}(end)-1);
            Analysis_Data_Flight{i,18}(isnan(Analysis_Data_Flight{i,18}))=0;
            Analysis_Data_Flight{i,18} = Analysis_Data_Flight{i,18}';
            
            
            %*************** NEW BINNING ARROUND OFFSET:
            %New Try of binning: until exact offset (3):
            Analysis_Data_Flight{i,19} = flight_offset_trigger(i-bout_counter)-analysis_window_pre_offset:binsize:flight_offset_trigger(i-bout_counter);
            
            %New Try of binning: from offset to end of time window (4):
            Analysis_Data_Flight{i,20} = flight_offset_trigger(i-bout_counter):binsize:flight_offset_trigger(i-bout_counter)+analysis_window_post_offset;
            
            %NEW: Spikes within these time windows
            %** BINNED (3)
            for z=1:length(Analysis_Data_Flight{i,19})-1
                %calculates number of spikes wihtin each bin
                Analysis_Data_Flight{i,21}(z) = nansum(spikes(Analysis_Data_Flight{i,19}(z):Analysis_Data_Flight{i,19}(z+1)));
            end
            %** (3) UNBINNED = SPIKEVECTOR
            Analysis_Data_Flight{i,22} = spikes(Analysis_Data_Flight{i,19}(1):Analysis_Data_Flight{i,19}(end)-1);
            Analysis_Data_Flight{i,22}(isnan(Analysis_Data_Flight{i,22}))=0;
            Analysis_Data_Flight{i,22} = Analysis_Data_Flight{i,22}';
            
            %** (4)
            for z=1:length(Analysis_Data_Flight{i,20})-1
                %calculates number of spikes wihtin each bin
                Analysis_Data_Flight{i,23}(z) = nansum(spikes(Analysis_Data_Flight{i,20}(z):Analysis_Data_Flight{i,20}(z+1)));
            end
            %** (4) UNBINNED = SPIKEVECTOR
            Analysis_Data_Flight{i,24} = spikes(Analysis_Data_Flight{i,20}(1):Analysis_Data_Flight{i,20}(end)-1);
            Analysis_Data_Flight{i,24}(isnan(Analysis_Data_Flight{i,24}))=0;
            Analysis_Data_Flight{i,24} = Analysis_Data_Flight{i,24}';
            
            if i == bout_counter+length(flight_onset_trigger) %Used later in averaged file to assign flyID
                idx_fly_counter_clean = idx_fly_counter_clean +1;
            end
            
            %Restore flying Onset Time
            if analysis_window_pre_onset_back == 0 %In case Fly is flying in the beginning
                analysis_window_pre_onset = analysis_window_pre_onset_back;
            end
            
            %Save Spike Trace
            Analysis_Data_Flight{i,25} = spikes;
            
            
            
            % NEW Calculate the membrane potential shortly before and after Onset
            %of LED Stim
            analysis_membrane_potential_prior_onset = analysis_membrane_potential_prior_onset_in_sec*sampling_rate;
            analysis_membrane_potential_post_onset = analysis_membrane_potential_post_onset_in_sec*sampling_rate;
            analysis_membrane_potential_prior_offset = analysis_membrane_potential_prior_offset_in_sec*sampling_rate;
            analysis_membrane_potential_post_offset = analysis_membrane_potential_post_offset_in_sec*sampling_rate;
            
            
            membrane_potential_prior_onset = median(intra_orig(flight_onset_trigger(i-bout_counter)-analysis_membrane_potential_prior_onset:flight_onset_trigger(i-bout_counter)));
            Analysis_Data_Flight{i,26} = membrane_potential_prior_onset;
            membrane_potential_post_onset = median(intra_orig(flight_onset_trigger(i-bout_counter):flight_onset_trigger(i-bout_counter)+analysis_membrane_potential_post_onset));
            Analysis_Data_Flight{i,27} = membrane_potential_post_onset;
            
            membrane_potential_prior_offset = median(intra_orig(flight_offset_trigger(i-bout_counter)-analysis_membrane_potential_prior_offset:flight_offset_trigger(i-bout_counter)));
            Analysis_Data_Flight{i,28} = membrane_potential_prior_offset;
            membrane_potential_post_offset = median(intra_orig(flight_offset_trigger(i-bout_counter):flight_offset_trigger(i-bout_counter)+analysis_membrane_potential_post_offset));
            Analysis_Data_Flight{i,29} = membrane_potential_post_offset;
            
            %New Intra Spike Traces before on and offset and corresponding
            %Velocity Vectors
            Analysis_Data_Flight{i,30} = intra_orig(flight_onset_trigger(i-bout_counter)-analysis_membrane_potential_prior_onset:flight_onset_trigger(i-bout_counter));%Prior Onset
            Analysis_Data_Flight{i,31} = intra_orig(flight_onset_trigger(i-bout_counter):flight_onset_trigger(i-bout_counter)+analysis_membrane_potential_post_onset);%Post Onset
            Analysis_Data_Flight{i,32} = intra_orig(flight_offset_trigger(i-bout_counter)-analysis_membrane_potential_prior_offset:flight_offset_trigger(i-bout_counter));%Prior Offset
            Analysis_Data_Flight{i,33} = intra_orig(flight_offset_trigger(i-bout_counter):flight_offset_trigger(i-bout_counter)+analysis_membrane_potential_post_offset);%Post Offset
            
            Analysis_Data_Flight{i,34} = tacho(flight_onset_trigger(i-bout_counter)-analysis_membrane_potential_prior_onset:flight_onset_trigger(i-bout_counter));%Prior Onset
            Analysis_Data_Flight{i,35} = tacho(flight_onset_trigger(i-bout_counter):flight_onset_trigger(i-bout_counter)+analysis_membrane_potential_post_onset);%Post Onset
            Analysis_Data_Flight{i,36} = tacho(flight_offset_trigger(i-bout_counter)-analysis_membrane_potential_prior_offset:flight_offset_trigger(i-bout_counter));%Prior Offset
            Analysis_Data_Flight{i,37} = tacho(flight_offset_trigger(i-bout_counter):flight_offset_trigger(i-bout_counter)+analysis_membrane_potential_post_offset);%Post Offset
            
            Analysis_Data_Flight{i,38} = 0;%LED_Stim(flight_onset_trigger-analysis_membrane_potential_prior_onset:flight_onset_trigger);%Prior Onset; was LED stim
            Analysis_Data_Flight{i,39} = 0;%LED_Stim(flight_onset_trigger:flight_onset_trigger+analysis_membrane_potential_post_onset);%Post Onset;  was LED stim
            Analysis_Data_Flight{i,40} = 0;%LED_Stim(flight_offset_trigger-analysis_membrane_potential_prior_offset:flight_offset_trigger);%Prior Offset;  was LED stim
            Analysis_Data_Flight{i,41} = 0;%LED_Stim(flight_offset_trigger:flight_offset_trigger+analysis_membrane_potential_post_offset);%Post Offset;  was LED stim
            
            Analysis_Data_Flight{i,42} = 0;  %was LED stim
            Analysis_Data_Flight{i,43} = intra_orig;
            Analysis_Data_Flight{i,44} = flight_onset_trigger(i-bout_counter);
            Analysis_Data_Flight{i,45} = flight_offset_trigger(i-bout_counter);
            Analysis_Data_Flight{i,46} = tacho;
            Analysis_Data_Flight{i,47} = tacho;
            Analysis_Data_Flight{i,48} = flightthreshold;
            %filtered intra trace
            Analysis_Data_Flight{i,49} = medfilt1(intra_orig, smoothing_factor_membrane);
            

            Analysis_Data_Flight{i,50} = 0;
            %fly_activity_offset_clean_corresponding_activity_offset
            Analysis_Data_Flight{i,51} = 0;
            clearvars nearest_onset nearest_offset


        end
        
        %******************** NEW SUM OF SPIKES FOR ALL TRIALS
        spike_vector_per_flyID{idx_fly_counter_clean, 1} = flyID;
        %******** Putting together the single spike vectors PRE
        for i = bout_counter+1:length(Analysis_Data_Flight(:,1))
            temp_pre_onset(i-bout_counter,1:length(Analysis_Data_Flight{i,16}(1:end))) = Analysis_Data_Flight{i,16};
            temp_pre_offset(i-bout_counter,1:length(Analysis_Data_Flight{i,22}(1:end))) = Analysis_Data_Flight{i,22};
        end
        %******** Putting together the single spike vectors POST
        for i = bout_counter+1:length(Analysis_Data_Flight(:,1))
            temp_post_onset(i-bout_counter,1:length(Analysis_Data_Flight{i,18}(1:end))) = Analysis_Data_Flight{i,18};
            temp_post_offset(i-bout_counter,1:length(Analysis_Data_Flight{i,24}(1:end))) = Analysis_Data_Flight{i,24};
        end
        spike_vector_per_flyID{idx_fly_counter_clean,2} = temp_pre_onset;
        spike_vector_per_flyID{idx_fly_counter_clean,3} = temp_post_onset;
        spike_vector_per_flyID{idx_fly_counter_clean,4} = temp_pre_offset;
        spike_vector_per_flyID{idx_fly_counter_clean,5} = temp_post_offset;
        % SUM
        spike_vector_per_flyID{idx_fly_counter_clean,6} = sum(temp_pre_onset,1);
        spike_vector_per_flyID{idx_fly_counter_clean,7} = sum(temp_post_onset,1);
        spike_vector_per_flyID{idx_fly_counter_clean,8} = sum(temp_pre_offset,1);
        spike_vector_per_flyID{idx_fly_counter_clean,9} = sum(temp_post_offset,1);
        
        clearvars temp_pre_onset temp_pre_offset temp_post_onset temp_post_offset
        
        %********************* Current average Spikes per Bin ****************
        %15 = pre_onset, 17 = post_onset; 21 = pre_offset, 23 = post_offset
        %USING THE OLD BINNING DATA
        
        for i = bout_counter+1:length(Analysis_Data_Flight(:,1))
            temp_binned_average_spikes_pre_onset(i,1:length(Analysis_Data_Flight{i,15}(1:end))) = Analysis_Data_Flight{i,15};
            temp_binned_average_spikes_post_onset(i,1:length(Analysis_Data_Flight{i,17}(1:end))) = Analysis_Data_Flight{i,17};
            temp_binned_average_spikes_pre_offset(i,1:length(Analysis_Data_Flight{i,21}(1:end))) = Analysis_Data_Flight{i,21};
            temp_binned_average_spikes_post_offset(i,1:length(Analysis_Data_Flight{i,23}(1:end))) = Analysis_Data_Flight{i,23};
        end
        average_spikes_flyID{idx_fly_counter_clean, 1} = flyID;
        average_spikes_flyID{idx_fly_counter_clean, 2} = (sum(temp_binned_average_spikes_pre_onset,1))./length(bout_counter+1:length(Analysis_Data_Flight(:,1)));%./length(bout_counter+1:length(Analysis_Data_Flight(:,1)))  to average in the number of walking bouts
        average_spikes_flyID{idx_fly_counter_clean, 3} = (sum(temp_binned_average_spikes_post_onset,1))./length(bout_counter+1:length(Analysis_Data_Flight(:,1)));%./length(bout_counter+1:length(Analysis_Data_Flight(:,1)))*2;
        average_spikes_flyID{idx_fly_counter_clean, 4} = (sum(temp_binned_average_spikes_pre_offset,1))./length(bout_counter+1:length(Analysis_Data_Flight(:,1)));%./length(bout_counter+1:length(Analysis_Data_Flight(:,1)))*2;
        average_spikes_flyID{idx_fly_counter_clean, 5} = (sum(temp_binned_average_spikes_post_offset,1))./length(bout_counter+1:length(Analysis_Data_Flight(:,1)));%./length(bout_counter+1:length(Analysis_Data_Flight(:,1)))*2;
        
        clearvars temp_binned_average_spikes_pre_onset temp_binned_average_spikes_post_onset temp_binned_average_spikes_pre_offset temp_binned_average_spikes_post_offset
        
        
        
        
        
        
        
        
        
        
        
        
        if checkbox_manual_evaluation == 1 %ONLY EXECUTED WHEN SET TO FULLY AUTOMATED ANALYSIS
            %Simple Plotting Plotting Results of walking onset
            h = figure('name',[flyID ' Fly Onset'], 'Color', figure_plotting_color);
            set(h, 'position', [1, 600, 450, 450]);
            hold on
            for z=bout_counter+1:length(Analysis_Data_Flight(:,1))
                %Plot Data
                x_axis = (1:length(Analysis_Data_Flight{z,4}(1:end)))/sampling_rate;
                plot(x_axis, Analysis_Data_Flight{z,4}(1:end))
                %Plot walking onset
                plot(analysis_window_pre_onset./sampling_rate,min(Analysis_Data_Flight{z,4}):max(Analysis_Data_Flight{z,4}), '|', 'MarkerSize',10, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
            end
            hold off
            
            %Simple plotting of Results of walking offset
            h = figure('name',[flyID ' Fly Offset']);
            set(h, 'position', [1, 600, 450, 450]);
            hold on
            for z=bout_counter+1:length(Analysis_Data_Flight(:,1))
                %Plot Data
                x_axis = (1:length(Analysis_Data_Flight{z,8}(1:end)))/sampling_rate;
                plot(x_axis, Analysis_Data_Flight{z,8}(1:end))
                %Plot walking offset
                plot(analysis_window_pre_offset./sampling_rate,min(Analysis_Data_Flight{z,8}):max(Analysis_Data_Flight{z,8}), '|', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1], 'LineWidth', 4);
                
            end
            hold off
        end
        
        %% ************* Plotting All Results and Discard Trials? *************
        if checkbox_manual_evaluation == 1 %ONLY EXECUTED WHEN SET TO FULLY AUTOMATED ANALYSIS
            indicator_analysis_correct = 0;
            promt_showed_up = 0;
            %BACKUP FOR RETURN
            Analysis_Data_Flight_Back = Analysis_Data_Flight;
            while indicator_analysis_correct == 0
                %In Case All Trials have been deleted:
                if bout_counter+1 > length(Analysis_Data_Flight(:,1))
                    break
                end
                
                %Plotting of all fly bouts (onset) of this MAT file under another
                h = figure('name',[flyID ' All walking Bouts (ONSET)']);
                set(h, 'position', [1, 600, 450, 450]);
                t = tiledlayout(length(bout_counter+1:length(Analysis_Data_Flight(:,1))),1);
                % Variables needed
                total_min = [];
                total_max = [];
                
                for u = bout_counter+1:length(Analysis_Data_Flight(:,1))
                    %Determine max Y Value in the files to plot for yLim
                    if u == bout_counter+1
                        total_max = max(Analysis_Data_Flight{u,4});
                        total_min = min(Analysis_Data_Flight{u,4});
                        for  z = bout_counter+1:length(Analysis_Data_Flight(:,1))
                            current_min = min(Analysis_Data_Flight{z,4});
                            current_max = max(Analysis_Data_Flight{z,4});
                            if current_min < total_min ==1
                                total_min = current_min;
                            end
                            
                            if current_max > total_max ==1
                                total_max = current_max;
                            end
                        end
                    end
                    
                    nexttile
                    %Plot Intra Trace
                    x_axis = (1:length(Analysis_Data_Flight{u,4}(1:end)))/sampling_rate;
                    plot(x_axis, Analysis_Data_Flight{u,4}(1:end))
                    hold on
                    %Plot Start of walking
                    plot(analysis_window_pre_onset./sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
                    %Plot end of walking if in range
                    if flight_offset_trigger(u-bout_counter)-flight_onset_trigger(u-bout_counter)+analysis_window_pre_onset<length(Analysis_Data_Flight{u,4})
                        %                 plot((fly_flight_offset_trigger(u-bout_counter)-fly_flight_onset_trigger(u-bout_counter)+analysis_window_pre_stiming)/sampling_rate, min(Analysis_Data_Flight{u,4}):max(Analysis_Data_Flight{u,4}), '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1], 'LineWidth', 2);
                        plot(abs(flight_offset_trigger(u-bout_counter)-flight_onset_trigger(u-bout_counter)+analysis_window_pre_onset)/sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1], 'LineWidth', 2);
                        
                    end
                    ax = gca;
                    if u==length(Analysis_Data_Flight(:,1))
                        set(ax,'YTick', []);
                    else
                        set(ax,'XTick',[], 'YTick', []);
                    end
                    
                    %Clicking on the subplot using ginput later returns number of
                    %subplot to delete data
                    set(gca,'tag',num2str(u));
                    
                    ylim([total_min,total_max]);
                    
                end
                clearvars total_min total_max current_min current_max
                
                %%Decrease space between subplots
                t.Padding = 'none';
                t.TileSpacing = 'none';
                hold off
                
                % Variables needed
                total_min = [];
                total_max = [];
                current_max = [];
                current_min = [];
                
                
                %Plotting of all fly bouts (OFFSET) of this MAT file under another
                f = figure('name',[flyID ' All walking Bouts (OFFSET)']);
                
                tf = tiledlayout(length(bout_counter+1:length(Analysis_Data_Flight(:,1))),1);
                for u = bout_counter+1:length(Analysis_Data_Flight(:,1))
                    nexttile
                    
                    %Determine max Y Value in the files to plot for yLim
                    if u == bout_counter+1
                        total_max = max(Analysis_Data_Flight{u,8});
                        total_min = min(Analysis_Data_Flight{u,8});
                        for  z = bout_counter+1:length(Analysis_Data_Flight(:,1))
                            current_min = min(Analysis_Data_Flight{z,8});
                            current_max = max(Analysis_Data_Flight{z,8});
                            if current_min < total_min ==1
                                total_min = current_min;
                            end
                            
                            if current_max > total_max ==1
                                total_max = current_max;
                            end
                        end
                    end
                    
                    %Plot Intra Trace
                    x_axis = (1:length(Analysis_Data_Flight{u,8}(1:end)))/sampling_rate;
                    plot(x_axis, Analysis_Data_Flight{u,8}(1:end))
                    hold on
                    %Plot end of walking
                    %           plot(x_axis(end)/2, min(Analysis_Data_Flight{u,8}):max(Analysis_Data_Flight{u,8}), '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
                    plot(analysis_window_pre_offset./sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1], 'LineWidth', 2);
                    
                    %Plot start of walking if in range
                    %                     if fly_flight_offset_trigger(u-bout_counter)-walking_onset_trigger(u-bout_counter)+analysis_window_pre_stiming<length(Analysis_Data_Flight{u,8})
                    %                         %                 plot((walking_onset_trigger(u-bout_counter)-fly_flight_offset_trigger(u-bout_counter)+analysis_window_pre_stiming)/sampling_rate, min(Analysis_Data_Flight{u,8}):max(Analysis_Data_Flight{u,8}), '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1]);
                    %                         plot((abs((walking_onset_trigger(u-bout_counter)-fly_flight_offset_trigger(u-bout_counter)+analysis_window_pre_stiming)/sampling_rate)), total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
                    %
                    %                     end
                    
                    
                    
                    if flight_offset_trigger(u-bout_counter)-flight_onset_trigger(u-bout_counter)<analysis_window_pre_offset
                        %                 plot((fly_flight_onset_trigger(u-bout_counter)-fly_flight_offset_trigger(u-bout_counter)+analysis_window_pre_stiming)/sampling_rate, min(Analysis_Data_Flight{u,8}):max(Analysis_Data_Flight{u,8}), '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1]);
                        %                         plot(((fly_flight_onset_trigger(u-bout_counter)-fly_flight_offset_trigger(u-bout_counter)+analysis_window_pre_stiming)/sampling_rate), total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
                        plot(analysis_window_pre_offset./sampling_rate-(flight_offset_trigger(u-bout_counter)-flight_onset_trigger(u-bout_counter))/sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
                    end
                    ax = gca;
                    if u==length(Analysis_Data_Flight(:,1))
                        set(ax,'YTick', []);
                    else
                        set(ax,'XTick',[], 'YTick', []);
                    end
                    
                    %Clicking on the subplot using ginput later returns number of
                    %subplot to delete data
                    set(gca,'tag',num2str(u));
                    
                    ylim([total_min,total_max]);
                end
                clearvars total_min total_max current_min current_max
                
                %%Decrease space between subplots
                tf.Padding = 'none';
                tf.TileSpacing = 'none';
                hold off
                
                %% ************* Exclude data from Dataset *************
                
                indicator_discard = 0;
                indicator_undo = 0;
                
                % Promt to ask whether the bouts shall be deleted or included in analysis
                dlgTitle    = 'Discard Trials from the Dataset?';
                dlgQuestion = 'Do you wish to ACCEPT, DISCARD Trials from the Analysis or UNDO Changes?';
                choice = questdlg(dlgQuestion,dlgTitle,'ACCEPT','DISCARD', 'UNDO', 'ACCEPT');
                if length(choice)==6  %Analysis is correct
                    indicator_analysis_correct=1;
                elseif length(choice)==7 %Discard Trials
                    indicator_discard=1;
                elseif length(choice)==4 %Add Trials
                    indicator_undo=1;
                end
                
                promt_showed_up = 1;% Used later on to be able to discard trials during last plotting
                
                if promt_showed_up == 1 % Used later on to be able to discard trials during this last plotting
                    %******************** DISCARDING TRIALS**************************%
                    if indicator_discard == 1
                        %Select Trials by clicking at the beginning and end of a walking bout,  return key for
                        %confirm
                        idx_counter_trial = 0;
                        while 1 == 1
                            w = waitforbuttonpress;
                            switch w
                                case 1 % keyboard
                                    key = get(gcf,'currentcharacter');
                                    if key==27 % (the Esc key)
                                        try; delete(h); end
                                        break
                                    end
                                case 0 % mouse click
                                    mousept = get(gca,'currentPoint');
                                    x = mousept(1,1);
                                    y = mousept(1,2);
                                    %               try; delete(h); end
                                    h = text(x,y,get(gca,'tag'),'vert','middle','horiz','center');
                                    idx_delete_trial(idx_counter_trial+1) = str2num(h.String);
                                    idx_counter_trial = idx_counter_trial + 1;
                            end
                        end
                        
                        close(gcf);
                        close(gcf);
                        pause(1);
                        
                        if exist('idx_delete_trial') == 1
                            Analysis_Data_Flight(idx_delete_trial,:) = [];
                            clearvars idx_delete_trial
                        end
                    end
                    
                    %******************** UNDO **************************%
                    if indicator_undo == 1
                        Analysis_Data_Flight = Analysis_Data_Flight_Back;
                    end
                end
                
            end
            close all
        end
        
        if input_save_figures == 0
            close all
        end
        
    end
end
% % % % % % end


%% LEARING WORKSPACE
clearvars current_orig velocX velocY velocZ clean_activity flight fly_activity intra_orig intra_smooth LED_Stim spikes tacho Velocity_Data velocXYZ


%%
%********************* Calculate Median Membrane Amplitude before and after offset
membrane_potential_prior_onset = [Analysis_Data_Flight{:,26}];
membrane_potential_post_onset = [Analysis_Data_Flight{:,27}];
membrane_potential_prior_offset = [Analysis_Data_Flight{:,28}];
membrane_potential_post_offset = [Analysis_Data_Flight{:,29}];

%percentile for plotting 
lowerprctl_membrane_potential_prior_onset = prctile(membrane_potential_prior_onset, 25, 'all');
upper_prctlmembrane_potential_prior_onset = prctile(membrane_potential_prior_onset, 75, 'all');
lowerprctl_membrane_potential_post_onset = prctile(membrane_potential_post_onset, 25, 'all');
upper_prctl_membrane_potential_post_onset = prctile(membrane_potential_post_onset, 75, 'all');

lowerprctl_membrane_potential_prior_offset = prctile(membrane_potential_prior_offset, 25, 'all');
upper_prctl_membrane_potential_prior_offset = prctile(membrane_potential_prior_offset, 75, 'all');
lowerprctl_membrane_potential_post_offset = prctile(membrane_potential_post_offset, 25, 'all');
upper_prctl_membrane_potential_post_offset = prctile(membrane_potential_post_offset, 75, 'all');

membranepotential_figure = figure('Color', figure_plotting_color);
hold on

%plotting Median
plot([-0.4 0.4], [median(membrane_potential_prior_onset) median(membrane_potential_prior_onset)], 'LineWidth', 2, 'Color', 'r')
plot([0.6 1.4], [median(membrane_potential_post_onset) median(membrane_potential_post_onset)], 'LineWidth', 2, 'Color', 'r')

plot([1.6 2.4], [median(membrane_potential_prior_offset) median(membrane_potential_prior_offset)], 'LineWidth', 2, 'Color', 'r')
plot([2.6 3.4], [median(membrane_potential_post_offset) median(membrane_potential_post_offset)], 'LineWidth', 2, 'Color', 'r')

%Plot Percentile %NEW
plot([-0.4 0.4], [lowerprctl_membrane_potential_prior_onset lowerprctl_membrane_potential_prior_onset], 'LineWidth', 2, 'Color', 'b')
plot([-0.4 0.4], [upper_prctlmembrane_potential_prior_onset upper_prctlmembrane_potential_prior_onset], 'LineWidth', 2, 'Color', 'b')
plot([0.6 1.4], [lowerprctl_membrane_potential_post_onset lowerprctl_membrane_potential_post_onset], 'LineWidth', 2, 'Color', 'b')
plot([0.6 1.4], [upper_prctl_membrane_potential_post_onset upper_prctl_membrane_potential_post_onset], 'LineWidth', 2, 'Color', 'b')
plot([1.6 2.4], [lowerprctl_membrane_potential_prior_offset lowerprctl_membrane_potential_prior_offset], 'LineWidth', 2, 'Color', 'b')
plot([1.6 2.4], [upper_prctl_membrane_potential_prior_offset upper_prctl_membrane_potential_prior_offset], 'LineWidth', 2, 'Color', 'b')
plot([2.6 3.4], [lowerprctl_membrane_potential_post_offset lowerprctl_membrane_potential_post_offset], 'LineWidth', 2, 'Color', 'b')
plot([2.6 3.4], [upper_prctl_membrane_potential_post_offset upper_prctl_membrane_potential_post_offset], 'LineWidth', 2, 'Color', 'b')


temp1 = [zeros(length(membrane_potential_prior_onset),1) ones(length(membrane_potential_prior_onset),1)];
temp2 = [membrane_potential_prior_onset' membrane_potential_post_onset'];

temp3 = [ones(length(membrane_potential_prior_offset),1)*2 ones(length(membrane_potential_prior_offset),1)*3];
temp4 = [membrane_potential_prior_offset' membrane_potential_post_offset'];

for i = 1:2
    if i == 1
        color = 'b';
    else
        color = 'm';
    end
    swarmchart(temp1(:,i), temp2(:,i), 30, color, 'filled');
    alpha(0.5)
end

for i = 1:2
    if i == 1
        color = 'b';
    else
        color = 'm';
    end
    swarmchart(temp3(:,i), temp4(:,i), 30, color, 'filled');
    alpha(0.5)
end

X = categorical({'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'});
X = reordercats(X,{'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'});
set(gca,'XTick',0:3,'XTickLabel',X)

%setting axes and enviroment of figur %2021_08_13
ylabel('Vm (mV)');
if plotting_figures_in_black == 1
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',12);
    ax_membranepotential_figure = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_membranepotential_figure,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_membranepotential_figure,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
else
    ax_membranepotential_figure = gca;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',12);
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_membranepotential_figure,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_membranepotential_figure,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
end


[A, h] = signrank(membrane_potential_prior_onset, membrane_potential_post_onset);
[B, h] = signrank(membrane_potential_prior_offset, membrane_potential_post_offset);
% Create text file with stats
T = table(A, B, 'VariableNames', { 'Prior-Onset vs. Post-Onset', 'Prior-Offset vs. Post-Offset'} )
% Write data to text file
writetable(T, 'Stats_Delta_Membranepotential.txt')

fid = fopen('Stats_Delta_Membranepotential2.txt','w');
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_prior_onset: ', num2str(median(membrane_potential_prior_onset)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_post_onset: ', num2str(median(membrane_potential_post_onset)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_prior_offset: ', num2str(median(membrane_potential_prior_offset)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_post_offset: ', num2str(median(membrane_potential_post_offset)));

fprintf(fid, '%6s %12s\r\n', 'lowerprctl_membrane_potential_prior_onset: ', num2str(lowerprctl_membrane_potential_prior_onset));
fprintf(fid, '%6s %12s\r\n', 'upper_prctlmembrane_potential_prior_onset: ', num2str(upper_prctlmembrane_potential_prior_onset));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_membrane_potential_post_onset: ', num2str(lowerprctl_membrane_potential_post_onset));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_membrane_potential_post_onset: ', num2str(upper_prctl_membrane_potential_post_onset));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_membrane_potential_prior_offset: ', num2str(lowerprctl_membrane_potential_prior_offset));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_membrane_potential_prior_offset: ', num2str(upper_prctl_membrane_potential_prior_offset));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_membrane_potential_post_offset: ', num2str(lowerprctl_membrane_potential_post_offset));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_membrane_potential_post_offset: ', num2str(upper_prctl_membrane_potential_post_offset));

fclose(fid);

filename = [file_modifier '_Delta_Membranepotential.tif'];
saveas(gcf, filename)
filename = [file_modifier '_Delta_Membranepotential.eps'];
print(filename, '-depsc2', '-tiff', '-r300', '-painters')
pause(1)
close gcf

%**************Same with Normalization
membrane_potential_prior_onset_norm = ones(length([Analysis_Data_Flight{:,26}]), 1)';
membrane_potential_post_onset_norm = membrane_potential_post_onset./membrane_potential_prior_onset;

membrane_potential_prior_offset_norm = membrane_potential_prior_offset./membrane_potential_prior_onset;
membrane_potential_post_offset_norm = membrane_potential_post_offset./membrane_potential_prior_onset;

%percentile for plotting 
lowerprctl_membrane_potential_prior_onset_norm = prctile(membrane_potential_prior_onset_norm, 25, 'all');
upper_prctl_membrane_potential_prior_onset_norm = prctile(membrane_potential_prior_onset_norm, 75, 'all');
lowerprctl_membrane_potential_post_onset_norm = prctile(membrane_potential_post_onset_norm, 25, 'all');
upper_prctl_membrane_potential_post_onset_norm = prctile(membrane_potential_post_onset_norm, 75, 'all');

lowerprctl_membrane_potential_prior_offset_norm = prctile(membrane_potential_prior_offset_norm, 25, 'all');
upper_prctl_membrane_potential_prior_offset_norm = prctile(membrane_potential_prior_offset_norm, 75, 'all');
lowerprctl_membrane_potential_post_offset_norm = prctile(membrane_potential_post_offset_norm, 25, 'all');
upper_prctl_membrane_potential_post_offset_norm = prctile(membrane_potential_post_offset_norm, 75, 'all');

membranepotential_figure = figure('Color', figure_plotting_color);
hold on

%plotting Median
plot([-0.4 0.4], [median(membrane_potential_prior_onset_norm) median(membrane_potential_prior_onset_norm)], 'LineWidth', 2, 'Color', 'r')
plot([0.6 1.4], [median(membrane_potential_post_onset_norm) median(membrane_potential_post_onset_norm)], 'LineWidth', 2, 'Color', 'r')

plot([1.6 2.4], [median(membrane_potential_prior_offset_norm) median(membrane_potential_prior_offset_norm)], 'LineWidth', 2, 'Color', 'r')
plot([2.6 3.4], [median(membrane_potential_post_offset_norm) median(membrane_potential_post_offset_norm)], 'LineWidth', 2, 'Color', 'r')

%Plot Percentile 
plot([-0.4 0.4], [lowerprctl_membrane_potential_prior_onset_norm lowerprctl_membrane_potential_prior_onset_norm], 'LineWidth', 2, 'Color', 'b')
plot([-0.4 0.4], [upper_prctl_membrane_potential_prior_onset_norm upper_prctl_membrane_potential_prior_onset_norm], 'LineWidth', 2, 'Color', 'b')
plot([0.6 1.4], [lowerprctl_membrane_potential_post_onset_norm lowerprctl_membrane_potential_post_onset_norm], 'LineWidth', 2, 'Color', 'b')
plot([0.6 1.4], [upper_prctl_membrane_potential_post_onset_norm upper_prctl_membrane_potential_post_onset_norm], 'LineWidth', 2, 'Color', 'b')
plot([1.6 2.4], [lowerprctl_membrane_potential_prior_offset_norm lowerprctl_membrane_potential_prior_offset_norm], 'LineWidth', 2, 'Color', 'b')
plot([1.6 2.4], [upper_prctl_membrane_potential_prior_offset_norm upper_prctl_membrane_potential_prior_offset_norm], 'LineWidth', 2, 'Color', 'b')
plot([2.6 3.4], [lowerprctl_membrane_potential_post_offset_norm lowerprctl_membrane_potential_post_offset_norm], 'LineWidth', 2, 'Color', 'b')
plot([2.6 3.4], [upper_prctl_membrane_potential_post_offset_norm upper_prctl_membrane_potential_post_offset_norm], 'LineWidth', 2, 'Color', 'b')


temp1 = [zeros(length(membrane_potential_prior_onset_norm),1) ones(length(membrane_potential_prior_onset_norm),1)];
temp2 = [membrane_potential_prior_onset_norm' membrane_potential_post_onset_norm'];

temp3 = [ones(length(membrane_potential_prior_offset_norm),1)*2 ones(length(membrane_potential_prior_offset_norm),1)*3];
temp4 = [membrane_potential_prior_offset_norm' membrane_potential_post_offset_norm'];

for i = 1:2
    if i == 1
        color = 'b';
    else
        color = 'm';
    end
    swarmchart(temp1(:,i), temp2(:,i), 30, color, 'filled');
    alpha(0.5)
end

for i = 1:2
    if i == 1
        color = 'b';
    else
        color = 'm';
    end
    swarmchart(temp3(:,i), temp4(:,i), 30, color, 'filled');
    alpha(0.5)
end
axis ij
X = categorical({'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'});
X = reordercats(X,{'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'});
set(gca,'XTick',0:3,'XTickLabel',X)
% set(gca,'Ydir','normal')

%setting axes and enviroment of figur 
ylabel('Delta Vm (mV)');
if plotting_figures_in_black == 1
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',12);
    ax_membranepotential_figure = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_membranepotential_figure,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_membranepotential_figure,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
else
    ax_membranepotential_figure = gca;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',12);
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_membranepotential_figure,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_membranepotential_figure,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
end


[A, h] = signrank(membrane_potential_prior_onset_norm, membrane_potential_post_onset_norm);
[B, h] = signrank(membrane_potential_prior_offset_norm, membrane_potential_post_offset_norm);
% Create text file with stats
T = table(A, B, 'VariableNames', { 'Prior-Onset vs. Post-Onset', 'Prior-Offset vs. Post-Offset'} )
% Write data to text file
writetable(T, 'Stats_Delta_Membranepotential_norm.txt')

fid = fopen('Stats_Delta_Membranepotential_norm2.txt','w');
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_prior_onset_norm: ', num2str(median(membrane_potential_prior_onset_norm)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_post_onset_norm: ', num2str(median(membrane_potential_post_onset_norm)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_prior_offset_norm: ', num2str(median(membrane_potential_prior_offset_norm)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_post_offset_norm: ', num2str(median(membrane_potential_post_offset_norm)));


fprintf(fid, '%6s %12s\r\n', 'lowerprctl_membrane_potential_prior_onset_norm: ', num2str(lowerprctl_membrane_potential_prior_onset_norm));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_membrane_potential_prior_onset_norm: ', num2str(upper_prctl_membrane_potential_prior_onset_norm));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_membrane_potential_post_onset_norm: ', num2str(lowerprctl_membrane_potential_post_onset_norm));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_membrane_potential_post_onset_norm: ', num2str(upper_prctl_membrane_potential_post_onset_norm));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_membrane_potential_prior_offset_norm: ', num2str(lowerprctl_membrane_potential_prior_offset_norm));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_membrane_potential_prior_offset_norm: ', num2str(upper_prctl_membrane_potential_prior_offset_norm));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_membrane_potential_post_offset_norm: ', num2str(lowerprctl_membrane_potential_post_offset_norm));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_membrane_potential_post_offset_norm: ', num2str(upper_prctl_membrane_potential_post_offset_norm));


fclose(fid);

filename = [file_modifier '_Delta_Membranepotential_norm.tif'];
saveas(gcf, filename)
filename = [file_modifier '_Delta_Membranepotential_norm.eps'];
print(filename, '-depsc2', '-tiff', '-r300', '-painters')
pause(1)
close gcf

%% BASELINE SUBSTRACTION
membrane_potential_prior_onset_norm2 = zeros(length([Analysis_Data_Flight{:,26}]), 1)';
membrane_potential_post_onset_norm2 = membrane_potential_post_onset-membrane_potential_prior_onset;
membrane_potential_prior_offset_norm2 = membrane_potential_prior_offset-membrane_potential_prior_onset;
membrane_potential_post_offset_norm2 = membrane_potential_post_offset-membrane_potential_prior_onset;

%percentile for plotting
lowerprctl_prior_onset_norm2 = prctile(membrane_potential_prior_onset_norm2, 25, 'all');
upper_prctl_prior_onset_norm2 = prctile(membrane_potential_prior_onset_norm2, 75, 'all');
lowerprctl_post_onset_norm2 = prctile(membrane_potential_post_onset_norm2, 25, 'all');
upper_prctl_post_onset_norm2 = prctile(membrane_potential_post_onset_norm2, 75, 'all');

lowerprctl_prior_offset_norm2 = prctile(membrane_potential_prior_offset_norm2, 25, 'all');
upper_prctl_prior_offset_norm2 = prctile(membrane_potential_prior_offset_norm2, 75, 'all');
lowerprctl_post_offset_norm2 = prctile(membrane_potential_post_offset_norm2, 25, 'all');
upper_prctl_post_offset_norm2 = prctile(membrane_potential_post_offset_norm2, 75, 'all');

figure
hold on
%Plot Median
plot([-0.4 0.4], [median(membrane_potential_prior_onset_norm2) median(membrane_potential_prior_onset_norm2)], 'LineWidth', 2, 'Color', 'r')
plot([0.6 1.4], [median(membrane_potential_post_onset_norm2) median(membrane_potential_post_onset_norm2)], 'LineWidth', 2, 'Color', 'r')
plot([1.6 2.4], [median(membrane_potential_prior_offset_norm2) median(membrane_potential_prior_offset_norm2)], 'LineWidth', 2, 'Color', 'r')
plot([2.6 3.4], [median(membrane_potential_post_offset_norm2) median(membrane_potential_post_offset_norm2)], 'LineWidth', 2, 'Color', 'r')

%Plot Percentile 
plot([-0.4 0.4], [lowerprctl_prior_onset_norm2 lowerprctl_prior_onset_norm2], 'LineWidth', 2, 'Color', 'b')
plot([-0.4 0.4], [upper_prctl_prior_onset_norm2 upper_prctl_prior_onset_norm2], 'LineWidth', 2, 'Color', 'b')
plot([0.6 1.4], [lowerprctl_post_onset_norm2 lowerprctl_post_onset_norm2], 'LineWidth', 2, 'Color', 'b')
plot([0.6 1.4], [upper_prctl_post_onset_norm2 upper_prctl_post_onset_norm2], 'LineWidth', 2, 'Color', 'b')
plot([1.6 2.4], [lowerprctl_prior_offset_norm2 lowerprctl_prior_offset_norm2], 'LineWidth', 2, 'Color', 'b')
plot([1.6 2.4], [upper_prctl_prior_offset_norm2 upper_prctl_prior_offset_norm2], 'LineWidth', 2, 'Color', 'b')
plot([2.6 3.4], [lowerprctl_post_offset_norm2 lowerprctl_post_offset_norm2], 'LineWidth', 2, 'Color', 'b')
plot([2.6 3.4], [upper_prctl_post_offset_norm2 upper_prctl_post_offset_norm2], 'LineWidth', 2, 'Color', 'b')

temp1 = [zeros(length(membrane_potential_prior_onset_norm2),1) ones(length(membrane_potential_prior_onset_norm2),1)];
temp2 = [membrane_potential_prior_onset_norm2' membrane_potential_post_onset_norm2'];
temp3 = [ones(length(membrane_potential_prior_offset_norm2),1)*2 ones(length(membrane_potential_prior_offset_norm2),1)*3];
temp4 = [membrane_potential_prior_offset_norm2' membrane_potential_post_offset_norm2'];
for i = 1:2
if i == 1
color = 'b';
else
color = 'm';
end
swarmchart(temp1(:,i), temp2(:,i), 30, color, 'filled');
alpha(0.5)
end
for i = 1:2
if i == 1
color = 'b';
else
color = 'm';
end
swarmchart(temp3(:,i), temp4(:,i), 30, color, 'filled');
alpha(0.5)
end
axis ij
X = categorical({'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'});
X = reordercats(X,{'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'});
set(gca,'XTick',0:3,'XTickLabel',X)
% set(gca,'Ydir','normal') 

[A, h] = signrank(membrane_potential_prior_onset_norm2, membrane_potential_post_onset_norm2);
[B, h] = signrank(membrane_potential_prior_offset_norm2, membrane_potential_post_offset_norm2);
% Create text file with stats
T = table(A, B, 'VariableNames', { 'Prior-Onset vs. Post-Onset', 'Prior-Offset vs. Post-Offset'} )
% Write data to text file
writetable(T, 'Stats_Delta_Membranepotential_baselinesubstracted.txt')

fid = fopen('Stats_Delta_Membranepotential_baselinesubstracted2.txt','w');
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_prior_onset_norm2: ', num2str(median(membrane_potential_prior_onset_norm2)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_prior_onset: ', num2str(median(membrane_potential_prior_onset)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_post_onset_norm2: ', num2str(median(membrane_potential_post_onset_norm2)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_prior_offset_norm2: ', num2str(median(membrane_potential_prior_offset_norm2)));
fprintf(fid, '%6s %12s\r\n', 'membrane_potential_post_offset_norm2: ', num2str(median(membrane_potential_post_offset_norm2)));

fprintf(fid, '%6s %12s\r\n', 'lowerprctl_prior_onset_norm2: ', num2str(lowerprctl_prior_onset_norm2));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_prior_onset_norm2: ', num2str(upper_prctl_prior_onset_norm2));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_post_onset_norm2: ', num2str(lowerprctl_post_onset_norm2));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_post_onset_norm2: ', num2str(upper_prctl_post_onset_norm2));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_prior_offset_norm2: ', num2str(lowerprctl_prior_offset_norm2));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_prior_offset_norm2: ', num2str(upper_prctl_prior_offset_norm2));
fprintf(fid, '%6s %12s\r\n', 'lowerprctl_post_offset_norm2: ', num2str(lowerprctl_post_offset_norm2));
fprintf(fid, '%6s %12s\r\n', 'upper_prctl_post_offset_norm2: ', num2str(upper_prctl_post_offset_norm2));

fclose(fid);


filename = [file_modifier '_Delta_Membranepotential_baselinesubstracted.tif'];
saveas(gcf, filename)
filename = [file_modifier '_Delta_Membranepotential_baselinesubstracted.eps'];
print(filename, '-depsc2', '-tiff', '-r300', '-painters')
pause(1)
close gcf


%% Overlay of Intra Traces before and after Onset/Offset + Overlay of corresponding XYZ Velocity traces

%HARDCODED FOR THE NEXT TWO SECTIONS
y_axis_Vm_min = -40;
y_axis_Vm_max = -35;

% Selection of Examples for plotting in first subplot
selection_for_plotting = [2, 5, 7, 9, 15];

plotting_offset = 0;



% create colorvector for plotting
color = colormap(bone(length(Analysis_Data_Flight(:, 44))+20));
%cut off white and black extrems for better vizualisation
color = color(10:end-10, :);
close gcf

for z2 = 1:2

%************************Onset
membrane_figure_onset = figure('name', 'Membrane Figure Onset', 'Color', figure_plotting_color);
set(membrane_figure_onset, 'position', [1, 1, 250, 1500]);
clearvars y std_y
hold on
for i = 1:length(Analysis_Data_Flight(:, 44))
    sp1 = subplot(4,1,1);
    hold on
    %Prior Onset
    if any(selection_for_plotting(:)==i)
        x = (1:length(Analysis_Data_Flight{i, 43}(Analysis_Data_Flight{i, 44}-analysis_membrane_potential_prior_onset:Analysis_Data_Flight{i, 44})))/sampling_rate;
        y = Analysis_Data_Flight{i, 43}(Analysis_Data_Flight{i, 44}-analysis_membrane_potential_prior_onset:Analysis_Data_Flight{i, 44});
        y_back_for_text = y;
        
        %Adding plotting offset
        y = y + plotting_offset;
        plot(x,y, 'Color', color(i, :))
        
        %plot membrane pot of example
            xli = xlim;
            text(xli(1)+0.3, y(1),num2str(y_back_for_text(1)), 'HorizontalAlignment','right','VerticalAlignment','top', 'Color', 'k')%plot membrane pot
        
        %Post Onset
        x = (analysis_membrane_potential_prior_onset+1:analysis_membrane_potential_prior_onset+analysis_membrane_potential_post_onset+1)/sampling_rate;
        y = Analysis_Data_Flight{i, 43}(Analysis_Data_Flight{i, 44}:Analysis_Data_Flight{i, 44}+analysis_membrane_potential_post_onset);
        y_back_for_text = y;
        
        %Adding plotting offset 2021_08_25
        y = y + plotting_offset;
        plotting_offset = plotting_offset + 30;
        plot(x,y, 'Color', color(i, :))
        
        %plot membrane pot of example
            xli = xlim;
            text(xli(2)+0.3, y(end),num2str(y_back_for_text(end)), 'HorizontalAlignment','right','VerticalAlignment','top', 'Color', 'k')%plot membrane pot
    
    end
    %plot mean of intra
    %-----------Prior ONSET
    if i == length(Analysis_Data_Flight(:, 44))
        x = (1:length(mean(Analysis_Data_Flight{i, 43}(Analysis_Data_Flight{i, 44}-analysis_membrane_potential_prior_onset:Analysis_Data_Flight{i, 44}),2)))/sampling_rate;
        clear y
        clear y2
        
        %Mean Closeup of medfitted Curve with STD as Patch
        sp2 = subplot(4,1,2);
        hold on
        for z = 1:length(Analysis_Data_Flight(:, 49))
            if z2 == 1%medfilt
            y(:, z) = Analysis_Data_Flight{z, 49}(Analysis_Data_Flight{z, 44}-analysis_membrane_potential_prior_onset:Analysis_Data_Flight{z, 44});
            else %MEAN OF INTRA ORIG
             y(:, z) = Analysis_Data_Flight{z, 43}(Analysis_Data_Flight{z, 44}-analysis_membrane_potential_prior_onset:Analysis_Data_Flight{z, 44});   
            end
            temp_mean_membrane_potential_prior_onset(z) = mean(y(:, z), 1);
            median_onset(z) = median(y(end-0.5*sampling_rate:end, z));
            %*** plot Data of single Trials? NEW 2021_09_24
            y(:, z) = y(:, z)-median_onset(z);
%             plot(x, y(:,z), 'Color', [color(z, :) 0.3]) %NEW 2021_08_25
        end
        
        %Calculating Means (between trials and of all trials) and reduce plotting data
        mean_between_trials = mean(y, 2); %between
        mean_of_all_trials = mean(mean_between_trials, 1); %all trials
        mean_of_all_trials_onset = mean_of_all_trials;%"HERE" EVERYTHING with statement: has been subsequently added for substraction of baseline mean Vm
        
        
        %Calculating STD arround mean 
        for z = 1:length(y(:,1))
            std_y(z) = std(y(z,:));
        end
        std_y_downsample_x = linspace(x(1), x(end), 1000);
        std_y_downsample_y = interp1(x, std_y-mean_of_all_trials_onset, std_y_downsample_x);%HERE
        
        
        
        %Plotting Means
        xnew = linspace(x(1), x(end), 1000);
        ynew = interp1(x, mean_between_trials, xnew)-mean_of_all_trials_onset;%HERE
        plot(xnew, ynew, 'LineWidth', 3, 'Color', 'g')
        
        %dotted line as mean
        y2(1:length(x)) = mean_of_all_trials-mean_of_all_trials_onset;%HERE
        plot(x,y2, '--', 'Color', 'g', 'LineWidth', 3);
        clear y2
        
        %plotting std arround mean
        patchx =  [std_y_downsample_x,fliplr(std_y_downsample_x)];
        patchy = [ynew + std_y_downsample_y,fliplr(ynew)-fliplr(std_y_downsample_y)];
        p = patch(patchx, patchy, 'red', 'LineStyle','none');
        alpha(0.2)
        
        %setting axis
        ylim([min(mean_between_trials)-2 max(mean_between_trials)+2 ])
        
        %Save means and STD for writing it into txt file later on
        mean_membrane_potential_prior_onset = mean(temp_mean_membrane_potential_prior_onset, 2);
        mean_membrane_potential_prior_onset_STD =  mean(std_y, 2);
        
        mean_of_all_trials_prior_onset = mean_of_all_trials;
        
        clearvars patchx patchy y2 xnew ynew std_y_downsample_x std_y_downsample_y mean_between_trials mean_of_all_trials temp_mean_membrane_potential_prior_onset
        
        
        %------------POST ONSET
        sp1 = subplot(4,1,1);
        hold on
        x = (analysis_membrane_potential_prior_onset+1:analysis_membrane_potential_prior_onset+analysis_membrane_potential_post_onset+1)/sampling_rate;
        clear y
        clear y2
        
        
        %Closeup - only mean line
        sp2 = subplot(4,1,2);
        hold on
        for z = 1:length(Analysis_Data_Flight(:, 49))
            if z2 == 1%medfilt
            y(:, z) = Analysis_Data_Flight{z, 49}(Analysis_Data_Flight{z, 44}:Analysis_Data_Flight{z, 44}+analysis_membrane_potential_post_onset);
            else %MEAN OF INTRA ORIG
            y(:, z) = Analysis_Data_Flight{z, 43}(Analysis_Data_Flight{z, 44}:Analysis_Data_Flight{z, 44}+analysis_membrane_potential_post_onset);
            end
            temp_mean_membrane_potential_post_onset(z) = mean(y(:, z), 1);
                        %*** plot Data of single Trials? NEW 2021_09_24
            y(:, z) = y(:, z)-median_onset(z);
%             plot(x, y(:,z), 'Color', [color(z, :) 0.3]) %NEW 2021_08_25
        end
        
        
        %Calculating Means (between trials and of all trials) and reduce plotting data
        mean_between_trials = mean(y, 2); %between
        mean_of_all_trials = mean(mean_between_trials, 1); %all trials
        
        
        %Calculating STD arround mean 
        for z = 1:length(y(:,1))
            std_y(z) = std(y(z,:));
        end
        std_y_downsample_x = linspace(x(1), x(end), 1000);
        std_y_downsample_y = interp1(x, std_y-mean_of_all_trials_onset, std_y_downsample_x);%HERE
        
        
        
        %Plotting Means
        xnew = linspace(x(1), x(end), 1000);
        ynew = interp1(x, mean_between_trials-mean_of_all_trials_onset, xnew);%HERE
        plot(xnew, ynew, 'LineWidth', 3, 'Color', 'r')
        
        %dotted line as mean
        y2(1:length(x)) = mean_of_all_trials-mean_of_all_trials_onset;%HERE
        plot(x,y2, '--', 'Color', 'r', 'LineWidth', 3);
        clear y2
        
        %plotting std arround mean
        patchx =  [std_y_downsample_x,fliplr(std_y_downsample_x)];
        patchy = [ynew + std_y_downsample_y,fliplr(ynew)-fliplr(std_y_downsample_y)];
        p = patch(patchx, patchy, 'red', 'LineStyle','none');
        alpha(0.2)
        
        %setting axis
        %         ylim([min(mean_between_trials)-2 max(mean_between_trials)+2 ])
        ylim([-15 15])
        
        %Save means and STD for writing it into txt file later on
        mean_membrane_potential_post_onset = mean(temp_mean_membrane_potential_post_onset, 2);
        mean_membrane_potential_post_onset_STD =  mean(std_y, 2);
        diff_in_Vm_onset = mean_membrane_potential_prior_onset-mean_membrane_potential_post_onset;%HERE
        mean_of_all_trials_post_onset = mean_of_all_trials;
        
        clearvars patchx patchy y2 xnew ynew std_y std_y_downsample_x std_y_downsample_y mean_between_trials mean_of_all_trials temp_mean_membrane_potential_post_onset
        
        
    end
    
    %Tacho
    sp3 = subplot(4,1,3);
    hold on
    x = (1:length(Analysis_Data_Flight{i, 47}(Analysis_Data_Flight{i, 44}-analysis_membrane_potential_prior_onset:Analysis_Data_Flight{i, 44}+analysis_membrane_potential_post_onset)))/sampling_rate;
    y = smooth(Analysis_Data_Flight{i, 47}(Analysis_Data_Flight{i, 44}-analysis_membrane_potential_prior_onset:Analysis_Data_Flight{i, 44}+analysis_membrane_potential_post_onset), smoothing_factor_tacho);
    %Setting beginning of Tacho Values to Zero
    y = y - y(1);
    plot(x,y, 'Color', color(i, :))
    
    %plot mean of Tacho
    if i == length(Analysis_Data_Flight(:, 47))
        

x = [analysis_membrane_potential_prior_onset/sampling_rate analysis_membrane_potential_prior_onset/sampling_rate];
y = [min(ylim) max(ylim)];
plot(x, y, '--', 'Color', 'k')

    end

end
%setting axes and enviroment of figur 
sp1 = subplot(4,1,1);
if plotting_figures_in_black == 1
    ax1 = gca;
    ax1.YLabel.String = 'Membrane Potential [mV]';
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    set(gca,'xtick',[])
else
    ax1 = gca;
    ax1.YLabel.String = 'Membrane Potential [mV]';
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    set(gca,'xtick',[])
end

sp2 = subplot(4,1,2);
if plotting_figures_in_black == 1
    ax2 = gca;
    ax2.YLabel.String = 'Membrane Potential [mV]';
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax2,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax2,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    set(gca,'xtick',[])
    title('Closeup - medfilt Vm');
else
    ax2 = gca;
    ax2.YLabel.String = 'Membrane Potential [mV]';
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax2,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax2,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    set(gca,'xtick',[])
    title('Closeup - medfilt Vm');
end

sp3 = subplot(4,1,3);
if plotting_figures_in_black == 1
    ax3 = gca;
    ax3.YLabel.String = 'Fly Activity';
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax3,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax3,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    %set(gca,'xtick',[])
else
    ax3 = gca;
    ax3.YLabel.String = 'Fly Activity';
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax3,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax3,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    %set(gca,'xtick',[])
end

ax3.XLabel.String = 'time (s)';

linkaxes([sp1, sp2, sp3], 'x'); %sp4
xlim([0 analysis_membrane_potential_prior_onset_in_sec+analysis_membrane_potential_post_onset_in_sec]);

if input_save_figures == 1
    if z2 == 1%medfilt
    filename = [file_modifier '_Membranepotential_onset_medfilt.tif'];
    saveas(gcf, filename)
    filename = [file_modifier '_Membranepotential_onset_medfilt.eps'];
    print(filename, '-depsc2', '-tiff', '-r300', '-painters')
    pause(1)
    else
    filename = [file_modifier '_Membranepotential_onset_meanfilt.tif'];
    saveas(gcf, filename)
    filename = [file_modifier '_Membranepotential_onset_meanfilt.eps'];
    print(filename, '-depsc2', '-tiff', '-r300', '-painters')
    pause(1)    
    end
end
% close gcf


%CLOSEUP FIGURE
    set(membrane_figure_onset, 'position', [1, 1, 500, 1500]);
    xlim([analysis_membrane_potential_prior_onset_in_sec-0.3 analysis_membrane_potential_prior_onset_in_sec+0.4])
    sp2 = subplot(4,1,2);
    ylim([-3 3]);
    if input_save_figures == 1
        if z2 == 1%medfilt
            filename = [file_modifier '_Membranepotential_onset_medfilt_CLOSEUP.tif'];
            saveas(gcf, filename)
            filename = [file_modifier '_Membranepotential_onset_medfilt_CLOSEUP.eps'];
            print(filename, '-depsc2', '-tiff', '-r300', '-painters')
            pause(1)
        else
            filename = [file_modifier '_Membranepotential_onset_meanfilt_CLOSEUP.tif'];
            saveas(gcf, filename)
            filename = [file_modifier '_Membranepotential_onset_meanfilt_CLOSEUP.eps'];
            print(filename, '-depsc2', '-tiff', '-r300', '-painters')
            pause(1)
        end
    end


clear sp1 sp2 sp3 sp4


 if z2 == 1
        fid = fopen([file_modifier '_Membranepotential_onset_LED_medfilt.txt'],'w');
        fprintf(fid, '%6s %12s\r\n', 'mean_of_all_trials_prior_onset_LED: ', num2str(mean_of_all_trials_prior_onset));
        fprintf(fid, '%6s %12s\r\n', 'mean_of_all_trials_post_onset_LED: ', num2str(mean_of_all_trials_post_onset));
        fprintf(fid, '%6s %12s\r\n', 'Dif_in_VM_ONSET_LED: ', num2str(diff_in_Vm_onset));%HERE
        fclose(fid);
    else
        fid = fopen([file_modifier '_Membranepotential_onset_LED_meanfilt.txt'],'w');
        fprintf(fid, '%6s %12s\r\n', 'mean_of_all_trials_prior_onset_LED: ', num2str(mean_of_all_trials_prior_onset));
        fprintf(fid, '%6s %12s\r\n', 'mean_of_all_trials_post_onset_LED: ', num2str(mean_of_all_trials_post_onset));
        fprintf(fid, '%6s %12s\r\n', 'Dif_in_VM_ONSET_LED: ', num2str(diff_in_Vm_onset));%HERE
        fclose(fid);
    end








% ********************** OFFSET
membrane_figure_offset = figure('name', 'Membrane Figure Offset', 'Color', figure_plotting_color);
set(membrane_figure_offset, 'position', [1, 1, 250, 1500]);

plotting_offset = 0; 

hold on
for i = 1:length(Analysis_Data_Flight(:, 45))
    sp1 = subplot(4,1,1);
    hold on
    %Prior Offset
    if any(selection_for_plotting(:)==i)
        x = (1:length(Analysis_Data_Flight{i, 43}(Analysis_Data_Flight{i, 45}-analysis_membrane_potential_prior_offset:Analysis_Data_Flight{i, 45})))/sampling_rate;
        y = Analysis_Data_Flight{i, 43}(Analysis_Data_Flight{i, 45}-analysis_membrane_potential_prior_offset:Analysis_Data_Flight{i, 45});
        %Adding plotting offset
        y_back_for_text = y;
        y = y + plotting_offset;
        plot(x,y, 'Color', color(i, :))
        
        %plot membrane pot of example
            xli = xlim;
            text(xli(1)+0.3, y(1),num2str(y_back_for_text(1)), 'HorizontalAlignment','right','VerticalAlignment','top', 'Color', 'k')%plot membrane pot
        
        %Post Offset
        x = (analysis_membrane_potential_prior_offset+1:analysis_membrane_potential_prior_offset+analysis_membrane_potential_post_offset+1)/sampling_rate;
        y = Analysis_Data_Flight{i, 43}(Analysis_Data_Flight{i, 45}:Analysis_Data_Flight{i, 45}+analysis_membrane_potential_post_offset);
        %Adding plotting offset
        y_back_for_text = y;
        y = y + plotting_offset;
        plotting_offset = plotting_offset + 30;
        plot(x,y, 'Color', color(i, :))
        
        %plot membrane pot of example
            xli = xlim;
            text(xli(2)+0.3, y(end),num2str(y_back_for_text(end)), 'HorizontalAlignment','right','VerticalAlignment','top', 'Color', 'k')%plot membrane pot
    
    end
    
    %plot mean of intra
    %-----------Prior OFFSET
    if i == length(Analysis_Data_Flight(:, 45))
        x = (1:length(mean(Analysis_Data_Flight{i, 43}(Analysis_Data_Flight{i, 45}-analysis_membrane_potential_prior_offset:Analysis_Data_Flight{i, 45}),2)))/sampling_rate;
        clear y
        clear y2
        
        % Mean Closeup of medfitted Curve with STD as Patch
        sp2 = subplot(4,1,2);
        hold on
        for z = 1:length(Analysis_Data_Flight(:, 49))
            if z2 == 1%medfilt
            y(:, z) = Analysis_Data_Flight{z, 49}(Analysis_Data_Flight{z, 45}-analysis_membrane_potential_prior_offset:Analysis_Data_Flight{z, 45});
            else %MEAN OF INTRA ORIG
            y(:, z) = Analysis_Data_Flight{z, 43}(Analysis_Data_Flight{z, 45}-analysis_membrane_potential_prior_offset:Analysis_Data_Flight{z, 45});
            end
            temp_mean_membrane_potential_prior_offset(z) = mean(y(:, z), 1);
            median_offset(z) = median(y(end-0.5*sampling_rate:end, z));
              %*** plot Data of single Trials? NEW 2021_09_24
            y(:, z) = y(:, z)-median_offset(z);
%             plot(x, y(:,z), 'Color', [color(z, :) 0.3]) %NEW 2021_08_25
        end
        
        
        %Calculating Means (between trials and of all trials) and reduce plotting data
        mean_between_trials = mean(y, 2); %between
        mean_of_all_trials = mean(mean_between_trials, 1); %all trials
        mean_of_all_trials_offset = mean_of_all_trials;%HERE
        
        
        %Calculating STD arround mean 
        clear std_y
        for z = 1:length(y(:,1))
            std_y(z) = std(y(z,:));
        end
        std_y_downsample_x = linspace(x(1), x(end), 1000);
        std_y_downsample_y = interp1(x, std_y-mean_of_all_trials_offset, std_y_downsample_x);%HERE
        
        
        
        %Plotting Means
        xnew = linspace(x(1), x(end), 1000);
        ynew = interp1(x, mean_between_trials-mean_of_all_trials_offset, xnew);%HERE
        plot(xnew, ynew, 'LineWidth', 3, 'Color', 'g')
        
        %dotted line as mean
        y2(1:length(x)) = mean_of_all_trials-mean_of_all_trials_offset;%HERE
        plot(x,y2, '--', 'Color', 'g', 'LineWidth', 3);
        clear y2
        
        %plotting std arround mean
        patchx =  [std_y_downsample_x,fliplr(std_y_downsample_x)];
        patchy = [ynew + std_y_downsample_y,fliplr(ynew)-fliplr(std_y_downsample_y)];
        p = patch(patchx, patchy, 'red', 'LineStyle','none');
        alpha(0.2)
        
        %setting axis
        ylim([min(mean_between_trials)-2 max(mean_between_trials)+2 ])
        
        %Save means and STD for writing it into txt file later on
        mean_membrane_potential_prior_offset = mean(temp_mean_membrane_potential_prior_offset, 2);
        mean_membrane_potential_prior_offset_STD =  mean(std_y, 2);
        
        mean_of_all_trials_prior_offset = mean_of_all_trials;
        
        clearvars patchx patchy y2 xnew ynew std_y_downsample_x std_y_downsample_y mean_between_trials mean_of_all_trials temp_mean_membrane_potential_prior_offset
        
        
        
        
        %------------POST OFFSET
        sp1 = subplot(4,1,1);
        clear y
        x = (analysis_membrane_potential_prior_offset+1:analysis_membrane_potential_prior_offset+analysis_membrane_potential_post_offset+1)/sampling_rate;
        
        %Mean Closeup
        sp2 = subplot(4,1,2);
        hold on
        for z = 1:length(Analysis_Data_Flight(:, 49))
            if z2 ==1 %medfilt
            y(:, z) = Analysis_Data_Flight{z, 49}(Analysis_Data_Flight{z, 45}:Analysis_Data_Flight{z, 45}+analysis_membrane_potential_post_offset);
            else
            y(:, z) = Analysis_Data_Flight{z, 43}(Analysis_Data_Flight{z, 45}:Analysis_Data_Flight{z, 45}+analysis_membrane_potential_post_offset);
            end
            temp_mean_membrane_potential_post_offset(z) = mean(y(:, z), 1);
              %*** plot Data of single Trials? NEW 2021_09_24
            y(:, z) = y(:, z)-median_offset(z);
%             plot(x, y(:,z), 'Color', [color(z, :) 0.3]) %NEW 2021_08_25
        end
        
        
        %Calculating Means (between trials and of all trials) and reduce plotting data
        mean_between_trials = mean(y, 2); %between
        mean_of_all_trials = mean(mean_between_trials, 1); %all trials
        
        
        %Calculating STD arround mean 
        for z = 1:length(y(:,1))
            std_y(z) = std(y(z,:));
        end
        std_y_downsample_x = linspace(x(1), x(end), 1000);
        std_y_downsample_y = interp1(x, std_y-mean_of_all_trials_offset, std_y_downsample_x);%HERE
        
        
        
        %Plotting Means
        xnew = linspace(x(1), x(end), 1000);
        ynew = interp1(x, mean_between_trials-mean_of_all_trials_offset, xnew);%HERE
        plot(xnew, ynew, 'LineWidth', 3, 'Color', 'r')
        
        %dotted line as mean
        y2(1:length(x)) = mean_of_all_trials-mean_of_all_trials_offset;%HERE
        plot(x,y2, '--', 'Color', 'r', 'LineWidth', 3);
        clear y2
        
        %plotting std arround mean
        patchx =  [std_y_downsample_x,fliplr(std_y_downsample_x)];
        patchy = [ynew + std_y_downsample_y,fliplr(ynew)-fliplr(std_y_downsample_y)];
        p = patch(patchx, patchy, 'red', 'LineStyle','none');
        alpha(0.2)
        
        %setting axis
        %         ylim([min(mean_between_trials)-2 max(mean_between_trials)+2 ])
        ylim([-15 15])
        
        %Save means and STD for writing it into txt file later on
        mean_membrane_potential_post_offset = mean(temp_mean_membrane_potential_post_offset, 2);
        mean_membrane_potential_post_offset_STD =  mean(std_y, 2);
        diff_in_Vm_offset = mean_membrane_potential_prior_offset-mean_membrane_potential_post_offset;%HERE
        mean_of_all_trials_post_offset = mean_of_all_trials;
        
        clearvars patchx patchy y2 xnew ynew std_y_downsample_x std_y_downsample_y mean_between_trials mean_of_all_trials temp_mean_membrane_potential_post_offset
        
    end
    
    %Tacho
    sp3 = subplot(4,1,3);
    hold on
    x = (1:length(Analysis_Data_Flight{i, 47}(Analysis_Data_Flight{i, 45}-analysis_membrane_potential_prior_offset:Analysis_Data_Flight{i, 45}+analysis_membrane_potential_post_offset)))/sampling_rate;
    y = smooth(Analysis_Data_Flight{i, 47}(Analysis_Data_Flight{i, 45}-analysis_membrane_potential_prior_offset:Analysis_Data_Flight{i, 45}+analysis_membrane_potential_post_offset), smoothing_factor_tacho);
    %Setting beginning of Tacho Values to Zero
    y = y - y(end);
    plot(x,y, 'Color', color(i, :))
    
    if i == length(Analysis_Data_Flight(:, 47))

x = [analysis_membrane_potential_prior_offset/sampling_rate analysis_membrane_potential_prior_offset/sampling_rate];
y = [min(ylim) max(ylim)];
plot(x, y, '--', 'Color', 'r')

    end
    

end

sp1 = subplot(4,1,1);
if plotting_figures_in_black == 1
    ax1 = gca;
    ax1.YLabel.String = 'Membrane Potential [mV]';
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    set(gca,'xtick',[])
else
    ax1 = gca;
    ax1.YLabel.String = 'Membrane Potential [mV]';
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    set(gca,'xtick',[])
end

sp2 = subplot(4,1,2);
if plotting_figures_in_black == 1
    ax2 = gca;
    ax2.YLabel.String = 'Membrane Potential [mV]';
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax2,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax2,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    set(gca,'xtick',[])
    title('Closeup - medfilt Vm');
else
    ax2 = gca;
    ax2.YLabel.String = 'Membrane Potential [mV]';
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax2,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax2,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    set(gca,'xtick',[])
    title('Closeup - medfilt Vm');
end

sp3 = subplot(4,1,3);
if plotting_figures_in_black == 1
    ax3 = gca;
    ax3.YLabel.String = 'Fly Activity';
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax3,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax3,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    %set(gca,'xtick',[])
else
    ax3 = gca;
    ax3.YLabel.String = 'Fly Activity';
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax3,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax3,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    a = get(gca,'XTickLabel');
    set(gca,'FontName','Helvetica','fontsize',12);
    %set(gca,'xtick',[])
end

ax3.XLabel.String = 'time (s)';

linkaxes([sp1, sp2, sp3], 'x'); %sp4
xlim([0 analysis_membrane_potential_prior_offset_in_sec+analysis_membrane_potential_post_offset_in_sec]);

if input_save_figures == 1
           if z2 == 1
    filename = [file_modifier '_Membranepotential_offset_medfilt.tif'];
    saveas(gcf, filename)
    filename = [file_modifier '__Membranepotential_offset_medfilt.eps'];
    print(filename, '-depsc2', '-tiff', '-r300', '-painters')
    pause(1)
           else
               filename = [file_modifier '_Membranepotential_offset_meanfilt.tif'];
    saveas(gcf, filename)
    filename = [file_modifier '__Membranepotential_offset_meanfilt.eps'];
    print(filename, '-depsc2', '-tiff', '-r300', '-painters')
    pause(1)    
           end
end

% close gcf

 %CLOSEUP FIGURE
    set(membrane_figure_offset, 'position', [1, 1, 500, 1500]);
     xlim([analysis_membrane_potential_prior_onset_in_sec-0.3 analysis_membrane_potential_prior_onset_in_sec+0.4])
    sp2 = subplot(4,1,2);
    ylim([-3 3]);
    if input_save_figures == 1
        if z2 == 1%medfilt
            filename = [file_modifier '_Membranepotential_offset_medfilt_CLOSEUP.tif'];
            saveas(gcf, filename)
            filename = [file_modifier '_Membranepotential_offset_medfilt_CLOSEUP.eps'];
            print(filename, '-depsc2', '-tiff', '-r300', '-painters')
            pause(1)
        else
            filename = [file_modifier '_Membranepotential_offset_meanfilt_CLOSEUP.tif'];
            saveas(gcf, filename)
            filename = [file_modifier '_Membranepotential_offset_meanfilt_CLOSEUP.eps'];
            print(filename, '-depsc2', '-tiff', '-r300', '-painters')
            pause(1)
        end
    end
    
    
clear sp1 sp2 sp3 sp4


if z2 == 1
        fid = fopen([file_modifier '_Membranepotential_offset_LED_medfilt.txt'],'w');
        fprintf(fid, '%6s %12s\r\n', 'mean_of_all_trials_prior_offset: ', num2str(mean_of_all_trials_prior_offset));
        fprintf(fid, '%6s %12s\r\n', 'mean_of_all_trials_post_offset: ', num2str(mean_of_all_trials_post_offset));
        fprintf(fid, '%6s %12s\r\n', 'Dif_in_VM_OffSET_LED: ', num2str(diff_in_Vm_offset));%HERE
        fclose(fid);
    else
        fid = fopen([file_modifier '_Membranepotential_offset_LED_meanfilt.txt'],'w');
        fprintf(fid, '%6s %12s\r\n', 'mean_of_all_trials_prior_offset: ', num2str(mean_of_all_trials_prior_offset));
        fprintf(fid, '%6s %12s\r\n', 'mean_of_all_trials_post_offset: ', num2str(mean_of_all_trials_post_offset));
        fprintf(fid, '%6s %12s\r\n', 'Dif_in_VM_OffSET_LED: ', num2str(diff_in_Vm_offset));%HERE
        fclose(fid);
    end

end

%Save Mean Values and STD as TXT
fid = fopen('Membranepotential_Figure_onset_and_offset.txt','w');
fprintf(fid, '%6s %12s\r\n', 'mean_membrane_potential_prior_onset: ', num2str(mean_membrane_potential_prior_onset));
fprintf(fid, '%6s %12s\r\n', 'STD_membrane_potential_prior_onset: ', num2str(mean_membrane_potential_prior_onset_STD));
fprintf(fid, '%6s %12s\r\n', 'mean_membrane_potential_post_onset: ', num2str(mean_membrane_potential_post_onset));
fprintf(fid, '%6s %12s\r\n', 'STD_membrane_potential_post_onset: ', num2str(mean_membrane_potential_post_onset_STD));
fprintf(fid, '%6s %12s\r\n', 'mean_membrane_potential_prior_offset: ', num2str(mean_membrane_potential_prior_offset));
fprintf(fid, '%6s %12s\r\n', 'STD_membrane_potential_prior_offset: ', num2str(mean_membrane_potential_prior_offset_STD));
fprintf(fid, '%6s %12s\r\n', 'mean_membrane_potential_post_offset: ', num2str(mean_membrane_potential_post_offset));
fprintf(fid, '%6s %12s\r\n', 'STD_membrane_potential_post_offset: ', num2str(mean_membrane_potential_post_offset_STD));
fclose(fid);


%% Determine Number Flies used (works only with flys below 100 and cells
%%below 100
number_of_flies = 1;
number_of_cells = 1;
for i = 2:length(spike_vector_per_flyID(:,1))
    
    number_of_flies = number_of_flies + 1;
    number_of_cells = number_of_cells + 1;
    %Find position of _Cell which can be used to count cells and flies
    k = strfind(spike_vector_per_flyID{i,1},'_CELL' );
    if isequal(spike_vector_per_flyID{i,1}(1:k-1), spike_vector_per_flyID{i-1,1}(1:k-1))
        number_of_flies = number_of_flies - 1; %23
        
        
        if isequal(spike_vector_per_flyID{i,1}(k:end-6), spike_vector_per_flyID{i-1,1}(k:end-6))
            number_of_cells = number_of_cells - 1; %~+10
        end
        
    end
end
disp(['Number of Flies analyzed: ' num2str(number_of_flies) '; Number of Cells analyzed: ' num2str(number_of_cells)]);



%% Keep individual Plots?
% Promt to ask whether all Figures shall be closd before showing the final
% results?
dlgTitle    = 'Close all Figures?';
dlgQuestion = 'Do you wish to Close all Figure before showing the summarizing results?';
choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
if length(choice)==3  %close all
    close all;
end

disp('Plotting Results. Depending in the amount of data, plots might appear minutes after "finish"...');

%% Save Summarizing Plots?
dlgTitle    = 'Save summarizing Results?';
dlgQuestion = 'Do you wish to save all upcoming result figures?';
save_choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
if length(save_choice)==3  %close all
    input_save_figures = 1;
else
    input_save_figures = 0;
end

%% ************* Plot All Spike Traces  *************
%For Onset
% Variables needed
total_min = [];
total_max = [];
bout_counter = 0;
if checkbox_show_all_results_onset == 1
    h = figure('name',['All Flight Bouts from ' Analysis_Data_Flight{1,1} ' to ' Analysis_Data_Flight{end,1}]);
    th = tiledlayout(length(Analysis_Data_Flight(:,1)),1);
    for u = 1:length(Analysis_Data_Flight(:,1))
        nexttile
        %Determine max Y Value in the files to plot for yLim
        if u == bout_counter+1
            total_max = max(Analysis_Data_Flight{u,4});
            total_min = min(Analysis_Data_Flight{u,4});
            for  z = bout_counter+1:length(Analysis_Data_Flight(:,1))
                current_min = min(Analysis_Data_Flight{z,4});
                current_max = max(Analysis_Data_Flight{z,4});
                if current_min < total_min ==1
                    total_min = current_min;
                end
                
                if current_max > total_max ==1
                    total_max = current_max;
                end
            end
        end
        %Plot Intra Trace
        x_axis = (1:length(Analysis_Data_Flight{u,4}(1:end)))/sampling_rate;
        plot(x_axis, Analysis_Data_Flight{u,4}(1:end))
        hold on
        %Plot Start of flight
        plot(analysis_window_pre_onset./sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
        %Plot end of flight if in range
        %         if Analysis_Data_Flight{u,10}-Analysis_Data_Flight{u,9}<analysis_window_pre_onset
        if Analysis_Data_Flight{u,10}-Analysis_Data_Flight{u,9}+analysis_window_pre_onset<analysis_window_post_offset
            %             plot((Analysis_Data_Flight{u,10}-Analysis_Data_Flight{u,9}+analysis_window_pre_onset)/sampling_rate, min(Analysis_Data_Flight{u,4}):max(Analysis_Data_Flight{u,4}), '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1], 'LineWidth', 2);
            %             plot((Analysis_Data_Flight{u,10}-Analysis_Data_Flight{u,9}+analysis_window_pre_onset)/sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1], 'LineWidth', 2);
            plot((Analysis_Data_Flight{u,10}./sampling_rate)-(Analysis_Data_Flight{u,9}./sampling_rate)+analysis_window_pre_onset./sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1], 'LineWidth', 2);
            
        end
        ax = gca;
        if u==length(Analysis_Data_Flight)
            set(ax,'YTick', []);
        else
            set(ax,'XTick',[], 'YTick', []);
        end
    end
    ylim([total_min, total_max]);
    th.Padding = 'none';
    th.TileSpacing = 'none';
    hold off
    %     if input_save_figures == 1
    %         pause(0.5)
    %     saveas(gcf,'Figure 6.tif')
    %     saveas(gcf,'Figure 6','epsc')
    %     end
    
    %For Offset
    % Variables needed
    total_min = [];
    total_max = [];
    current_max = [];
    current_min = [];
    h2 = figure('name',['All Flight Bouts from ' Analysis_Data_Flight{1,1} ' to ' Analysis_Data_Flight{end,1}]);
    th2 = tiledlayout(length(Analysis_Data_Flight(:,1)),1);
    for u = 1:length(Analysis_Data_Flight(:,1))
        nexttile
        %Determine max Y Value in the files to plot for yLim
        if u == bout_counter+1
            total_max = max(Analysis_Data_Flight{u,8});
            total_min = min(Analysis_Data_Flight{u,8});
            for  z = bout_counter+1:length(Analysis_Data_Flight(:,1))
                current_min = min(Analysis_Data_Flight{z,8});
                current_max = max(Analysis_Data_Flight{z,8});
                if current_min < total_min ==1
                    total_min = current_min;
                end
                
                if current_max > total_max ==1
                    total_max = current_max;
                end
            end
        end
        %Plot Intra Trace
        x_axis = (1:length(Analysis_Data_Flight{u,8}(1:end)))/sampling_rate;
        plot(x_axis, Analysis_Data_Flight{u,8}(1:end))
        hold on
        %Plot End of flight
        plot(analysis_window_pre_onset./sampling_rate,total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[1 0 1], 'LineWidth', 2);
        %Plot Start of this flight if in range
        if Analysis_Data_Flight{u,10}-Analysis_Data_Flight{u,9}<analysis_window_pre_onset
            %             plot((Analysis_Data_Flight{u,9}-Analysis_Data_Flight{u,10}+analysis_window_pre_onset)/sampling_rate, min(Analysis_Data_Flight{u,8}):max(Analysis_Data_Flight{u,8}), '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
            
            plot((Analysis_Data_Flight{u,9}-Analysis_Data_Flight{u,10}+analysis_window_pre_onset)/sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
        end
        % % %  DOESNT WORK YET        %Is there another flight start of the next bout in range (analysis_window_post_offset)
        % % %    if u ~= length(Analysis_Data_Flight(:,1)) && isequal(Analysis_Data_Flight{u+1,1}, Analysis_Data_Flight{u,1}) &&  Analysis_Data_Flight{u+1,9} - Analysis_Data_Flight{u,10} < analysis_window_post_offset
        % % %             %            plot((Analysis_Data_Flight{u+1,9}-Analysis_Data_Flight{u+1,10}+analysis_window_pre_onset)/sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
        % % %             plot((Analysis_Data_Flight{u,10}+Analysis_Data_Flight{u+1,9})/sampling_rate, total_min:total_max, '|', 'MarkerSize',5, 'MarkerEdgeColor',[0 0 1], 'LineWidth', 2);
        % % %         end
        
        ax = gca;
        if u==length(Analysis_Data_Flight)
            set(ax,'YTick', []);
        else
            set(ax,'XTick',[], 'YTick', []);
        end
    end
    ylim([total_min, total_max]);
    th2.Padding = 'none';
    th2.TileSpacing = 'none';
    hold off
    
    %     if input_save_figures == 1
    %         pause(0.5)
    %     saveas(gcf,'Figure 7.tif')
    %     saveas(gcf,'Figure 7','epsc')
    % end
    
end





%% ********** Ask whether number of spikes shall be determined per
%%%%%%*.Mat File or per Cell ID
dlgTitle    = 'Analyze number of Spikes per *.MAT File or per Cell-ID?';
dlgQuestion = 'Analyze number of Spikes per *.MAT File or per Cell-ID?';
choice = questdlg(dlgQuestion,dlgTitle,'Cell-ID','*.MAT', 'Cell-ID');
if length(choice)==7
    analysed_dataset = 'per cell';
    idx_counter = 0;
    clearvars spike_vector_per_flyID_temp
    for i = 1:length(spike_vector_per_flyID(:, 1))
        if i ~= 1 %Cause it can not be compared in the first cell
            
            if isequal(spike_vector_per_flyID{i,1}(1:end-6), spike_vector_per_flyID{i-1,1}(1:end-6)) == 0 %If it is not the same fly
                spike_vector_per_flyID_temp{i-idx_counter, 1} = spike_vector_per_flyID{i, 1}; % If it is not the same cell than use these values
                spike_vector_per_flyID_temp{i-idx_counter, 2} = spike_vector_per_flyID{i, 2};
                spike_vector_per_flyID_temp{i-idx_counter, 3} = spike_vector_per_flyID{i, 3};
                spike_vector_per_flyID_temp{i-idx_counter, 4} = spike_vector_per_flyID{i, 4};
                spike_vector_per_flyID_temp{i-idx_counter, 5} = spike_vector_per_flyID{i, 5};
                spike_vector_per_flyID_temp{i-idx_counter, 6} = spike_vector_per_flyID{i, 6};
                spike_vector_per_flyID_temp{i-idx_counter, 7} = spike_vector_per_flyID{i, 7};
                spike_vector_per_flyID_temp{i-idx_counter, 8} = spike_vector_per_flyID{i, 8};
                spike_vector_per_flyID_temp{i-idx_counter, 9} = spike_vector_per_flyID{i, 9};
                %         spike_vector_per_flyID_temp{i-idx_counter, 10} = spike_vector_per_flyID{i, 10};
                %         spike_vector_per_flyID_temp{i-idx_counter, 11} = spike_vector_per_flyID{i, 11};
                
                %         idx_counter = idx_counter - 1;
            else %If it it the same cell but has more than one trial ('Mat files') than calculate the mean from all entries
                
                spike_vector_per_flyID_temp{i-1-idx_counter, 1} = spike_vector_per_flyID{i-1, 1};
                spike_vector_per_flyID_temp{i-1-idx_counter, 2} = [spike_vector_per_flyID_temp{i-1-idx_counter, 2} ; spike_vector_per_flyID{i, 2}]; %NEW 2021_05_26 there was an error here + following lines in indices, fixed it. adapt to all other scripts!
                spike_vector_per_flyID_temp{i-1-idx_counter, 3} = [spike_vector_per_flyID_temp{i-1-idx_counter, 3} ; spike_vector_per_flyID{i, 3}];
                spike_vector_per_flyID_temp{i-1-idx_counter, 4} = [spike_vector_per_flyID_temp{i-1-idx_counter, 4} ; spike_vector_per_flyID{i, 4}];
                spike_vector_per_flyID_temp{i-1-idx_counter, 5} = [spike_vector_per_flyID_temp{i-1-idx_counter, 5} ; spike_vector_per_flyID{i, 5}];
                spike_vector_per_flyID_temp{i-1-idx_counter, 6} = [spike_vector_per_flyID_temp{i-1-idx_counter, 6} ; spike_vector_per_flyID{i, 6}];
                spike_vector_per_flyID_temp{i-1-idx_counter, 7} = [spike_vector_per_flyID_temp{i-1-idx_counter, 7} ; spike_vector_per_flyID{i, 7}];
                spike_vector_per_flyID_temp{i-1-idx_counter, 8} = [spike_vector_per_flyID_temp{i-1-idx_counter, 8} ; spike_vector_per_flyID{i, 8}];
                spike_vector_per_flyID_temp{i-1-idx_counter, 9} = [spike_vector_per_flyID_temp{i-1-idx_counter, 9} ; spike_vector_per_flyID{i, 9}];
                %        spike_vector_per_flyID_temp{i-idx_counter, 10} = [spike_vector_per_flyID_temp{i-1, 10} ; spike_vector_per_flyID{i-1-idx_counter, 10}];
                %        spike_vector_per_flyID_temp{i-idx_counter, 11} = [spike_vector_per_flyID_temp{i-1, 11} ; spike_vector_per_flyID{i-1-idx_counter, 11}];
                idx_counter = idx_counter + 1;
                
                
            end
            
        else %for the first cell
            spike_vector_per_flyID_temp{i-idx_counter, 1} = spike_vector_per_flyID{i, 1}; % If it is not the same cell than use these values
            spike_vector_per_flyID_temp{i-idx_counter, 2} = spike_vector_per_flyID{i-idx_counter, 2};
            spike_vector_per_flyID_temp{i-idx_counter, 3} = spike_vector_per_flyID{i-idx_counter, 3};
            spike_vector_per_flyID_temp{i-idx_counter, 4} = spike_vector_per_flyID{i-idx_counter, 4};
            spike_vector_per_flyID_temp{i-idx_counter, 5} = spike_vector_per_flyID{i-idx_counter, 5};
            spike_vector_per_flyID_temp{i-idx_counter, 6} = spike_vector_per_flyID{i-idx_counter, 6};
            spike_vector_per_flyID_temp{i-idx_counter, 7} = spike_vector_per_flyID{i-idx_counter, 7};
            spike_vector_per_flyID_temp{i-idx_counter, 8} = spike_vector_per_flyID{i-idx_counter, 8};
            spike_vector_per_flyID_temp{i-idx_counter, 9} = spike_vector_per_flyID{i-idx_counter, 9};
            %         spike_vector_per_flyID_temp{i-idx_counter, 10} = spike_vector_per_flyID{i-idx_counter, 10};
            %         spike_vector_per_flyID_temp{i-idx_counter, 11} = spike_vector_per_flyID{i-idx_counter, 11};
        end
        
    end
    clearvars spike_vector_per_flyID
    spike_vector_per_flyID = spike_vector_per_flyID_temp;
else
    analysed_dataset = 'per Mat';
end











%% ************* Plot all Spikes in spike diagram + Calculating Number of Spikes (unbinnend) within that timeframe
% % spike_vector_per_flyID{idx_fly_counter,2} = temp_pre_onset;
% % spike_vector_per_flyID{idx_fly_counter,3} = temp_post_onset;
% % spike_vector_per_flyID{idx_fly_counter,4} = temp_pre_offset;
% % spike_vector_per_flyID{idx_fly_counter,5} = temp_post_offset;

% ********* NEW PLOTTING SPIKE VECTOR SORTED BY NUMBER OF SPIKES
dlgTitle    = 'Sort Flies in ascending order by number of spikes?';
dlgQuestion = 'Sort Flies in ascending order by number of spikes?';
choice = questdlg(dlgQuestion,dlgTitle,'YES','NO', 'NO');
if length(choice)==3
    dlgTitle    = 'Sort Flies in ascending order by number of spikes after ONSET or OFFSET?';
    dlgQuestion = 'Sort Flies in ascending order by number of spikes after ONSET or OFFSET?';
    choice = questdlg(dlgQuestion,dlgTitle,'ONSET','OFFSET', 'OFFSET');
    if length(choice)==6%=offset
        for i = 1:length(spike_vector_per_flyID(:,5))
            spike_vector_per_flyID{i,11} = length(find(spike_vector_per_flyID{i,5}(:,:)==1))/(analysis_window_post_onset/sampling_rate); %NUMBER OF SPIKES POST ONSET PER MAT
        end
        sorted_by = 'sorted by offset';
    else
        for i = 1:length(spike_vector_per_flyID(:,3))
            spike_vector_per_flyID{i,11} = length(find(spike_vector_per_flyID{i,3}(:,:)==1))/(analysis_window_post_offset/sampling_rate); %NUMBER OF SPIKES POST ONSET PER MAT
        end
        sorted_by = 'sorted by onset';
    end
    % If sorted by onset --> sorted in 'ascend' order
    if length(choice)==6%=offset
        spike_vector_per_flyID_sorted_by_Spikes = sortrows(spike_vector_per_flyID, 11, 'ascend');
    else
        spike_vector_per_flyID_sorted_by_Spikes = sortrows(spike_vector_per_flyID, 11, 'descend');
    end
    spike_vector_per_flyID_back = spike_vector_per_flyID;
    spike_vector_per_flyID = spike_vector_per_flyID_sorted_by_Spikes;
else
    spike_vector_per_flyID_back = spike_vector_per_flyID;
end

clearvars spike_vector_per_flyID_back spike_vector_per_flyID_sorted_by_Spikes spike_vector_per_flyID_temp


% Calculating Spike Frequency before plotting
%*****************PRIOR ONSET
clearvars pre_onset_relative_spike_times temp_pre_onset_spike_numbers temp_pre_onset_spike_frequency average_number_of_spikes_pre_onset spike_frequency_per_mat_pre_onset avaraged_number_of_spikes_per_mat_pre_onset all_average_number_of_spikes_pre_onset

for i=1:length(spike_vector_per_flyID(:,1)) %for all Cells
    
    for z = 1:length(spike_vector_per_flyID{i,2}(:,1))%for all bouts recorded in that cell
        %Pre Onset (spike_vector_per_flyID{:,2}), analysis window is: analysis_window_pre_onset
        
        
        %Finding Spike times in that vector
        pre_onset_relative_spike_times{i,z} = find(spike_vector_per_flyID{i,2}(z,1:end)==1)./sampling_rate;
        %Extracting number of spikes wihtin analysis window
        temp_pre_onset_spike_numbers(z) = length(pre_onset_relative_spike_times{i,z});
        temp_pre_onset_spike_frequency(z) = temp_pre_onset_spike_numbers(z)/(analysis_window_pre_onset/sampling_rate);
        
        %At the end of each recording, avarage the spike frequency over all
        %bouts recorded in that mat file
        if z == length(spike_vector_per_flyID{i,2}(:,1))
            spike_frequency_per_mat_pre_onset(i) = mean(temp_pre_onset_spike_frequency);
            avaraged_number_of_spikes_per_mat_pre_onset(i) = mean(temp_pre_onset_spike_numbers);
            clearvars temp_pre_onset_spike_numbers(z) temp_pre_onset_spike_frequency
        end
    end
    %At the end of all recorded mat files determine the average spike frequency and
    %numer of spikes for all mat files (pre_onset)
    if i == length(spike_vector_per_flyID(:,1))
        all_average_spike_frequency_pre_onset =  mean(spike_frequency_per_mat_pre_onset);
        all_average_number_of_spikes_pre_onset =  mean(avaraged_number_of_spikes_per_mat_pre_onset);
    end
    
    
end

%*****************POST ONSET
clearvars post_onset_relative_spike_times temp_post_onset_spike_numbers temp_post_onset_spike_frequency average_number_of_spikes_post_onset spike_frequency_per_mat_post_onset avaraged_number_of_spikes_per_mat_post_onset all_average_number_of_spikes_post_onset

for i=1:length(spike_vector_per_flyID(:,1)) %for all Cells
    
    for z = 1:length(spike_vector_per_flyID{i,3}(:,1))%for all bouts recorded in that cell
        %post Onset (spike_vector_per_flyID{:,2}), analysis window is: analysis_window_post_onset
        
        %Finding Spike times in that vector
        post_onset_relative_spike_times{i,z} = find(spike_vector_per_flyID{i,3}(z,1:end)==1)./sampling_rate;
        %Extracting number of spikes wihtin analysis window
        temp_post_onset_spike_numbers(z) = length(post_onset_relative_spike_times{i,z});
        temp_post_onset_spike_frequency(z) = temp_post_onset_spike_numbers(z)/(analysis_window_post_onset/sampling_rate);
        
        %At the end of each recording, avarage the spike frequency over all
        %bouts recorded in that mat file
        if z == length(spike_vector_per_flyID{i,3}(:,1))
            spike_frequency_per_mat_post_onset(i) = mean(temp_post_onset_spike_frequency);
            avaraged_number_of_spikes_per_mat_post_onset(i) = mean(temp_post_onset_spike_numbers);
            clearvars temp_post_onset_spike_numbers(z) temp_post_onset_spike_frequency
        end
    end
    %At the end of all recorded mat files determine the average spike frequency and
    %numer of spikes for all mat files (post_onset)
    if i == length(spike_vector_per_flyID(:,1))
        all_average_spike_frequency_post_onset =  mean(spike_frequency_per_mat_post_onset);
        all_average_number_of_spikes_post_onset =  mean(avaraged_number_of_spikes_per_mat_post_onset);
    end
    
    
end


%*****************PRIOR OFFSET
clearvars pre_offset_relative_spike_times temp_pre_offset_spike_numbers temp_pre_offset_spike_frequency average_number_of_spikes_pre_offset spike_frequency_per_mat_pre_offset avaraged_number_of_spikes_per_mat_pre_offset all_average_number_of_spikes_pre_offset

for i=1:length(spike_vector_per_flyID(:,1)) %for all Cells
    
    for z = 1:length(spike_vector_per_flyID{i,4}(:,1))%for all bouts recorded in that cell
        %Pre offset (spike_vector_per_flyID{:,2}), analysis window is: analysis_window_pre_offset
        
        %Finding Spike times in that vector
        pre_offset_relative_spike_times{i,z} = find(spike_vector_per_flyID{i,4}(z,1:end)==1)./sampling_rate;
        %Extracting number of spikes wihtin analysis window
        temp_pre_offset_spike_numbers(z) = length(pre_offset_relative_spike_times{i,z});
        temp_pre_offset_spike_frequency(z) = temp_pre_offset_spike_numbers(z)/(analysis_window_pre_offset/sampling_rate);
        
        %At the end of each recording, avarage the spike frequency over all
        %bouts recorded in that mat file
        if z == length(spike_vector_per_flyID{i,4}(:,1))
            spike_frequency_per_mat_pre_offset(i) = mean(temp_pre_offset_spike_frequency);
            avaraged_number_of_spikes_per_mat_pre_offset(i) = mean(temp_pre_offset_spike_numbers);
            clearvars temp_pre_offset_spike_numbers(z) temp_pre_offset_spike_frequency
        end
    end
    %At the end of all recorded mat files determine the average spike frequency and
    %numer of spikes for all mat files (pre_offset)
    if i == length(spike_vector_per_flyID(:,1))
        all_average_spike_frequency_pre_offset =  mean(spike_frequency_per_mat_pre_offset);
        all_average_number_of_spikes_pre_offset =  mean(avaraged_number_of_spikes_per_mat_pre_offset);
    end
    
    
end

%*****************POST OFFSET
clearvars post_offset_relative_spike_times temp_post_offset_spike_numbers temp_post_offset_spike_frequency average_number_of_spikes_post_offset spike_frequency_per_mat_post_offset avaraged_number_of_spikes_per_mat_post_offset all_average_number_of_spikes_post_offset

for i=1:length(spike_vector_per_flyID(:,1)) %for all Cells
    
    for z = 1:length(spike_vector_per_flyID{i,5}(:,1))%for all bouts recorded in that cell
        %post offset (spike_vector_per_flyID{:,2}), analysis window is: analysis_window_post_offset
        
        %Finding Spike times in that vector
        post_offset_relative_spike_times{i,z} = find(spike_vector_per_flyID{i,5}(z,1:end)==1)./sampling_rate;
        %Extracting number of spikes wihtin analysis window
        temp_post_offset_spike_numbers(z) = length(post_offset_relative_spike_times{i,z});
        temp_post_offset_spike_frequency(z) = temp_post_offset_spike_numbers(z)/(analysis_window_post_offset/sampling_rate);
        
        %At the end of each recording, avarage the spike frequency over all
        %bouts recorded in that mat file
        if z == length(spike_vector_per_flyID{i,5}(:,1))
            spike_frequency_per_mat_post_offset(i) = mean(temp_post_offset_spike_frequency);
            avaraged_number_of_spikes_per_mat_post_offset(i) = mean(temp_post_offset_spike_numbers);
            clearvars temp_post_offset_spike_numbers(z) temp_post_offset_spike_frequency
        end
    end
    %At the end of all recorded mat files determine the average spike frequency and
    %numer of spikes for all mat files (post_offset)
    if i == length(spike_vector_per_flyID(:,1))
        all_average_spike_frequency_post_offset =  mean(spike_frequency_per_mat_post_offset);
        all_average_number_of_spikes_post_offset =  mean(avaraged_number_of_spikes_per_mat_post_offset);
    end
    
    
end


%% _PLOT BAR CHART PRE/POST FLIGHT
X = categorical({'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'});
X = reordercats(X,{'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'});
Y = [all_average_spike_frequency_pre_onset all_average_spike_frequency_post_onset all_average_spike_frequency_pre_offset all_average_spike_frequency_post_offset];
pre_post_flight = figure('color', figure_plotting_color);
hold on
ax_sp1 = gca;
if plotting_figures_in_black == 1
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
else
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
end
b = bar(X,Y);
hold on
b.FaceColor = 'flat';
b.CData(1,:) = [0.3010 0.7450 0.9330];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.9290 0.6940 0.1250];
b.CData(4,:) = [0.8500 0.3250 0.0980];

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18);
if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 10.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 10.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 10.eps'
    end
end

%% PLOT Swarmchart PRE/POST FLIGHT ON and OFF
swarm_figure = figure('color', figure_plotting_color);
hold on
if plotting_figures_in_black == 1
    ax_sp1 = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
else
    ax_sp1 = gca;
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
end
%Pre Onset
swarmchart(ones(length(spike_frequency_per_mat_pre_onset(1,:)),1), spike_frequency_per_mat_pre_onset(1,:), [], [0.3010 0.7450 0.9330])
hold on
plot(0.6:0.001:1.4,all_average_spike_frequency_pre_onset, 'Marker', '|', 'MarkerEdgeColor', [0.3010 0.7450 0.9330])
%Post Onset
swarmchart(ones(length(spike_frequency_per_mat_post_onset(1,:)),1)+1, spike_frequency_per_mat_post_onset(1,:), [], [0 0.4470 0.7410])
hold on
plot(1.6:0.001:2.4,all_average_spike_frequency_post_onset, 'Marker', '|', 'MarkerEdgeColor', [0 0.4470 0.7410])
%Pre Offset
swarmchart(ones(length(spike_frequency_per_mat_pre_offset(1,:)),1)+2, spike_frequency_per_mat_pre_offset(1,:), [], [0.9290 0.6940 0.1250])
hold on
plot(2.6:0.001:3.4,all_average_spike_frequency_pre_offset, 'Marker', '|', 'MarkerEdgeColor', [0.9290 0.6940 0.1250])
%Post Offset
swarmchart(ones(length(spike_frequency_per_mat_post_offset(1,:)),1)+3, spike_frequency_per_mat_post_offset(1,:), [], [0.8500 0.3250 0.0980])
hold on
plot(3.6:0.001:4.4,all_average_spike_frequency_post_offset, 'Marker', '|', 'MarkerEdgeColor', [0.8500 0.3250 0.0980])

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',12);
xticks([1 2 3 4])
xticklabels({'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'})
hold off
if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 11_MEAN.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 11_MEAN.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 11_MEAN.eps'
    end
end








%% PLOT Swarmchart PRE/POST FLIGHT ON and OFF MEDIAN not MEAN
swarm_figure = figure('color', figure_plotting_color);
hold on
if plotting_figures_in_black == 1
    ax_sp1 = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
else
    ax_sp1 = gca;
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
end
%Pre Onset
swarmchart(ones(length(spike_frequency_per_mat_pre_onset(1,:)),1), spike_frequency_per_mat_pre_onset(1,:), [], [0.3010 0.7450 0.9330])
hold on
plot(0.6:0.001:1.4,median(spike_frequency_per_mat_pre_onset), 'Marker', '|', 'MarkerEdgeColor', [0.3010 0.7450 0.9330])
%Post Onset
swarmchart(ones(length(spike_frequency_per_mat_post_onset(1,:)),1)+1, spike_frequency_per_mat_post_onset(1,:), [], [0 0.4470 0.7410])
hold on
plot(1.6:0.001:2.4,median(spike_frequency_per_mat_post_onset), 'Marker', '|', 'MarkerEdgeColor', [0 0.4470 0.7410])
%Pre Offset
swarmchart(ones(length(spike_frequency_per_mat_pre_offset(1,:)),1)+2, spike_frequency_per_mat_pre_offset(1,:), [], [0.9290 0.6940 0.1250])
hold on
plot(2.6:0.001:3.4,median(spike_frequency_per_mat_pre_offset), 'Marker', '|', 'MarkerEdgeColor', [0.9290 0.6940 0.1250])
%Post Offset
swarmchart(ones(length(spike_frequency_per_mat_post_offset(1,:)),1)+3, spike_frequency_per_mat_post_offset(1,:), [], [0.8500 0.3250 0.0980])
hold on
plot(3.6:0.001:4.4,median(spike_frequency_per_mat_post_offset), 'Marker', '|', 'MarkerEdgeColor', [0.8500 0.3250 0.0980])

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',12);
xticks([1 2 3 4])
xticklabels({'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'})
hold off
if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 11_MEDIAN.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 11_MEDIAN.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 11_MEDIAN.eps'
    end
end


%% Creating a y_vector for plotting spikes in spike graph (below section)
clearvars y_vector_for_plotting_spike_events_graph
idx_line = 0;
for i=1:length(pre_onset_relative_spike_times(:,1)) %for all Cells
   	idx_line = idx_line -1;
    
    %of that fly, determine the maximum number of y_values needed for all
    %on and offset data
    [~, col1] = find(~cellfun('isempty', pre_onset_relative_spike_times(i, :))==1);
    length_of_current_spike_vector_pre_onset = max(col1); 
    
    [~, col2] = find(~cellfun('isempty', post_onset_relative_spike_times(i, :))==1);
    length_of_current_spike_vector_post_onset = max(col2); 
    
    [~, col3] = find(~cellfun('isempty', pre_offset_relative_spike_times(i, :))==1);
    length_of_current_spike_vector_pre_offset = max(col3); 
    
    [~, col4] = find(~cellfun('isempty', post_offset_relative_spike_times(i, :))==1);
    length_of_current_spike_vector_post_offset = max(col4); 
    
    %determine the maxium number of rows in later spike event graph needed for
    %that fly (max of all 4 values)
%     I2 = unique([col1 col2 col3 col4]);
%     I2 = length(I2)
    [M,I] = max([length_of_current_spike_vector_pre_onset length_of_current_spike_vector_post_onset length_of_current_spike_vector_pre_offset length_of_current_spike_vector_post_offset]);

    %I is the one with the most rows needed (if 4 than post offset...). If
    %it is 4, than the current y vector for that fly needs to be at least 4
    %rows long
    for z1 = 1:M
    y_vector_for_plotting_spike_events_graph{i, z1}(1) = idx_line-0.4;
    y_vector_for_plotting_spike_events_graph{i, z1}(2) = idx_line+0.4;
    idx_line = idx_line - 1;
    
    end
    
end


%%  Plot Spike Event Graph Change this whole section
% selection_for_plotting = %[110, 94];%, 112];$NEW, plots these fly in magenta, those flis are later plotted in more detail
color = colormap(hsv(length(spike_vector_per_flyID(:,1))));
idx_line = +1;% Used later to plot the Y- Values correctly, must be getting minus to look like the plot before (plotten from up to down)

%ONSET
Spike_Event_Figure_onset = figure('name', 'Spike Event Figure (Onset)', 'Color', figure_plotting_color);
set(Spike_Event_Figure_onset, 'position', [1, 1, 800, 1500]);
hold on

%Prior Onset
for i=1:length(pre_onset_relative_spike_times(:,1)) %for all Cells
    idx_line = idx_line-1;
    
    %Determine length of non empty values in spike times (to exclude empty
    %ones that are created because other rows might have more spikes)
    [~, col] = find(~cellfun('isempty', y_vector_for_plotting_spike_events_graph(i, :))==1); %NEW 2021_09_11
    length_of_current_spike_vector = max(col); %sum(~cellfun('isempty', pre_onset_relative_spike_times(i, :)));
    
    for z = 1:length_of_current_spike_vector%for all bouts recorded in that cell
        idx_line = idx_line-1;
        x_vector = [pre_onset_relative_spike_times{i,z}]; %position of spikes
        x_vector(2,:)= x_vector; % create a vector to vizualize
        
        
        y_vector = zeros(length(x_vector(1,:)),2)';
        y_vector(1,:) = y_vector_for_plotting_spike_events_graph{i, z}(1); %NEW 2021_09_11 idx_line-0.4;
        y_vector(2,:) = y_vector_for_plotting_spike_events_graph{i, z}(2);%idx_line+0.4;
        plot(x_vector, y_vector, 'Color', color(i,1:end)) 
        %Save Color Code for this fly in Spikevector Fly ID for later plot
        spike_vector_per_flyID{i,10} = color(i,1:end);
    end
    
end
clearvars x_vector y_vector length_of_current_spike_vector


%Post Onset
idx_line = +1;
for i=1:length(post_onset_relative_spike_times(:,1)) %for all Cells
    idx_line = idx_line-1;
    
    %Determine length of non empty values in spike times (to exclude empty
    %ones that are created because other rows might have more spikes)
    [~, col] = find(~cellfun('isempty', y_vector_for_plotting_spike_events_graph(i, :))==1); %NEW 2021_09_11
    length_of_current_spike_vector = max(col); % sum(~cellfun('isempty', post_onset_relative_spike_times(i, :)));
    
    for z = 1:length_of_current_spike_vector%for all bouts recorded in that cell
        idx_line = idx_line-1;
        x_vector = [post_onset_relative_spike_times{i,z}] + analysis_window_pre_onset/sampling_rate; %position of spikes
        x_vector(2,:)= x_vector; % create a vector to vizualize
        
        
        y_vector = zeros(length(x_vector(1,:)),2)';
%         y_vector(1,:) = idx_line-0.4;
%         y_vector(2,:) = idx_line+0.4;
        y_vector(1,:) = y_vector_for_plotting_spike_events_graph{i, z}(1); %NEW 2021_09_11 idx_line-0.4;
        y_vector(2,:) = y_vector_for_plotting_spike_events_graph{i, z}(2);%idx_line+0.4;
        plot(x_vector, y_vector, 'Color', color(i,1:end))
        %Save Color Code for this fly in Spikevector Fly ID for later plot
        spike_vector_per_flyID{i,10} = color(i,1:end);
        
        
        
        %Plot Onset
        plot(analysis_window_pre_onset/sampling_rate, idx_line, 'Marker', '|', 'MarkerEdgeColor', 'blue')
        if z == length_of_current_spike_vector
            plot((analysis_window_pre_onset/sampling_rate), idx_line-1, 'Marker', '|', 'MarkerEdgeColor', 'blue')
        end
    end
end

if plotting_figures_in_black == 1
    ax_spikevector_plot = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_spikevector_plot,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_spikevector_plot,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    ylabel('Walking Bouts  (n)');
    xlabel('Time (s)');
else
    ax_spikevector_plot = gca;
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_spikevector_plot,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_spikevector_plot,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    ylabel('Walking Bouts');
    xlabel('Time (s)');
end
a = get(gca,'XTickLabel');
set(gca,'FontName','Helvetica','fontsize',12);

if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 8.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 8.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 8.eps'
    end
end


clearvars x_vector y_vector length_of_current_spike_vector






%OFFSET
idx_line = +1;% Used later to plot the Y- Values correctly, must be getting minus to look like the plot before (plotten from up to down)

Spike_Event_Figure_offset = figure('name', 'Spike Event Figure (Offset)', 'Color', figure_plotting_color);
set(Spike_Event_Figure_offset, 'position', [1, 1, 800, 1500]);
hold on

%Prior Offset
for i=1:length(pre_offset_relative_spike_times(:,1)) %for all Cells
    idx_line = idx_line-1;
    
    %Determine length of non empty values in spike times (to exclude empty
    %ones that are created because other rows might have more spikes)
    [~, col] = find(~cellfun('isempty', y_vector_for_plotting_spike_events_graph(i, :))==1); %NEW 2021_09_11
    length_of_current_spike_vector = max(col); %sum(~cellfun('isempty', pre_onset_relative_spike_times(i, :)));
    
    for z = 1:length_of_current_spike_vector%for all bouts recorded in that cell
        idx_line = idx_line-1;
        x_vector = [pre_offset_relative_spike_times{i,z}]; %position of spikes
        x_vector(2,:)= x_vector; % create a vector to vizualize
        
        
        y_vector = zeros(length(x_vector(1,:)),2)';
        y_vector(1,:) = y_vector_for_plotting_spike_events_graph{i, z}(1); %NEW 2021_09_11 idx_line-0.4;
        y_vector(2,:) = y_vector_for_plotting_spike_events_graph{i, z}(2);%idx_line+0.4;
%         y_vector(1,:) = idx_line-0.4;
%         y_vector(2,:) = idx_line+0.4;
        plot(x_vector, y_vector, 'Color', color(i,1:end))
        %Save Color Code for this fly in Spikevector Fly ID for later plot
        spike_vector_per_flyID{i,10} = color(i,1:end);
    end
    
end
clearvars x_vector y_vector length_of_current_spike_vector


%Post Offset
idx_line = +1;
for i=1:length(post_offset_relative_spike_times(:,1)) %for all Cells
    idx_line = idx_line-1;
    
    %Determine length of non empty values in spike times (to exclude empty
    %ones that are created because other rows might have more spikes)
    [~, col] = find(~cellfun('isempty', y_vector_for_plotting_spike_events_graph(i, :))==1); %NEW 2021_09_11
    length_of_current_spike_vector = max(col); % sum(~cellfun('isempty', post_onset_relative_spike_times(i, :)));
    
    for z = 1:length_of_current_spike_vector%for all bouts recorded in that cell
        idx_line = idx_line-1;
        x_vector = [post_offset_relative_spike_times{i,z}] + analysis_window_pre_offset/sampling_rate; %position of spikes
        x_vector(2,:)= x_vector; % create a vector to vizualize
        
        
        y_vector = zeros(length(x_vector(1,:)),2)';
        y_vector(1,:) = y_vector_for_plotting_spike_events_graph{i, z}(1); %NEW 2021_09_11 idx_line-0.4;
        y_vector(2,:) = y_vector_for_plotting_spike_events_graph{i, z}(2);%idx_line+0.4;
%         y_vector(1,:) = idx_line-0.4;
%         y_vector(2,:) = idx_line+0.4;
        plot(x_vector, y_vector, 'Color', color(i,1:end))
        %Save Color Code for this fly in Spikevector Fly ID for later plot
        spike_vector_per_flyID{i,10} = color(i,1:end);
        
        
        
        %Plot Onset
        plot(analysis_window_pre_offset/sampling_rate, idx_line, 'Marker', '|', 'MarkerEdgeColor', 'blue')
        if z == length_of_current_spike_vector
            plot((analysis_window_pre_offset/sampling_rate), idx_line-1, 'Marker', '|', 'MarkerEdgeColor', 'blue')
        end
    end
end


if plotting_figures_in_black == 1
    ax_spikevector_plot = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_spikevector_plot,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_spikevector_plot,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    ylabel('Walking Bouts');
    xlabel('Time (s)');
else
    ax_spikevector_plot = gca;
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_spikevector_plot,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_spikevector_plot,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    ylabel('Walking Bouts');
    xlabel('Time (s)');
end
a = get(gca,'XTickLabel');
set(gca,'FontName','Helvetica','fontsize',12);

if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 9.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 9.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 9.eps'
    end
end

clearvars x_vector y_vector length_of_current_spike_vector





%% Plotting lines of indiv. Mat files into that Graph
color = colormap(hsv(length(spike_frequency_per_mat_pre_onset(1,:))));
spikerate_fig=figure('name','Spikerate Figure', 'color', figure_plotting_color);
set(spikerate_fig, 'position', [1, 1, 800, 1500]);
hold on
for i = 1:length(spike_frequency_per_mat_pre_onset(1,:))
    %creating vectors for plotting)
    tmp_plot1_x = [1 2];
    tmp_plot1_y =  [spike_frequency_per_mat_pre_onset(1,i) spike_frequency_per_mat_post_onset(1,i)];
    plot(tmp_plot1_x, tmp_plot1_y, 'Color', color(i,1:end))
    hold on
    tmp_plot2_x = [3 4];
    tmp_plot2_y =  [spike_frequency_per_mat_pre_offset(1,i) spike_frequency_per_mat_post_offset(1,i)];
    plot(tmp_plot2_x, tmp_plot2_y, 'Color', color(i,1:end))
    hold on
    % Plot the in between
    tmp_plot3_x = [2 3];
    tmp_plot3_y =  [spike_frequency_per_mat_post_onset(1,i) spike_frequency_per_mat_pre_offset(1,i)];
    plot(tmp_plot3_x, tmp_plot3_y, ':' ,'Color', [0.25 0.25 0.25])
    %     plot(tmp_plot3_x, tmp_plot3_y, ':' ,'Color', color(i,1:end))
    xlim([0.5 4.5]);
    ylim([-0.1 inf]);
end

if plotting_figures_in_black == 1
    ax_sp1 = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
    xlabel('State');
else
    ax_sp1 = gca;
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
    xlabel('State');
end

if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 12.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 12.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 12.eps'
    end
end

%% New Summarizing Graph
%Color Code for Flies
color = colormap(gray(length(spike_frequency_per_mat_pre_onset(1,:)) + length(spike_frequency_per_mat_pre_onset(1,:))/2));%devided by 6 to get rid of the white values
color = color(length(spike_frequency_per_mat_pre_onset(1,:))/2:end,1:3); % for black background
% color = color(length(number_of_averaged_spikes_per_mat(:,1))/6:end,1:3);
% %for white and put 'ength(number_of_averaged_spikes_per_mat(:,1))/2))' to
% 6
color = color(randperm(size(color,1)),:); %shuffle them
% color = colormap(gray(length(number_of_averaged_spikes_per_mat(:,1))));
%Figure
spikerate_fig_all=figure('name','Spikerate Figure', 'Color', figure_plotting_color);
%Setting Axes
hold on
ax = gca;
ax.YLabel.String = 'Mean Frequency [Hz]';
ax.XLabel.String = 'Walking State';
if plotting_figures_in_black == 1
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
else
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out')
end
set(spikerate_fig_all, 'position', [1, 1, 800, 1500]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18);
%Y Values for Bar
Y = [all_average_spike_frequency_pre_onset all_average_spike_frequency_post_onset all_average_spike_frequency_pre_offset all_average_spike_frequency_post_offset];
%Setting Color of Bars
b = bar([1 2 3 4],Y, 0.4, 'EdgeColor', 'k', 'FaceAlpha', 0.8, 'LineWidth', 2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.3010 0.7450 0.9330];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.9290 0.6940 0.1250];
b.CData(4,:) = [0.8500 0.3250 0.0980];
title('Summarizing Result');
hold on


%Actual plotting
for i = 1:length(spike_frequency_per_mat_pre_onset(1,:))
    %creating vectors for plotting)
    tmp_plot1_x = [1 2];
    tmp_plot1_y =  [spike_frequency_per_mat_pre_onset(1,i) spike_frequency_per_mat_post_onset(1,i)];
    plot(tmp_plot1_x, tmp_plot1_y, 'Color', color(i,1:end), 'LineWidth',1.5)
    hold on
    tmp_plot2_x = [3 4];
    tmp_plot2_y =  [spike_frequency_per_mat_pre_offset(1,i) spike_frequency_per_mat_post_offset(1,i)];
    plot(tmp_plot2_x, tmp_plot2_y, 'Color', color(i,1:end), 'LineWidth',1.5)
    hold on
    % Plot the in between
    tmp_plot3_x = [2 3];
    tmp_plot3_y =  [spike_frequency_per_mat_post_onset(1,i) spike_frequency_per_mat_pre_offset(1,i)];
    %     plot(tmp_plot3_x, tmp_plot3_y, '--' ,'Color', [0.25 0.25 0.25], 'LineWidth',1)
    plot(tmp_plot3_x, tmp_plot3_y, ':' ,'Color', '#D3D3D3', 'LineWidth',1)
    xlim([0.5 4.5]);
    ylim([-0.1 inf]);
end
%Plot Mean
mean_spike = [all_average_spike_frequency_pre_onset all_average_spike_frequency_post_onset all_average_spike_frequency_pre_offset all_average_spike_frequency_post_offset];%mean(number_of_averaged_spikes_per_mat/(input_time_width_average_window*binsize_factor),1);
plot([1 2 3 4], mean_spike, 'Color', 'r', 'LineWidth',3);
%Changing Tick Labels of X
xticks([1 2 3 4]);
set(gca, 'XTickLabel',{'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'})
%Plot N and n values
txt = ['Number of Mat Files: ' num2str(length(spike_vector_per_flyID(:,1))) ', number of flight bouts: ' num2str( length(Analysis_Data_Flight(:,1)))];
t = text(1,-0.3,txt);
set(t, 'horizontalAlignment', 'left')
txt = ['Number of Flies: ' num2str(number_of_flies) ', number of cells: ' num2str(number_of_cells)];
t = text(1,-0.4,txt);
set(t, 'horizontalAlignment', 'left')
if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 13.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 13.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 13.eps'
    end
end









%% Summarizing Graph Median and IQR
spikerate_fig_all = figure('name','Spikerate Figure', 'Color', figure_plotting_color);
hold on
set(spikerate_fig_all, 'position', [1, 1, 800, 1500]);
%Setting Axes
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18);
title('Summarizing Result');

%------- Boxplot
hBox = boxplot2([spike_frequency_per_mat_pre_onset', spike_frequency_per_mat_post_onset', spike_frequency_per_mat_pre_offset', spike_frequency_per_mat_post_offset'], [1, 2, 3, 4]);
hold on
%Calculating upper an lower prctl
lowerprctl(1) = prctile(spike_frequency_per_mat_pre_onset, 5, 'all');
upperprtcl(1) = prctile(spike_frequency_per_mat_pre_onset, 95, 'all');
lowerprctl(2) = prctile(spike_frequency_per_mat_post_onset, 5, 'all');
upperprtcl(2) = prctile(spike_frequency_per_mat_post_onset, 95, 'all');
lowerprctl(3) = prctile(spike_frequency_per_mat_pre_offset, 5, 'all');
upperprtcl(3) = prctile(spike_frequency_per_mat_pre_offset, 95, 'all');
lowerprctl(4) = prctile(spike_frequency_per_mat_post_offset, 5, 'all');
upperprtcl(4) = prctile(spike_frequency_per_mat_post_offset, 95, 'all');

%plotting
for j = 1:4
    y = get(hBox.uwhis(j),'YData');
    set(hBox.uwhis(j),'YData',[y(1),upperprtcl(j)]);
    y = get(hBox.lwhis(j),'YData');
    set(hBox.lwhis(j),'YData',[lowerprctl(j),y(2)]);
    
    set(hBox.uadj(j),'YData',[upperprtcl(j),upperprtcl(j)]);
    
    set(hBox.ladj(j),'YData',[lowerprctl(j),lowerprctl(j)]);
end

if plotting_figures_in_black == 1
    ax_sp1 = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency (Hz)');
    xlabel('State');
else
    ax_sp1 = gca;
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    ylabel('Spike Frequency  (Hz)');
    xlabel('State');
end

if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 13_boxplot.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 13_boxplot.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 13_boxplot.eps'
    end
end

fid = fopen('Figure 13_boxplot.txt','w');
fprintf(fid, '%6s %12s\r\n', 'Median (Pre On/Post On/Pre Off/Post Off)', [num2str(median(spike_frequency_per_mat_pre_onset)) ' ' num2str(median(spike_frequency_per_mat_post_onset)) ' ' num2str(median(spike_frequency_per_mat_pre_offset)) ' ' num2str(median(spike_frequency_per_mat_post_offset))]);
fprintf(fid, '%6s %12s\r\n', 'Mean (Pre On/Post On/Pre Off/Post Off)', [num2str(mean(spike_frequency_per_mat_pre_onset)) ' ' num2str(mean(spike_frequency_per_mat_post_onset)) ' ' num2str(mean(spike_frequency_per_mat_pre_offset)) ' ' num2str(mean(spike_frequency_per_mat_post_offset))]);
fprintf(fid, '%6s %12s\r\n', 'Max pre on:', num2str(max(spike_frequency_per_mat_pre_onset)));
fprintf(fid, '%6s %12s\r\n', 'Max post on:', num2str(max(spike_frequency_per_mat_post_onset)));
fprintf(fid, '%6s %12s\r\n', 'Max pre off:', num2str(max(spike_frequency_per_mat_pre_offset)));
fprintf(fid, '%6s %12s\r\n', 'Max post off:', num2str(max(spike_frequency_per_mat_post_offset)));
fprintf(fid, '%6s %12s\r\n', 'Number of Mat Files', num2str(length(spike_vector_per_flyID(:,1))));
fprintf(fid, '%6s %12s\r\n', 'Number of flight bouts:', num2str( length(Analysis_Data_Flight(:,1))));
fprintf(fid, '%6s %12s\r\n', 'Number of cells:', num2str(number_of_cells));

temp = Analysis_Data_Flight(:,1);
clear temp2
for z = 1:length(temp)
temp2{z} = temp{z}(1:6);
end
[d, id] = unique(temp2, 'stable');
length(id);

fprintf(fid, '%6s %12s\r\n', 'Number of Flies:', num2str(length(id)));
fclose(fid);


% Figure with single cells

color = colormap(hsv(length(spike_frequency_per_mat_pre_onset)));
% color = [0.75 0.75 0.75];
figure
%Actual plotting
for i = 1:length(spike_frequency_per_mat_pre_onset(1,:))
    %creating vectors for plotting)
    tmp_plot1_x = [1 2];
    tmp_plot1_y =  [spike_frequency_per_mat_pre_onset(1,i) spike_frequency_per_mat_post_onset(1,i)];
    plot(tmp_plot1_x, tmp_plot1_y, 'Color', color(i,1:end), 'LineWidth',1.5)
    hold on
    tmp_plot2_x = [3 4];
    tmp_plot2_y =  [spike_frequency_per_mat_pre_offset(1,i) spike_frequency_per_mat_post_offset(1,i)];
    plot(tmp_plot2_x, tmp_plot2_y, 'Color', color(i,1:end), 'LineWidth',1.5)
    hold on
    % Plot the in between
    tmp_plot3_x = [2 3];
    tmp_plot3_y =  [spike_frequency_per_mat_post_onset(1,i) spike_frequency_per_mat_pre_offset(1,i)];
    %     plot(tmp_plot3_x, tmp_plot3_y, '--' ,'Color', [0.25 0.25 0.25], 'LineWidth',1)
    plot(tmp_plot3_x, tmp_plot3_y, ':' ,'Color', '#D3D3D3', 'LineWidth',1)
    xlim([0.5 4.5]);
    ylim([0 25]);
end
%Plot Mean
mean_spike = [all_average_spike_frequency_pre_onset all_average_spike_frequency_post_onset all_average_spike_frequency_pre_offset all_average_spike_frequency_post_offset];%mean(number_of_averaged_spikes_per_mat/(input_time_width_average_window*binsize_factor),1);
plot([1 2 3 4], mean_spike, 'Color', 'r', 'LineWidth',3);
%Changing Tick Labels of X
xticks([1 2 3 4]);
set(gca, 'XTickLabel',{'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'})
%Plot N and n values
txt = ['Number of Mat Files: ' num2str(length(spike_vector_per_flyID(:,1))) ', number of flight bouts: ' num2str( length(Analysis_Data_Flight(:,1)))];
t = text(1,-0.3,txt);
set(t, 'horizontalAlignment', 'left')
txt = ['Number of Flies: ' num2str(number_of_flies) ', number of cells: ' num2str(number_of_cells)];
t = text(1,-0.4,txt);
set(t, 'horizontalAlignment', 'left')

if plotting_figures_in_black == 1
    ax_sp1 = gca;
    set(gca,'color',[0 0 0]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','w', 'XColor','w')
    set(gca,'TickDir','out');
    ylabel('Mean Spike Frequency per Cell (Hz)');
    xlabel('State');
else
    ax_sp1 = gca;
    set(gca,'color',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    set(gca, 'Xcolor', 'k');
    set(gcf,'inverthardcopy','off');
    set(ax_sp1,'YColor','k', 'XColor','k')
    set(gca,'TickDir','out');
    ylabel('Mean Spike Frequency per Cell (Hz)');
    xlabel('State');
end


if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 13.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 13.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 13.eps'
    end
end






%% Normalization to pre-onsett (in percent)
clearvars tmp_plot1_x tmp_plot1_y tmp_plot2_y tmp_plot2_x
number_of_averaged_spikes_per_mat_normalized_percent = [spike_frequency_per_mat_pre_onset; spike_frequency_per_mat_post_onset; spike_frequency_per_mat_pre_offset; spike_frequency_per_mat_post_offset];
number_of_averaged_spikes_per_mat_normalized_dif = [spike_frequency_per_mat_pre_onset; spike_frequency_per_mat_post_onset; spike_frequency_per_mat_pre_offset; spike_frequency_per_mat_post_offset];

for i = 1:length(number_of_averaged_spikes_per_mat_normalized_percent(1,:))
    number_of_averaged_spikes_per_mat_normalized_dif(2, i) = number_of_averaged_spikes_per_mat_normalized_dif(2, i)-number_of_averaged_spikes_per_mat_normalized_dif(1, i);
    number_of_averaged_spikes_per_mat_normalized_dif(3, i) = number_of_averaged_spikes_per_mat_normalized_dif(3, i)-number_of_averaged_spikes_per_mat_normalized_dif(1, i);
    number_of_averaged_spikes_per_mat_normalized_dif(4, i) = number_of_averaged_spikes_per_mat_normalized_dif(4, i)-number_of_averaged_spikes_per_mat_normalized_dif(1, i);
    number_of_averaged_spikes_per_mat_normalized_dif(1, i) = number_of_averaged_spikes_per_mat_normalized_dif(1, i)-number_of_averaged_spikes_per_mat_normalized_dif(1, i);
    
    number_of_averaged_spikes_per_mat_normalized_percent(2, i) = number_of_averaged_spikes_per_mat_normalized_percent(2, i)/number_of_averaged_spikes_per_mat_normalized_percent(1, i);
    number_of_averaged_spikes_per_mat_normalized_percent(3, i) = number_of_averaged_spikes_per_mat_normalized_percent(3, i)/number_of_averaged_spikes_per_mat_normalized_percent(1, i);
    number_of_averaged_spikes_per_mat_normalized_percent(4, i) = number_of_averaged_spikes_per_mat_normalized_percent(4, i)/number_of_averaged_spikes_per_mat_normalized_percent(1, i);
    number_of_averaged_spikes_per_mat_normalized_percent(1, i) = number_of_averaged_spikes_per_mat_normalized_percent(1, i)/number_of_averaged_spikes_per_mat_normalized_percent(1, i);
end
%Remove inf and NaN
number_of_averaged_spikes_per_mat_normalized_percent(isinf(number_of_averaged_spikes_per_mat_normalized_percent)|isnan(number_of_averaged_spikes_per_mat_normalized_percent)) = 0;


%*******************PLOTTING NORM PERCENT
spikerate_fig_all_norm_percent=figure('name','Spikerate Figure - Normalization %');
set(spikerate_fig_all_norm_percent, 'position', [1, 1, 800, 1500]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18);
%Y Values for Bar
Y = [all_average_spike_frequency_pre_onset/all_average_spike_frequency_pre_onset all_average_spike_frequency_post_onset/all_average_spike_frequency_pre_onset all_average_spike_frequency_pre_offset/all_average_spike_frequency_pre_onset all_average_spike_frequency_post_offset/all_average_spike_frequency_pre_onset];
%Setting Color of Bars
b = bar([1 2 3 4],Y, 0.4, 'EdgeColor', 'k', 'FaceAlpha', 0.5, 'LineWidth', 2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.3010 0.7450 0.9330];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.9290 0.6940 0.1250];
b.CData(4,:) = [0.8500 0.3250 0.0980];
title('Summarizing Result - Norm %');
hold on
%Setting Axes
ax = gca;
ax.YLabel.String = 'Mean Frequency [Hz]';
ax.XLabel.String = 'Flight State';
%Actual plotting
for i = 1:length(number_of_averaged_spikes_per_mat_normalized_percent(1,:))
    %creating vectors for plotting)
    tmp_plot1_x = [1 2];
    tmp_plot1_y =  [number_of_averaged_spikes_per_mat_normalized_percent(1, i) number_of_averaged_spikes_per_mat_normalized_percent(2, i)];
    plot(tmp_plot1_x, tmp_plot1_y, 'Color', color(i,1:end), 'LineWidth',1)
    hold on
    tmp_plot2_x = [3 4];
    tmp_plot2_y =  [number_of_averaged_spikes_per_mat_normalized_percent(3, i) number_of_averaged_spikes_per_mat_normalized_percent(4, i)];
    plot(tmp_plot2_x, tmp_plot2_y, 'Color', color(i,1:end), 'LineWidth',1)
    hold on
    % Plot the in between
    tmp_plot3_x = [2 3];
    tmp_plot3_y =  [number_of_averaged_spikes_per_mat_normalized_percent(2, i) number_of_averaged_spikes_per_mat_normalized_percent(3, i)];
    %     plot(tmp_plot3_x, tmp_plot3_y, '--' ,'Color', [0.25 0.25 0.25], 'LineWidth',1)
    plot(tmp_plot3_x, tmp_plot3_y, ':' ,'Color', '#D3D3D3', 'LineWidth',0.5)
    xlim([0.5 4.5]);
    %     ylim([-0.1 inf]);
end
%Plot Mean
mean_spike = mean(number_of_averaged_spikes_per_mat_normalized_percent,2);

plot([1 2 3 4], mean_spike, 'Color', 'black', 'LineWidth',3);
%Changing Tick Labels of X
set(gca, 'XTickLabel',{'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'})
%Plot N and n values
txt = ['Number of Mat Files: ' num2str(length(spike_vector_per_flyID(:,1))) ', number of flight bouts: ' num2str( length(Analysis_Data_Flight(:,1)))];
t = text(1,-0.3,txt);
set(t, 'horizontalAlignment', 'left')
txt = ['Number of Flies: ' num2str(number_of_flies) ', number of cells: ' num2str(number_of_cells)];
t = text(1,-0.4,txt);
set(t, 'horizontalAlignment', 'left')
if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 14.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 14.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 14.eps'
    end
end

%******************PLOTTING NOMR DIFF

spikerate_fig_all_norm_dif=figure('name','Spikerate Figure - Normalization diff');
set(spikerate_fig_all_norm_dif, 'position', [1, 1, 800, 1500]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18);
%Y Values for Bar
Y = [all_average_spike_frequency_pre_onset-all_average_spike_frequency_pre_onset all_average_spike_frequency_post_onset-all_average_spike_frequency_pre_onset all_average_spike_frequency_pre_offset-all_average_spike_frequency_pre_onset all_average_spike_frequency_post_offset-all_average_spike_frequency_pre_onset];
%Setting Color of Bars
b = bar([1 2 3 4],Y, 0.4, 'EdgeColor', 'k', 'FaceAlpha', 0.5, 'LineWidth', 2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.3010 0.7450 0.9330];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.9290 0.6940 0.1250];
b.CData(4,:) = [0.8500 0.3250 0.0980];
title('Summarizing Result - Norm Diff');
hold on
%Setting Axes
ax = gca;
ax.YLabel.String = 'Mean Frequency [Hz]';
ax.XLabel.String = 'Flight State';
%Actual plotting
for i = 1:length(number_of_averaged_spikes_per_mat_normalized_percent(1,:))
    %creating vectors for plotting)
    tmp_plot1_x = [1 2];
    tmp_plot1_y =  [number_of_averaged_spikes_per_mat_normalized_dif(1, i) number_of_averaged_spikes_per_mat_normalized_dif(2, i)];
    plot(tmp_plot1_x, tmp_plot1_y, 'Color', color(i,1:end), 'LineWidth',1)
    hold on
    tmp_plot2_x = [3 4];
    tmp_plot2_y =  [number_of_averaged_spikes_per_mat_normalized_dif(3, i) number_of_averaged_spikes_per_mat_normalized_dif(4,  i)];
    plot(tmp_plot2_x, tmp_plot2_y, 'Color', color(i,1:end), 'LineWidth',1)
    hold on
    % Plot the in between
    tmp_plot3_x = [2 3];
    tmp_plot3_y =  [number_of_averaged_spikes_per_mat_normalized_dif(2, i) number_of_averaged_spikes_per_mat_normalized_dif(3, i)];
    %     plot(tmp_plot3_x, tmp_plot3_y, '--' ,'Color', [0.25 0.25 0.25], 'LineWidth',1)
    plot(tmp_plot3_x, tmp_plot3_y, ':' ,'Color', '#D3D3D3', 'LineWidth',0.5)
    xlim([0.5 4.5]);
    %     ylim([-0.1 inf]);
end
%Plot Mean
mean_spike = mean(number_of_averaged_spikes_per_mat_normalized_dif, 2);
plot([1 2 3 4], mean_spike, 'Color', 'black', 'LineWidth',3);
%Changing Tick Labels of X
set(gca, 'XTickLabel',{'Pre-Onset','Post-Onset','Pre-Offset','Post-Offset'})
%Plot N and n values
txt = ['Number of Mat Files: ' num2str(length(spike_vector_per_flyID(:,1))) ', number of flight bouts: ' num2str( length(Analysis_Data_Flight(:,1)))];
t = text(1,-0.3,txt);
set(t, 'horizontalAlignment', 'left')
txt = ['Number of Flies: ' num2str(number_of_flies) ', number of cells: ' num2str(number_of_cells)];
t = text(1,-0.4,txt);
set(t, 'horizontalAlignment', 'left')
if input_save_figures == 1
    pause(0.5)
    filename = [file_modifier '_Figure 15.tif'];
    saveas(gcf,filename)
    if strcmp(file_modifier, 'unclean')
        print -depsc2 -tiff -r300 -painters 'unclean_Figure 15.eps'
    else
        print -depsc2 -tiff -r300 -painters 'clean_Figure 15.eps'
    end
    %     print(gcf,'Figure 13.emf', '-dmeta','-painters')
    %     print(gcf,'Figure 13.epsc','-dmeta','-r300')
end


%% Creating Figure with selected trials for plotting spike plot for puplication
selection_for_plotting = [1:15];%, 112];
%selection_for_plotting = [109, 94];%, 110 112];
color = colormap(bone(length(selection_for_plotting)*2));

% figure
% hold on
offset_for_plotting = 15;

for i = 1:length(selection_for_plotting)
    figure
hold on

sp1 = subplot(1,2,1);
hold on
%Prior Onset
 x = (1:length(Analysis_Data_Flight{selection_for_plotting(i), 43}(Analysis_Data_Flight{selection_for_plotting(i), 44}-analysis_window_pre_onset:Analysis_Data_Flight{selection_for_plotting(i), 44})))/sampling_rate;
        y = Analysis_Data_Flight{selection_for_plotting(i), 43}(Analysis_Data_Flight{selection_for_plotting(i), 44}-analysis_window_pre_onset:Analysis_Data_Flight{selection_for_plotting(i), 44});
%         y = y+(offset_for_plotting*i);
     plot(x, y, 'color', color(i, :))
     
        %Post Onset
        x = (analysis_window_pre_onset+1:analysis_window_pre_onset+analysis_window_post_onset+1)/sampling_rate;
        y = Analysis_Data_Flight{selection_for_plotting(i), 43}(Analysis_Data_Flight{selection_for_plotting(i), 44}:Analysis_Data_Flight{selection_for_plotting(i), 44}+analysis_window_post_onset);
%        y = y+(offset_for_plotting*i);
     plot(x, y, 'color', color(i, :))
     xli = xlim;
     text(xli(2)+0.3, y(end),num2str(y(1)), 'HorizontalAlignment','right','VerticalAlignment','top', 'Color', 'k')%plot membrane pot
        
    axis tight
    ylim([-60 20])
    
    

    sp2 = subplot(1,2,2);
    hold on
    %prior offset
     x = (1:length(Analysis_Data_Flight{selection_for_plotting(i), 43}(Analysis_Data_Flight{selection_for_plotting(i), 45}-analysis_window_pre_offset:Analysis_Data_Flight{selection_for_plotting(i), 45})))/sampling_rate;
     y = Analysis_Data_Flight{selection_for_plotting(i), 43}(Analysis_Data_Flight{selection_for_plotting(i), 45}-analysis_window_pre_offset:Analysis_Data_Flight{selection_for_plotting(i), 45});
%      y = y+(offset_for_plotting*i);
     plot(x, y, 'color', color(i, :))
 
      %Post Offset
        x = (analysis_window_pre_offset+1:analysis_window_pre_offset+analysis_window_post_offset+1)/sampling_rate;
        y = Analysis_Data_Flight{selection_for_plotting(i), 43}(Analysis_Data_Flight{selection_for_plotting(i), 45}:Analysis_Data_Flight{selection_for_plotting(i), 45}+analysis_window_post_offset);
%         y = y+(offset_for_plotting*i);
     plot(x, y, 'color', color(i, :))
         xli = xlim;
     text(xli(2)+0.3, y(end),num2str(y(1)), 'HorizontalAlignment','right','VerticalAlignment','top', 'Color', 'k')%plot membrane pot
        
    axis tight
     ylim([-60 20])
     

     filename = [num2str(Analysis_Data_Flight{selection_for_plotting(i), 1}) '_' num2str(Analysis_Data_Flight{selection_for_plotting(i), 9}/sampling_rate) '_Figure 17.tif'];
    saveas(gcf,filename)
    
    filename = [num2str(Analysis_Data_Flight{selection_for_plotting(i), 1}) '_' num2str(Analysis_Data_Flight{selection_for_plotting(i), 9}/sampling_rate)  '_Figure_17.eps'];
            print(filename, '-depsc2', '-tiff', '-r300', '-painters')
%     if strcmp(file_modifier, 'unclean')
%         print -depsc2 -tiff -r300 -painters 'Figure 17_1.eps'
%     else
%         print -depsc2 -tiff -r300 -painters 'clean_Figure 15.eps'
%     end
end



%% Stats
statistics_input = [spike_frequency_per_mat_pre_onset; spike_frequency_per_mat_post_onset; spike_frequency_per_mat_pre_offset; spike_frequency_per_mat_post_offset];

[A, h] = signrank(statistics_input(1, :), statistics_input(2, :));
% disp(['Prior-Onset vs. Post-Onset: ' num2str(A)]);

[B, h] = signrank(statistics_input(2, :), statistics_input(3, :));
% disp(['Post-Onset vs. Prior-Offset: ' num2str(B)]);

[C, h] = signrank(statistics_input(3, :), statistics_input(4, :));
% disp(['Prior-Offset vs. Post-Offset: ' num2str(C)]);

[D, h] = signrank(statistics_input(1, :), statistics_input(4, :));
% disp(['Prior-Onset vs. Post-Offset: ' num2str(D)]);

% Create text file with stats
T = table(A, B, C, D, 'VariableNames', { 'Prior-Onset vs. Post-Onset', 'Post-Onset vs. Prior-Offset', 'Prior-Offset vs. Post-Offset', 'Prior-Onset vs. Post-Offset'} )
% Write data to text file
writetable(T, 'Stats.txt')


%% ************* Save results for current dataset to mat file *************
%'-v7.3' flag needed for files larger than two GB
disp('Saving results....');

% save(['Analysis_TACHO' '.mat'],'flyID', 'current_File', 'Analysis_Data_Flight', 'sampling_rate', 'current', 'spike_vector_per_flyID', 'tacho', 'spikes', 'intra_smooth', 'average_spikes_flyID', 'average_averag_spike_all_pre_onset', 'average_averag_spike_all_pre_offset', 'average_averag_spike_all_post_onset', 'average_averag_spike_all_post_offset', 'averag_spike_all_pre_onset', 'averag_spike_all_pre_offset', 'averag_spike_all_post_onset', 'averag_spike_all_post_offset', 'average_number_of_spikes_pre_onset', 'average_number_of_spikes_post_onset', 'average_number_of_spikes_pre_offset', 'average_number_of_spikes_post_offset','number_of_averaged_spikes_per_mat','-v7.3');
save(['Analysis_fly_activity_' file_modifier '.mat'],'flyID', 'current_File', 'Analysis_Data_Flight', 'sampling_rate', 'current', 'spike_vector_per_flyID', 'average_spikes_flyID', 'spike_frequency_per_mat_post_offset', 'spike_frequency_per_mat_post_onset', 'spike_frequency_per_mat_pre_offset', 'spike_frequency_per_mat_pre_onset', 'avaraged_number_of_spikes_per_mat_pre_onset', 'avaraged_number_of_spikes_per_mat_post_offset', 'avaraged_number_of_spikes_per_mat_post_onset', 'avaraged_number_of_spikes_per_mat_pre_offset', 'all_average_number_of_spikes_post_offset', 'all_average_number_of_spikes_post_onset', 'all_average_number_of_spikes_pre_offset', 'all_average_number_of_spikes_pre_onset', 'all_average_spike_frequency_post_offset', 'all_average_spike_frequency_post_onset', 'all_average_spike_frequency_pre_offset', 'all_average_spike_frequency_pre_onset', '-v7.3'); %'average_averag_spike_all_pre_onset', 'average_averag_spike_all_pre_offset', 'average_averag_spike_all_post_onset', 'average_averag_spike_all_post_offset', 'averag_spike_all_pre_onset', 'averag_spike_all_pre_offset', 'averag_spike_all_post_onset', 'averag_spike_all_post_offset',

disp('FINISHED saving MAT file');


%% Saving Analysis Paramter
disp('Saving Analysis Paramter')
scriptversion = mfilename('fullpath');

% f = fopen('myFile.txt', 'wt');      % 'wt' - write file, text mode
% formspec = '%f  %f\n';            % two values in a row (\n - line break)
% % formspec = '%f  %f  %f\n';      % tree values in a row
%
%     fprintf(f, formspec, analysis_window_post_offset);
%     sprintf('%0.1f', formspec, num2str(analysis_window_post_offset));
% fclose(f);

fid = fopen('Analysis_Paramter_Main_Analysis.txt','w');
fprintf(fid, '%6s %12s\r\n', 'analysed dataset', analysed_dataset);
fprintf(fid, '%6s %12s\r\n', 'analysis_window_pre_onset', num2str(analysis_window_pre_onset/sampling_rate));
fprintf(fid, '%6s %12s\r\n', 'analysis_window_post_onset', num2str(analysis_window_post_onset/sampling_rate));
fprintf(fid, '%6s %12s\r\n', 'analysis_window_pre_offset', num2str(analysis_window_pre_offset/sampling_rate));
fprintf(fid, '%6s %12s\r\n', 'analysis_window_post_offset', num2str(analysis_window_post_offset/sampling_rate));

fprintf(fid, '%6s %12s\r\n', 'analysis_membrane_potential_prior_onset_in_sec', num2str(analysis_membrane_potential_prior_onset_in_sec));
fprintf(fid, '%6s %12s\r\n', 'analysis_membrane_potential_post_onset', num2str(analysis_membrane_potential_post_onset));
fprintf(fid, '%6s %12s\r\n', 'analysis_membrane_potential_prior_offset', num2str(analysis_membrane_potential_prior_offset));
fprintf(fid, '%6s %12s\r\n', 'analysis_membrane_potential_post_offset_in_sec', num2str(analysis_membrane_potential_post_offset_in_sec));


%fprintf(fid, '%6s %12s\r\n', 'Data in Spike plot are: ', sorted_by);

fprintf(fid, '%6s %12s\r\n', 'min_interflight_interval', num2str(min_interflight_interval));

fprintf(fid, '%6s %12s\r\n', 'min_flight_duration_in_sec (must at least be active this long to be evaluated)', num2str(min_flight_duration_in_sec));
fprintf(fid, '%6s %12s\r\n', 'minimum_spike_amplitude', num2str(minimum_spike_amplitude));
fprintf(fid, '%6s %12s\r\n', 'number_of_cells', num2str(number_of_cells));
fprintf(fid, '%6s %12s\r\n', 'number_of_flies', num2str(number_of_flies));
fprintf(fid, '%6s %12s\r\n', 'smoothing_factor_membrane', num2str(smoothing_factor_membrane));

fprintf(fid, '%6s %12s\r\n', 'Spikes are sorted by post offset');
fprintf(fid, '%6s %12s\r\n', 'scriptversion', scriptversion);


fclose(fid);

disp('FINISHED');






































