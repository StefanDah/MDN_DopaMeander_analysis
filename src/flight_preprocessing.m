%==========================================================================
% File Name: flight_preprocessing.m
%
% Description:
%   Detect flight bouts and save spike rates, membrane potential etc.
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
cd(filepath)
addpath(genpath(filepath));

%% Introducing Global Variables and predifined Paramteres
WingbeatChanName='WingBeat1';   % Str included in Channel name of Wingbeat Trace
CurrentChanName='I_MTest 1';    % Str included in Channel name of Voltage Trace
checkbox_safespike =1;          % Recommended value is '1' when 'checkbox_diff' = 0. When turned on, the detection of multiple spikes independend of 'refractory time' is omitted
checkbox_spikecheck =1;         % Exemplary shows three excerpts of the RAW voltage trace with overlaying detected spikes
input_checklength =30;          % Number of seconds used in 'spikecheck' to evaluate the spike detection. Automatically reduced in shorter trials
input_meanspike = 30;           % Datasets with spike amplitude less than this value might not be saved (promt ask whether or not)
checkbox_manual_evaluation = 1;
checkbox_meanspikeshape =1;     % Shows the avarage spike shape
nr_spikes_avg=30;               % Number of avaraged spikes during spike shape check

checkbox_diff=0;                % Recommended value is '0'. If set to 0, threshold is determined on RAW data. If set to 1, threshold is determined on 'diff' voltage trace
if checkbox_diff==0
    refractory=0.020;
else
    refractory=0.002;
end
checkbox_show_dif_and_raw = 1;
checkbox_flyID = 1;             % If value is '1', reads fly data from excel spreadsheet
if checkbox_flyID ==1
    %Read Fly Data from excel sheet (FlyID)
    [FileID, FlyIDENT, CellID, Trial] = fcn_Excelimport("X:\MATLAB\data_flight\paper_2025\Fly_ID.xlsx", "Sheet1", [1, 54]); %last input is the number of flies in Excel sheet
    
    %Shorten These Variables in case empty values need to be excluded
    FileID(find(isnan(Trial)==1)) = [];
    FlyIDENT(find(isnan(Trial)==1)) = [];
    CellID(find(isnan(Trial)==1)) = [];
    Trial(find(isnan(Trial)==1)) = [];
end
checkbox_truncate = 1;
checkbox_truncate_no=0;
checkbox_truncate_yes=0;
checkbox_truncate_split=0;


spiking=1;                      %set to 0 for nonspiker
tacho=[];
current=[];
medi_intra=[];                  %%?
idx_spike_missmatch = 0;        % Counter for falsy detected spikes
spiketimes = {};


%NEW 2021_08_27 For Thresholding and Analysis
analysis_window_pre_onset_in_sec = 2; % Major analysis input, determines the area in which the spikes are counted and the acctivity assesed before an LED puls
analysis_window_post_onset_in_sec = 2;
analysis_window_pre_offset_in_sec = 2;
analysis_window_post_offset_in_sec = 2;

%Flight Paramters
min_flight_freq = 300;              % when set to 150: is 133Hz at 20.000 Sampling Rate, 300 is 66 Hz
min_flight_duration_in_sec = 4;
artifact_during_flight = 1500;      % loss of signal during one cont. flight bout tollerance: 1000 is 50 ms at 20000 sampling rate
min_interflight_interval = 200000;  % only takes flight bouts into account in which 15s (when set to 300000 at sampling rate 20000) no flight is detected for and after one bout
max_preflight_activity_in_sec = 0.5; %Wingbeats with this durations are allowed before a trials is excluded (not "clean")

idx_excluded_flight_bouts = 0;

% input_time_prior_flight = 5*sampling_rate; % in seconds*sampling rate Bin Window arround flight onset
% input_time_post_flight = 15*sampling_rate; % Bin Window arround flight offset



%% Load Files
% Folder with Files
files=dir([filepath '\20*.abf']);
% files = dir([filepath '\2021_01_13_0005.abf'])
if isempty(files) == 1
    error('No files found')
end

% Sort out files that have already been analysed
files_aready_analysed = dir([filepath '\20*.mat']);
index_lines_deleted = [];
mean_baseline_previouse_analysis = [];
mean_spike_ampl_previouse_analysis = [];
datafile_previouse_analysis = [];
idx_previouse_analyis = 0; %later used as index to save some data
if length(files_aready_analysed) > 0
    dlgTitle    = 'Analyse data again?';
    dlgQuestion = [num2str(length(files_aready_analysed)) ' Analysed Datasets have been found. Analyse those again?'];
    choice = questdlg(dlgQuestion,dlgTitle,'NO','YES', 'NO');
    if length(choice)==2 %NO
        files_without_old = files;
        for i=1:length(files_aready_analysed)
            tmp = files_aready_analysed(i).name;
            for z=1:length(files)
                if isequal(files(z).name(1:15), tmp(1:15)) == 1
                    index_lines_deleted = [index_lines_deleted z]; %delte entry from file so that it is not analyzed
                    eval(['load ' tmp])%load old datafile and save mean spikamplitude and baseline into a tempory file, later merged with data from this analysis
                    mean_baseline_previouse_analysis(idx_previouse_analyis+1) = mean_baseline;
                    mean_spike_ampl_previouse_analysis(idx_previouse_analyis+1) = spike_amplitude;
                    datafile_previouse_analysis{idx_previouse_analyis+1} = [str2mat(tmp(1:end-4))];
                    current_previouse_analysis(idx_previouse_analyis+1) = current;
                    idx_previouse_analyis = idx_previouse_analyis +1;
                end
            end
            
        end
        files(index_lines_deleted) = [];
    end
end
if isempty(files) == 1
    error('No files found')
end

dlgTitle    = 'Retrieve Threshold and trunc information from existing files?';
        dlgQuestion = 'Retrieve Threshold and trunc information from existing file';
        choice = questdlg(dlgQuestion,dlgTitle, 'NO', 'YES', 'YES');
        if length(choice)==2 %NO
            checkbox_retrieve_results_prev_analysis = 0;
            skip_file = 0;
        else
            checkbox_retrieve_results_prev_analysis = 1;
            skip_file = 0;
        end



for fly=1:length(files)
    
    %% Retrieve Threshold and trunc information from already
    %existing Mat Files
    if isempty(files) == 0 && checkbox_retrieve_results_prev_analysis == 1
        dlgTitle    = 'Retrieve Threshold and trunc information from existing files?';
        dlgQuestion = ['Retrieve Threshold and trunc information from existing file: ' num2str(files(fly).name(1:end-4)) ' ?'];
        choice = questdlg(dlgQuestion,dlgTitle,'NO, do again','SKIP', 'YES', 'YES');
        if length(choice)==12 %NO
            retrieve_results_prev_analysis = 0;
            skip_file = 0;
        elseif length(choice)==4
            skip_file = 1;
        el
            retrieve_results_prev_analysis = 1;
            skip_file = 0;
            
        end
    else
        retrieve_results_prev_analysis = 0;
    end
    
    if skip_file == 0
    file=num2str(files(fly).name(1:end-4)); % File name without '.abf'
    [d1,h_1]=abfload_Sander(strcat(filepath, file, '.abf'));
    sampling_rate = 1/(h_1.si/1000000); %Sampling Rate in Hz from ABF FIle(stored as microseconds)
    
    if retrieve_results_prev_analysis == 1 && strcmp(files(fly).name(1:15), files_aready_analysed(fly).name(1:15)) %In case file has already been analysed and we want to retrieve flight treshold and trunctimes
        mat_filename = files_aready_analysed(fly).name;
        matObj = matfile(mat_filename); %since these files tend to be very big, create a handle to that file to access only certain variables
        flightthreshold = matObj.flightthreshold;
        flight = matObj.flight;
        intra = matObj.intra;
        intra_orig = matObj.intra_orig;
        sampling_rate = matObj.sampling_rate;
        truncate_end = matObj.truncate_end;
        truncate_start = matObj.truncate_start;
        filename_modifier = matObj.filename_modifier;
        current = matObj.current;
        flyID = matObj.flyID;
        tacho = matObj.tacho;
        current_File = matObj.current_File;
        if strcmp(filename_modifier, '_trunc') %if data was truncated
            filename_modifier = ['_trunc'];
            idx_split = 1;
        elseif isempty(filename_modifier) %not truncated or splitted
            filename_modifier=[];
            idx_split = 1;
            %NOT YET PROGRAMMED FOR SPLITTED DATA since I did not use them so
            %far
        end
    else %Do the analysis again or for the first time in case there is no such *mat file
        
        %% Find Wingbeat and Current Chan
        idx_wingchan=strcmp(h_1.recChNames, WingbeatChanName);
        tachoID=find(idx_wingchan==1); % Determines Wingbeat Channel from all Channels
        
        if sum(strcmp(h_1.recChNames, CurrentChanName))==0 %In case current trace has another name (old setup)
            CurrentChanName='I_MTest';
        end
        if sum(strcmp(h_1.recChNames, CurrentChanName))==0
            CurrentChanName='I_MTest 1';
        end
        idx_currentchan=strcmp(h_1.recChNames, CurrentChanName);
        currentID=find(idx_currentchan==1);% Determines Current Channel from all Channels
        
        %% Load Wingbeat Current and Intra Channel
        % current
        current=d1(:,currentID);
        % Intra
        intra_orig=d1(:,1);
        % Wingbeat
        tacho = d1(:,tachoID);
        
        %% Truncate Dataset?
        if checkbox_truncate==1
            truncfig=figure('name','Truncate Dataset');
            set(truncfig, 'position', [1, 600, 1900, 450]);
            hold on
            subplot(2,1,1);
            plot(intra_orig,'k');
            
            subplot(2,1,2);
            plot(tacho, 'k');
            hold off
            uiwait
            
            dlgTitle    = 'Truncate or Split Trials?';
            dlgQuestion = 'Do you wish to truncate or split the trial?';
            choice = questdlg(dlgQuestion,dlgTitle,'No','Truncate', 'Split', 'No');
            if length(choice)==2
                checkbox_truncate_no=1;
                truncate_start = 1;
                truncate_end = length(intra_orig);
            elseif length(choice)==8
                checkbox_truncate_yes=1;
            elseif length(choice)==5
                checkbox_truncate_split=1;
            end
            
            if checkbox_truncate_yes==1
                
                truncfig=figure('name','Truncate Trial');
                set(truncfig, 'position', [1, 600, 1900, 450]);
                
                hold on
                subplot(2,1,1);
                plot(intra_orig,'k');
                
                subplot(2,1,2);
                plot(tacho, 'k');
                hold off
                ax = gca;
                lim = axis;
                [truncate_time,y] = ginput(2);
                truncate_start = round(truncate_time(1,1));
                if truncate_start < 1
                    truncate_start = 1;
                end
                truncate_end = round(truncate_time(2,1));
                if truncate_end > length(intra_orig)
                    truncate_end = length(intra_orig);
                end
                close(gcf);
                pause(1);
                current = current(truncate_start:truncate_end);
                intra_orig = intra_orig(truncate_start:truncate_end);
                tacho = tacho(truncate_start:truncate_end);
                
                filename_modifier = ['_trunc'];
                idx_split = 1;
                checkbox_truncate_yes = 0;
            elseif checkbox_truncate_no==1
                filename_modifier=[];
                idx_split = 1;
            elseif checkbox_truncate_split == 1
                truncfig=figure('name','Truncate Trial');
                set(truncfig, 'position', [1, 600, 1900, 450]);
                
                hold on
                subplot(2,1,1);
                plot(intra_orig,'k');
                
                subplot(2,1,2);
                plot(tacho, 'k');
                hold off
                ax = gca;
                lim = axis;
                [truncate_time,y] = ginput(4);
                truncate_start1=round(truncate_time(1,1));
                truncate_end1=round(truncate_time(2,1));
                truncate_start2=round(truncate_time(3,1));
                truncate_end2=round(truncate_time(4,1));
                close(gcf);
                pause(1);
                current2=current(truncate_start2:truncate_end2);
                intra_orig2=intra_orig(truncate_start2:truncate_end2);
                tacho2=tacho(truncate_start2:truncate_end2);
                
                current=current(truncate_start1:truncate_end1);
                intra_orig=intra_orig(truncate_start1:truncate_end1);
                tacho=tacho(truncate_start1:truncate_end1);
                
                filename_modifier=['_part1'];
                idx_split = 2;
                
            else
                truncate_start = 1;
                truncate_end = length(intra_orig);
            end
            
            
        end
        
    end
    
    for z=1:idx_split
        %% If Data are splitted
        if z==2
            current=current2;
            intra_orig=intra_orig2;
            tacho=tacho2;
            filename_modifier=['_part2'];
        end
        
        
        
        
        %% Process Intra (smooth)
        intra_smooth=smooth(intra_orig,50, 'loess'); %7 is good
        intra=[intra_smooth];
        medi_intra=[medi_intra, median(intra)];
        
        % IF THESE DATA ARE NOT ALREADY PRESENT FROM PREVIOUSE ANALYSIS
        if retrieve_results_prev_analysis == 1 && strcmp(files(fly).name(1:15), files_aready_analysed(fly).name(1:15))
            %Do nothing since already loaded
        else
            
            %% Current Fly from Excel Sheet
            if checkbox_flyID==1
                idx_flyID = find(FileID==file);
                current_File = FileID(idx_flyID);
                current_FlyNumber = FlyIDENT(idx_flyID);
                current_CellNumber = CellID(idx_flyID);
                current_Trial = Trial(idx_flyID);
                if idx_split ==1
                    trialID = 1;
                elseif idx_split ==2
                    trialID = 2;
                end
                flyID = strcat(['FLY', num2str(current_FlyNumber), '_CELL', num2str(current_CellNumber), '_TRIAL', num2str(current_Trial)]);
                disp(['Current FlyID is:' flyID])
            else
                current_File = file;
                flyID = ['NO_ID'];
                temp_fly = files(fly).name;
                save([num2str(files(fly).name(1:end-4)), filename_modifier, '.mat'],'flight', 'intra','puff', 'spikes', 'sampling_rate', 'current', 'current_File', 'flyID', 'tacho', 'mean_baseline', 'spike_amplitude', 'filename_modifier', 'spike_thresh');
                
            end
            
            %% Determine Tach Threshold with ginput
            
            if length(intra)/sampling_rate > 800 %in case its a very long trial so its depicted in two subplots
                tachofig=figure('name','Truncate Trial');
                set(tachofig, 'position', [1, 100, 1900, 900]);
                sp1 = subplot(2,1,1);
                x = 1:length(tacho)/2;
                x = x/sampling_rate;
                y = tacho(1:length(tacho)/2);
                plot(x, y)
                axis tight
                sp2 = subplot(2,1,2);
                x = length(tacho)/2:length(tacho);
                x = x/sampling_rate;
                y = tacho(length(tacho)/2:length(tacho));
                plot(x, y)
                axis tight
            else
                tachofig=figure('name','Truncate Trial');
                set(tachofig, 'position', [1, 600, 1900, 450]);
                plot(tacho);
            end
            
            flightthreshold=ginput(1);
            flightthreshold=flightthreshold(2);
            close(gcf);
            pause(1);
            
            %***** wing beat analysis
            beat_tmp=zeros(size(tacho));
            beat_tmp(tacho<=flightthreshold)=1;
            beat_tmp=diff(beat_tmp);
            
            beat_tmp2=zeros(size(beat_tmp));
            beat_tmp2(beat_tmp==1)=1;
            clear beat_tmp
            
            beattimes=find(beat_tmp2==1);
            beatfreq1=diff(beattimes);
            beatfreq=nan(size(beat_tmp2));
            
            flight=zeros(size(intra));
            
            for k=1:length(beattimes)-1
                beatfreq(beattimes(k))=beatfreq1(k);
                flight(beattimes(k))=1;
            end
            
            
        end
        
        % Determine flight on and offsets of "clean" flight
        % bouts (meaning no flight activity prior to start or post end of
        % flight (analysis_window_pre_onset_in_sec...)
        
        analysis_window_pre_onset = analysis_window_pre_onset_in_sec*sampling_rate; % Major analysis input, determines the area in which the spikes are counted and the acctivity assesed before an LED puls
        analysis_window_post_onset = analysis_window_post_onset_in_sec*sampling_rate;
        analysis_window_pre_offset = analysis_window_pre_offset_in_sec*sampling_rate;
        analysis_window_post_offset = analysis_window_post_offset_in_sec*sampling_rate;
        min_flight_duration = sampling_rate*min_flight_duration_in_sec;
        
        
        
        % To Determine the exact onset and offset: 'clean_flight' is generated
        clean_flight=zeros(size(flight));
        
        for i=1:length(flight)-min_flight_freq
            if  flight(i)==1
                if isempty(find(flight(i:i+min_flight_freq)==1))==0             %In case there is another wingbeat withing min Flight Freq (min_flight_freq)
                    last_wingbeat=max(find(flight(i:i+min_flight_freq)==1));
                    clean_flight(i:i+last_wingbeat)=1;
                else % In case this is the last wingbeat
                    clean_flight(i)=1;
                end
            end
        end
        
        %In case the animal was already flying in the beginning or is
        %still flying in the end
        if sum(clean_flight(1:min_flight_freq),1) > 1
            clean_flight(2:min_flight_freq) = 1;
            clean_flight(1) = 0;
        end
        if sum(clean_flight(end-min_flight_freq:end),1) > 1
            clean_flight(end-min_flight_freq:end-1) = 1;
            clean_flight(end) = 0;
        end
        
        orig_flight_onset = find(diff(clean_flight)==1);
        orig_flight_offset = find(diff(~clean_flight)==1);
        
        %             %In case the fly was already flying in the beginning -> exlude
        %             %this on/offset
        %             if orig_flight_onset(1) > orig_flight_offset(1)
        %                 orig_flight_offset(1) = [];
        %             end
        %             if orig_flight_onset(end) > orig_flight_offset(end)
        %                orig_flight_onset(end) = [];
        %             end
        
        % ************* Search for gaps/disturbances in the flight tacho trace that leads to detection of wrong onset/offsets
        %Backup on and offset prior to analysis
        if artifact_during_flight > min_flight_duration
            artifact_during_flight = min_flight_duration;
            disp('WARNING: artifact time is bigger than min_flight_duration')
        end
        
        
        
        %Finding Gaps and Artifacts
        temp_flight_onset = orig_flight_onset;
        temp_flight_offset = orig_flight_offset;
        clean_flight_without_artifacts = clean_flight;
        marker_artifact_between_bouts = nan(length(clean_flight_without_artifacts),1); % Used in Plot to indicate position of these artifacts
        idx_artifacts_detected_in_flight = 0;
        
        %             if length(temp_flight_offset)<length(temp_flight_onset) %In case fly is still flying in the end of rec
        %                 temp=length(temp_flight_onset);
        %                 temp_flight_onset(temp)=length(tacho);
        %             end
        
        for i=1:length(temp_flight_offset)
            if i==length(temp_flight_offset)
                %???
            elseif abs(temp_flight_offset(i)-temp_flight_onset(i+1))<artifact_during_flight
                
                for z = temp_flight_offset(i):temp_flight_onset(i+1)
                    clean_flight_without_artifacts(z) = 1;
                    marker_artifact_between_bouts(z) = 1;                                     % Used in Plot to indicate position of these artifacts
                end
                idx_artifacts_detected_in_flight = idx_artifacts_detected_in_flight+1;
            end
        end
        
        % ************* Search for Bouts that are shorter than 'min_flight_duration', by default 90s
        clean_flight_without_short_bouts=clean_flight_without_artifacts;
        marker_to_short_bouts = nan(length(clean_flight_without_short_bouts),1); % Used as marker in plot to highlight bouts that are to short
        
        %New flight Off- and Onsets are defined because removing artifact stiched
        %flight bouts together
        tmp4 = ~clean_flight_without_short_bouts;
        temp_flight_offset = find(diff(tmp4)==1);
        tmp5 = clean_flight_without_short_bouts;
        temp_flight_onset = find(diff(tmp5)==1);
        idx_to_short_flight_bouts = 0;
        
        for i=1:length(temp_flight_offset)
            
            if temp_flight_offset(i)-temp_flight_onset(i)<min_flight_duration && temp_flight_offset(i)-temp_flight_onset(i) ~= 0
                
                for z =temp_flight_onset(i):temp_flight_offset(i)
                    
                    clean_flight_without_short_bouts(z) = 0;
                    marker_to_short_bouts(z) =1; % Used as marker in plot to highlight bouts that are to short
                    
                end
                idx_to_short_flight_bouts = idx_to_short_flight_bouts+1;
            end
        end
        
        % ************* Search for flight bouts that are to close to each other ('min_interflight_interval')
        
        clean_flights_without_short_bouts_and_artefacts = clean_flight_without_short_bouts;
        marker_to_short_interval = nan(length(clean_flights_without_short_bouts_and_artefacts),1); %Used as marker in plot to highlight bouts that are to close to each other and ed
        
        temp_flight_onset = find(diff(clean_flights_without_short_bouts_and_artefacts)==1);
        tmp3 = ~clean_flights_without_short_bouts_and_artefacts;
        temp_flight_offset = find(diff(tmp3)==1);
        
        
        for i=1:length(temp_flight_offset)
            
            if i==length(temp_flight_offset)
                
                %???
                
            elseif temp_flight_onset(i+1)-temp_flight_offset(i)<min_interflight_interval
                
                for z =temp_flight_onset(i):temp_flight_offset(i)
                    
                    clean_flights_without_short_bouts_and_artefacts(z) = 0;
                    marker_to_short_interval(z) = 1; % Used as marker in plot to highlight bouts that are to close to each other and discarded
                    
                end
                idx_excluded_flight_bouts = idx_excluded_flight_bouts+1;
            end
        end
        
        %Define new on and offsets after excluding
        temp_flight_onset = find(diff(clean_flights_without_short_bouts_and_artefacts)==1);
        tmp6 = ~clean_flights_without_short_bouts_and_artefacts;
        temp_flight_offset = find(diff(tmp6)==1);
        
        %Sort out flight bouts at the beginning and end of rec that do
        %not allow for analysis with analysis_window_pre_onset..
        idx_excluded_flight_bouts_analysis_window = 0;
        
        delete_temp_flight_onset = [];
        delete_temp_flight_offset = [];
        for i = 1:length(temp_flight_onset)
            if temp_flight_onset(i)-analysis_window_pre_onset <= 1
                delete_temp_flight_onset(i) = i;
                delete_temp_flight_offset(i) = i;
                idx_excluded_flight_bouts_analysis_window = idx_excluded_flight_bouts_analysis_window + 1;
            end
        end
        temp_flight_onset(find(delete_temp_flight_onset>0)) = [];
        temp_flight_offset(find(delete_temp_flight_offset>0)) = [];
        
        delete_temp_flight_onset = [];
        delete_temp_flight_offset = [];
        for i = 1:length(temp_flight_offset)
            if temp_flight_offset(i)+analysis_window_post_offset >= length(intra)
                delete_temp_flight_onset(i) = i;
                delete_temp_flight_offset(i) = i;
                idx_excluded_flight_bouts_analysis_window = idx_excluded_flight_bouts_analysis_window + 1;
            end
        end
        temp_flight_onset(find(delete_temp_flight_onset>0)) = [];
        temp_flight_offset(find(delete_temp_flight_offset>0)) = [];
        
        %Sort out flight bouts that are not "clean" depending on
        %analysis_window_pre_onset_in_sec..
        delete_temp_flight_onset = [];
        delete_temp_flight_offset = [];
        for i = 1:length(temp_flight_onset)
            if (sum(clean_flight_without_artifacts(temp_flight_onset(i)-analysis_window_pre_onset:temp_flight_onset(i)))/sampling_rate) >= max_preflight_activity_in_sec ... %for onset
                    || (sum(clean_flight_without_artifacts(temp_flight_offset(i):temp_flight_offset(i)+analysis_window_post_offset))/sampling_rate) >= max_preflight_activity_in_sec %for offset
                delete_temp_flight_onset(i) = i;
                delete_temp_flight_offset(i) = i;
                idx_excluded_flight_bouts = idx_excluded_flight_bouts + 1;
            end
        end
        
        temp_flight_onset(find(delete_temp_flight_onset>0)) = [];
        temp_flight_offset(find(delete_temp_flight_offset>0)) = [];
        
        %Delete these on and offsets from clean_flights_without_short_bouts_and_artefacts
        clean_flights_without_short_bouts_and_artefacts = zeros(length(clean_flights_without_short_bouts_and_artefacts),1);
        for i = 1:length(temp_flight_onset)
            clean_flights_without_short_bouts_and_artefacts(temp_flight_onset(i):temp_flight_offset(i)) = 1;
            
            %                 if clean_flights_without_short_bouts_and_artefacts(temp_flight_onset(i)+1) == 1 || clean_flights_without_short_bouts_and_artefacts(temp_flight_offset(i)+1) == 1
            %                    clean_flights_without_short_bouts_and_artefacts(temp_flight_onset(i):temp_flight_offset(i)) = 0;
            %                 end
        end
        
        
        %Display how many bouts have been excludede due to interflight interval
        disp([flyID ' from ' filepath files(fly).name])% Dsip name of the file
        if idx_excluded_flight_bouts > 0
            disp([num2str(idx_excluded_flight_bouts), ' Flight Trials have been exluded due to insufficient interflight intervals (magenta)'])
        end
        if idx_excluded_flight_bouts_analysis_window > 0
            disp([num2str(idx_excluded_flight_bouts_analysis_window), ' Flight Trials have been exluded due to analysis window conflict at the beginning and end of a trial'])
        end
        if idx_to_short_flight_bouts > 0
            disp([num2str(idx_to_short_flight_bouts), ' Flight Trials have been exluded because they are to short (red)'])
        end
        if idx_artifacts_detected_in_flight > 0
            disp([num2str(idx_artifacts_detected_in_flight), ' Artifacts have been idetified and corrected wihtin a flight bout (green)'])
        end
        
        
        %NEW 2021_10_12 IF IN PREVIOUSE ANALYSIS ON AND OFFSETS WERE DETERMINED PRINT
        %THEM HERE:
        if retrieve_results_prev_analysis == 1 && strcmp(files(fly).name(1:15), files_aready_analysed(fly).name(1:15)) %In case file has already been analysed and we want to retrieve flight treshold and trunctimes
            previouse_flight_onset = matObj.flight_onset;
            previouse_flight_offset = matObj.flight_offset;
            disp([num2str(length(previouse_flight_onset)) ' clean flight bouts detected from prev. analysis'])
            for z1 = 1:length(previouse_flight_onset)
                disp([num2str(z1) ' flight onset was at: ' num2str(previouse_flight_onset(z1)/sampling_rate) ' s'])
                disp([num2str(z1) ' flight offset was at: ' num2str(previouse_flight_offset(z1)/sampling_rate) ' s'])
                disp([])
            end
        end
        
        
        
        %% ************* Plotting Results *************
        if checkbox_manual_evaluation == 1 %ONLY EXECUTED WHEN SET TO FULLY AUTOMATED ANALYSIS
            %*Determine X-Achsis in s
            x_axis = (1:length(tacho))/sampling_rate;
            tachocheck=figure('name','Tach Analysis');
            set(tachocheck, 'position', [1, 200, 1900, 800]);
            hold on
            
            sp1 = subplot(2,1,1);
            plot(x_axis, tacho, 'b')
            hold on
            plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
            plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
            plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
            hold off
            
            sp2 = subplot(2,1,2);
            linkaxes([sp1, sp2], 'x')
            plot(x_axis, clean_flights_without_short_bouts_and_artefacts,'b')
            ylim([0.99,1.01])
            hold on
            plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
            plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
            plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
            hold off
            
            %% ************* Correct Results from Analysis manually, if needed *************
            uiwait % wait until figure is closed
            indicator_analysis_correct = 0;
            indicator_discard = 0;
            indicator_add = 0;
            idx_discard_onset = 0;
            idx_discard_offset = 0;
            
            while indicator_analysis_correct == 0
                
                % Promt to ask whether the bouts shall be deleted or included in analysis
                dlgTitle    = 'Discard Trials or Include Trials in Analysis?';
                dlgQuestion = 'Do you wish to DISCARD or ADD Trials in Analysis?';
                choice = questdlg(dlgQuestion,dlgTitle,'Analysis is correct','DISCARD', 'ADD', 'Analysis is correct');
                if length(choice)==19  %Analysis is correct
                    indicator_analysis_correct=1;
                elseif length(choice)==7 %Discard Trials
                    indicator_discard=1;
                elseif length(choice)==3 %Add Trials
                    indicator_add=1;
                end
                
                %BACKUP FOR RETURN
                clean_flights_without_short_bouts_and_artefacts_backup = clean_flights_without_short_bouts_and_artefacts;
                
                %******************** DISCARDING TRIALS**************************%
                if indicator_discard == 1
                    %Plotting
                    x_axis = (1:length(tacho))/sampling_rate;
                    tachocheck2=figure('name','Discard Trials from Tacho Analysis');
                    set(tachocheck2, 'position', [1, 200, 1900, 800]);
                    hold on
                    sp1 = subplot(2,1,1);
                    plot(x_axis, tacho, 'b')
                    hold on
                    plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
                    plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
                    plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
                    hold off
                    
                    sp2 = subplot(2,1,2);
                    linkaxes([sp1, sp2], 'x');
                    
                    plot(x_axis, clean_flights_without_short_bouts_and_artefacts,'b')
                    ylim([0.99,1.01])
                    hold on
                    
                    plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
                    plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
                    plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
                    
                    %Select Trials by clicking at the beginning and end of a flight bout,  return key for
                    %confirm
                    [idx_discard,y] = ginput;
                    idx_discard=idx_discard*sampling_rate;
                    close(gcf);
                    pause(1);
                    
                    %Convert input into two vectors
                    for z =1:length(idx_discard)
                        if floor(z/2)==z/2 % code for even = offset
                            idx_discard_offset(z-1) = idx_discard(z);
                            
                        else % code for odd = onset
                            idx_discard_onset(z) = idx_discard(z);
                            
                        end
                        
                    end
                    idx_discard_onset =idx_discard_onset(idx_discard_onset~=0); %remove zeros
                    idx_discard_offset =idx_discard_offset(idx_discard_offset~=0); %remove zeros
                    
                    %Use Ginput coordinates to find trials in tacho analysis
                    %(in: clean_flights_without_short_bouts_and_artefacts)
                    
                    for z = 1:length(idx_discard_offset)
                        % find nearest datapoint in "onset"
                        [minValue, closestIndex] = min(abs(temp_flight_onset - idx_discard_onset(z).'));
                        delete_onset = temp_flight_onset(closestIndex); %closestValue in onset
                        
                        % find nearest datapoint in "offset"
                        [minValue, closestIndex] = min(abs(temp_flight_offset - idx_discard_offset(z).'));
                        delete_offset = temp_flight_offset(closestIndex); %closestValue in offset
                        
                        clean_flights_without_short_bouts_and_artefacts(delete_onset:delete_offset)=0;
                        
                        discard_onset_to_accepted(z) = delete_onset;
                        discard_offset_to_accepted(z) = delete_offset;
                        
                        
                    end
                    %Summarizing Results
                    disp([num2str(z), ' Flight Bouts have been excluded: '])
                    for u = 1:length(idx_discard)
                        if floor(u/2)==u/2 % code for even = offset
                            disp([num2str(idx_discard(u)/sampling_rate), 's & ']);
                            
                        else % code for odd = onset
                            disp([num2str(idx_discard(u)/sampling_rate), 's - ']);
                            
                        end
                    end
                    
                    %Clearing Workspace
                    idx_discard_onset = 0;
                    idx_discard_offset = 0;
                    clearvars minValue closestIndex delete_onset delete_offset idx_discard
                    
                    % CORRECT OR RETURN TO PREVIOUSE RESULTS
                    %Plotting
                    x_axis = (1:length(clean_flights_without_short_bouts_and_artefacts_backup))/sampling_rate;
                    tachocheck3=figure('name','Corrections made correct?');
                    set(tachocheck3, 'position', [1, 200, 1900, 800]);
                    hold on
                    sp1 = subplot(2,1,1);
                    plot(x_axis, clean_flights_without_short_bouts_and_artefacts_backup, 'b')
                    ylim([0.99,1.01])
                    hold on
                    plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
                    plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
                    plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
                    hold off
                    sp2 = subplot(2,1,2);
                    linkaxes([sp1, sp2], 'x');
                    plot(x_axis, clean_flights_without_short_bouts_and_artefacts,'b')
                    ylim([0.99,1.01])
                    hold on
                    plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
                    plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
                    plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
                    
                    dlgTitle    = 'Corrections made are correct?';
                    dlgQuestion = 'Accept corrections made or undo changes?';
                    choice = questdlg(dlgQuestion,dlgTitle,'ACCEPT','UNDO', 'ACCEPT');
                    if length(choice)==6  %ACCEPT
                        clearvars clean_flights_without_short_bouts_and_artefacts_backup
                        temp_flight_onset = setdiff(temp_flight_onset, discard_onset_to_accepted');
                        temp_flight_offset = setdiff(temp_flight_offset, discard_offset_to_accepted');
                        close(gcf)
                        indicator_discard = 0;
                        
                        %Save this information in the original *mat File so that these corrections are automatically done in the next run of the analysis?
                        dlgTitle    = 'SAVE ANALYSIS TO ORIGINAL MAT?';
                        dlgQuestion = 'Save this information in the original *mat File so that these corrections are automatically done in the next run of the analysis?';
                        choice = questdlg(dlgQuestion,dlgTitle,'SAVE','DONT SAVE', 'SAVE');
                        if length(choice)==4  %SAVE
                            flight_onset_trigger =temp_flight_onset;
                            flight_offset_trigger = temp_flight_offset;
                            save([num2str(files(fly).name(1:end-4)), '.mat'],'flight', 'intra','puff', 'spikes', 'sampling_rate', 'current', 'current_File', 'flyID', 'tacho', 'mean_baseline', 'spike_amplitude', 'clean_flights_without_short_bouts_and_artefacts', 'flight_onset_trigger', 'flight_offset_trigger');
                        end
                        
                    elseif length(choice)==4 %UNDO
                        clearvars clean_flights_without_short_bouts_and_artefacts;
                        clean_flights_without_short_bouts_and_artefacts = clean_flights_without_short_bouts_and_artefacts_backup;
                        close(gcf)
                    end
                    
                    
                end
                %******************** ADD TRIALS**************************%
                
                if indicator_add == 1
                    
                    %Plotting
                    x_axis = (1:length(tacho))/sampling_rate;
                    tachocheck=figure('name','ADD Trials to Tacho Analysis');
                    set(tachocheck, 'position', [1, 200, 1900, 800]);
                    hold on
                    
                    sp1 = subplot(2,1,1);
                    plot(tacho, 'b')
                    hold on
                    plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
                    plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
                    plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
                    %                         hold off
                    
                    sp2 = subplot(2,1,2);
                    
                    plot(x_axis, clean_flights_without_short_bouts_and_artefacts,'b')
                    hold on
                    ylim([0.99,1.01])
                    %                         linkaxes([sp1, sp2], 'x');
                    
                    
                    plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
                    plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
                    plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
                    hold off
                    
                    %Select Trials by clicking at the beginning and end of a flight bout,  return key for
                    %confirm
                    [idx_add,y] = ginput;
                    idx_add=idx_add*sampling_rate;
                    close(gcf);
                    pause(1);
                    
                    %Convert input into two vectors
                    
                    for z =1:length(idx_add)
                        if floor(z/2)==z/2 % code for even = offset
                            idx_add_offset(z-1) = idx_add(z);
                            
                        else % code for odd = onset
                            idx_add_onset(z) = idx_add(z);
                            
                        end
                        
                    end
                    idx_add_onset =idx_add_onset(idx_add_onset~=0); %remove zeros
                    idx_add_offset =idx_add_offset(idx_add_offset~=0); %remove zeros
                    
                    %Use Ginput coordinates to find trials in tacho analysis
                    %(in: clean_flights_without_short_bouts_and_artefacts)
                    add_to_accepted = [];
                    
                    for z = 1:length(idx_add_offset)
                        % find nearest datapoint in "onset"
                        [minValue, closestIndex] = min(abs(orig_flight_onset - idx_add_onset(z).'));
                        add_onset = orig_flight_onset(closestIndex); %closestValue in onset
                        
                        % find nearest datapoint in "offset"
                        [minValue, closestIndex] = min(abs(orig_flight_offset - idx_add_offset(z).'));
                        add_offset = orig_flight_offset(closestIndex); %closestValue in offset
                        
                        clean_flights_without_short_bouts_and_artefacts(add_onset:add_offset)=1;
                        
                        add_onset_to_accepted(z) = add_onset;
                        add_offset_to_accepted(z) = add_offset;
                    end
                    
                    %Summarizing Results
                    disp([num2str(z), ' Flight Bouts have been added: '])
                    for u = 1:length(idx_add)
                        if floor(u/2)==u/2 % code for even = offset
                            disp([num2str(idx_add(u)/sampling_rate), 's & ']);
                            
                        else % code for odd = onset
                            disp([num2str(idx_add(u)/sampling_rate), 's - ']);
                            
                        end
                    end
                    
                    %Clearing Workspace
                    idx_add_onset = 0;
                    idx_add_offset = 0;
                    clearvars minValue closestIndex add_onset add_offset idx_add
                    
                    % CORRECT OR RETURN TO PREVIOUSE RESULTS
                    %Plotting
                    x_axis = (1:length(clean_flights_without_short_bouts_and_artefacts_backup))/sampling_rate;
                    tachocheck=figure('name','Corrections made correct?');
                    set(tachocheck, 'position', [1, 200, 1900, 800]);
                    hold on
                    sp1 = subplot(2,1,1);
                    plot(x_axis, clean_flights_without_short_bouts_and_artefacts_backup, 'b')
                    ylim([0.99,1.01])
                    hold on
                    plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
                    plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
                    plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
                    hold off
                    sp2 = subplot(2,1,2);
                    plot(x_axis, clean_flights_without_short_bouts_and_artefacts,'b')
                    ylim([0.99,1.01])
                    hold on
                    plot(x_axis, marker_to_short_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 0]);
                    plot(x_axis, marker_to_short_interval, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[1 0 1]);
                    plot(x_axis, marker_artifact_between_bouts, 'x', 'MarkerSize',10, 'MarkerEdgeColor',[0 1 0]);
                    linkaxes([sp1, sp2], 'x');
                    
                    dlgTitle    = 'Corrections made are correct?';
                    dlgQuestion = 'Accept corrections made or undo changes?';
                    choice = questdlg(dlgQuestion,dlgTitle,'ACCEPT','UNDO', 'ACCEPT');
                    if length(choice)==6  %ACCEPT
                        clearvars clean_flights_without_short_bouts_and_artefacts_backup
                        temp_flight_onset = union(temp_flight_onset, add_onset_to_accepted');
                        temp_flight_offset = union(temp_flight_offset, add_offset_to_accepted');
                        %NEW sort out onsets/offsets in between so that
                        %they have the same lengths
                        %NEW
                        idx_temp = length(temp_flight_onset)-length(temp_flight_offset); %if idx_temp is acchieved: both have equal size, continue
                        if length(temp_flight_onset) > length(temp_flight_offset)
                            disp('adjusting on and offsets');
                            while idx_temp ~= 0
                                for y=1:length(temp_flight_onset)
                                    if temp_flight_onset(y+1)<temp_flight_offset(y)
                                        temp_flight_onset_temp = y+1;
                                        idx_temp = idx_temp -1;
                                    end
                                    if idx_temp == 0
                                        break
                                    end
                                    
                                end
                            end
                            temp_flight_onset(temp_flight_onset_temp) = [];
                        end
                        %NEW END
                        
                        close(gcf)
                        indicator_add = 0;
                        
                        
                        %Save this information in the original *mat File so that these corrections are automatically done in the next run of the analysis?
                        dlgTitle    = 'SAVE ANALYSIS TO ORIGINAL MAT?';
                        dlgQuestion = 'Save this information in the original *mat File so that these corrections are automatically done in the next run of the analysis?';
                        choice = questdlg(dlgQuestion,dlgTitle,'SAVE','DONT SAVE', 'SAVE');
                        if length(choice)==4  %SAVE
                            flight_onset_trigger =temp_flight_onset;
                            flight_offset_trigger = temp_flight_offset;
                            save([num2str(files(fly).name(1:end-4)), '.mat'],'flight', 'intra','puff', 'spikes', 'sampling_rate', 'current', 'current_File', 'flyID', 'tacho', 'mean_baseline', 'spike_amplitude', 'clean_flights_without_short_bouts_and_artefacts', 'flight_onset_trigger', 'flight_offset_trigger');
                        end
                        
                        
                    elseif length(choice)==4 %UNDO
                        clearvars clean_flights_without_short_bouts_and_artefacts;
                        clean_flights_without_short_bouts_and_artefacts = clean_flights_without_short_bouts_and_artefacts_backup;
                        close(gcf)
                    end
                    
                    
                end
                
                
            end
        end
        
        
        %Producing final vectors
        flight_indicator = clean_flights_without_short_bouts_and_artefacts;
        flight_onset = temp_flight_onset;
        flight_offset = temp_flight_offset;
        
        
        
        %% Determine Spike Treshold with ginput
        
        checkbox_analysis_correct = 0;
        
        if spiking==1
            while checkbox_analysis_correct == 0 %Do it as long as checkbox_analysis_correct == 1
                
                if checkbox_show_dif_and_raw == 1
                    threshfig_raw=figure('name','RAW');
                    set(threshfig_raw, 'position', [1, 600, 1900, 450]);
                    plot(intra,'k');
                    threshfig_diff=figure('name','DIFF');
                    set(threshfig_diff, 'position', [0, 0, 1900, 450])
                    plot(diff(intra),'k');
                    
                    dlgTitle    = 'RAW or DIFF?';
                    dlgQuestion = 'Do you wish use RAW or DIFF to set the spike threshold?';
                    choice = questdlg(dlgQuestion,dlgTitle,'RAW','DIFF', 'RAW');
                    if length(choice)==3
                        checkbox_diff = 0;
                    else
                        checkbox_diff = 1;
                    end
                    
                    
                    close (gcf)
                    close (gcf)
                end
                
                threshfig=figure;
                set(threshfig, 'position', [1, 1, 1900, 1000])
                
                if checkbox_diff==1
                    plot(diff(intra),'k');
                else
                    plot(intra,'k');
                end
                set(threshfig, 'position', [1, 1, 2000, 1000])
                xlim([5*sampling_rate,150*sampling_rate])
                if checkbox_diff==1
                    ylim([0, 2])%max(diff(intra))])y
                end
                
                spike_thresh=ginput(1);
                spike_thresh=spike_thresh(1,2);
                
                close(threshfig);
                
                
                
                overthresh=zeros(size(intra));
                if checkbox_diff ==1
                    overthresh(diff(intra)>spike_thresh)=1;
                    tricky=diff(overthresh);
                else
                    overthresh(intra>spike_thresh)=1;
                    tricky=(overthresh);
                end
                
                
                spikes=nan(size(intra));
                spikes(tricky==1)=1;
                
                
                
                
                for s=3:length(spikes)
                    if spikes(s-1)==1
                        spikes(s:s+refractory*sampling_rate)=nan;
                    else
                    end
                end
                
                if checkbox_safespike ==1 %% Find and delete false positive spikes
                    for j=3:length(spikes)-2 %Does not take the first two and last two "spikes" into account
                        if spikes(j)==1
                            tmp_prespike=intra(j-1);
                            tmp_postspike=intra(j+1);
                            
                            if tmp_prespike>intra(j) && tmp_postspike<intra(j)
                                spikes(j)=nan;
                                idx_spike_missmatch = idx_spike_missmatch+1;
                            end
                            
                        end
                    end
                end
                
                if isempty(idx_spike_missmatch)&& checkbox_safespike ==1
                    disp('No Spikes missmatches detected')
                elseif checkbox_safespike ==1
                    disp([num2str(idx_spike_missmatch) ' Spikes missmatches detected'])
                end
                
                spiketimes=find(spikes==1);
                
                %% Checking Spikes
                if(checkbox_spikecheck==1)
                    
                    checklength=input_checklength*sampling_rate;
                    
                    if length(intra)/2-checklength < checklength
                        checklength=round(length(intra)/2-1);
                        disp('Checklength was automatically adjusted because the trial is too short')
                    end
                    
                    spikecheck=figure;
                    set(spikecheck, 'position', [0, 0, 1900, 450])
                    plot(intra(1:checklength), 'r')
                    hold all
                    plot(spikes(1:checklength)*medi_intra+10, '+r', 'MarkerSize', 15)
                    
                    plot(intra(length(intra)-checklength:(length(intra))), 'k')
                    plot(spikes(length(intra)-checklength:(length(intra)))*medi_intra+20, '+k', 'MarkerSize', 15)
                    
                    plot(intra(length(intra)/2-checklength:(length(intra)/2)), 'b')
                    plot(spikes(length(intra)/2-checklength:(length(intra)/2))*medi_intra+25, '+b', 'MarkerSize', 15)
                    
                    hold off
                    uiwait
                    % next, i want to see whether the spike sorting gave me something decent.
                    
                    dlgTitle    = 'Spike Sorting Correct?';
                    dlgQuestion = 'are you satisfied with the spike sorting, you pretty fuck?';
                    choice = questdlg(dlgQuestion,dlgTitle,'YES','NO', 'YES');
                    if length(choice)==2
                        %                     error('start over')
                        checkbox_analysis_correct = 0;
                        idx_spike_missmatch = 0;
                        clearvars thresh overthresh tricky spikes tmp_prespike tmp_postspike spiketimes
                    else
                        checkbox_analysis_correct = 1;
                    end
                    
                    
                    %here, i plot the average spike shape for the first bit, defined by
                    %nr_spikes_avg
                end
            end
            
            allspikes=[];
            meanspike=[];
            
            %% Plot mean spike shape
            if(checkbox_meanspikeshape==1)
                
                if length(spiketimes)<nr_spikes_avg
                    nr_spikes_avg = length(spiketimes);
                    disp(['Number of avaraged spikes was reduced to ' num2str(nr_spikes_avg) ' because there are not enought spikes'])
                end
                %Remove first spike if it is at the ebegining of the rec
                %and would result in -intra values
                if spiketimes(1)-0.02*sampling_rate<intra(1)
                    spiketimes(1) = [];
                end
                
                for i=2:nr_spikes_avg
                    allspikes(:,i-1)=intra(spiketimes(i)-0.02*sampling_rate:spiketimes(i)+0.02*sampling_rate);
                end
                meanspike=mean(allspikes');
                figure
                plot(meanspike, 'r');
                
            else
                disp('non-spiker');
                
            end
            
        end
        
        %% Print out means
        mean_spikerate=mean(diff(spiketimes))/sampling_rate;
        spike_interval=diff(spiketimes)/sampling_rate;
        mean_spike_interval=mean(spike_interval);
        mean_spike_rate=1/mean_spike_interval;
        
        disp(['The mean spike interval is ' num2str(mean_spike_interval) 'Hz'])
        disp(['The mean spike rate is ' num2str(mean_spike_rate) 'Hz'])
        
        %% Current Injected? How much?
        if exist('currentID')
            current_orig=d1(:,currentID);
            current=mean(current_orig);
        end
        disp([num2str(current) 'pA injected'])
        
        
        
        %% Calculating Baseline Resting Membrane potential at three different points
        region1 = median(intra(1:(floor((length(intra))/3))));
        region2 = median(intra((floor((length(intra))/3)):(floor((length(intra))/3)*2)));
        region3 = median(intra(floor(((length(intra))/3)*2):floor(((length(intra))/3)*3)));
        mean_baseline = mean([region1 region2 region3]);
        disp([num2str(mean_baseline) 'mV is the mean baseline activity']);
        disp([num2str(max(meanspike)-mean_baseline) 'mV is the average spike amplitude']);
        spike_amplitude = max(meanspike)-mean_baseline;
        
        
        %% Saving Results
        %First ask whether data shall be saved and if "meanspike" is larger
        %than the desired value (input_meanspike)
        if input_meanspike > max(meanspike)-mean_baseline
            dlgTitle    = 'Meanspike';
            disp([num2str(meanspike) 'mV is the mean spike-amplitude']);
            dlgQuestion = 'Save Analysis anyway?';
            choice = questdlg(dlgQuestion,dlgTitle,'YES','NO', 'YES');
            if length(choice)==2
                decission_meanspike = 0;
            else
                decission_meanspike = 1;
            end
        else
            decission_meanspike = 1;
        end
        
        if decission_meanspike == 1
            %*******Write Fly Data from excel sheet (FlyID) to mat file
            % Save the mean baseline and spike amplitude in a seperate
            % variable
            mean_baseline_all(fly) = mean_baseline;
            spike_amplitude_all(fly) = spike_amplitude;
            datafile{fly} = file;
            if isempty(current) ~= 1
                current_amplitude_all(fly) = current;
            else
                current_amplitude_all(fly) = 0;
            end
            
            %             if checkbox_flyID==1
            %                 idx_flyID = find(FileID==file);
            %                 current_File = FileID(idx_flyID);
            %                 current_FlyNumber = FlyIDENT(idx_flyID);
            %                 current_CellNumber = CellID(idx_flyID);
            %                 current_Trial = Trial(idx_flyID);
            %                 if idx_split ==1
            %                     trialID = 1;
            %                 elseif idx_split ==2
            %                     trialID = 2;
            %                 end
            %                 flyID = strcat(['FLY', num2str(current_FlyNumber), '_CELL', num2str(current_CellNumber), '_TRIAL', num2str(current_Trial)]);
            %                 disp(['Current FlyID is:' flyID])
            %% Save results for current dataset to mat file
            save([num2str(files(fly).name(1:end-4)), filename_modifier '.mat'],'flight', 'intra', 'spikes', 'sampling_rate', 'current', 'current_File', 'flyID', 'tacho', 'mean_baseline', 'spike_amplitude', 'filename_modifier', 'flight_indicator', 'flight_onset', 'flight_offset', 'flightthreshold', 'truncate_start', 'truncate_end', 'analysis_window_pre_onset', 'analysis_window_post_onset', 'analysis_window_pre_offset', 'analysis_window_post_offset', 'intra_orig', 'spike_thresh');
            
            %             else
            %                 current_File = file;
            %                 flyID = ['NO_ID'];
            %                 temp_fly = files(fly).name;
            %                 save([num2str(files(fly).name(1:end-4)), filename_modifier, '.mat'],'flight', 'intra','puff', 'spikes', 'sampling_rate', 'current', 'current_File', 'flyID', 'tacho', 'mean_baseline', 'spike_amplitude', 'filename_modifier');
            %             end
            
            if fly==length(files)
                disp('FINISH')
            end
        end
    end
end
end
%*******Write Fly Data from excel sheet (FlyID) to mat file
% Save the mean baseline and spike amplitude in a seperate
% variable

if isempty(mean_baseline_previouse_analysis) == 0
    tmp5 = nonzeros(mean_baseline_all);
    tmp6 = nonzeros(mean_baseline_previouse_analysis);
    clear mean_baseline_all
    mean_baseline_all = combine(tmp5, tmp6);
    
    tmp5 = nonzeros(spike_amplitude_all);
    tmp6 = nonzeros(mean_spike_ampl_previouse_analysis);
    clear spike_amplitude_all
    spike_amplitude_all = combine(tmp5, tmp6);
    
    
    tmp5 = nonzeros(current_amplitude_all);
    tmp6 = nonzeros(current_previouse_analysis);
    clear current_amplitude_all
    current_amplitude_all = combine(tmp5, tmp6);
    
    
    datafile = datafile(~cellfun('isempty',datafile));
    datafile_previouse_analysis = datafile_previouse_analysis(~cellfun('isempty',datafile_previouse_analysis));
    datafile = [datafile datafile_previouse_analysis] ;
    
end
save(['Analysis_SPIKE_MEAN' '.mat'],'datafile', 'mean_baseline_all', 'spike_amplitude_all', 'current_amplitude_all');

%% Histograms for SpikeAmp and Baseline
baseline_and_spike_figure = figure('Color', [0 0 0]);
set(baseline_and_spike_figure, 'position', [1, 100, 1900, 800]);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName','Times','fontsize',18);
hold on
sp1 = subplot(1,2,1);
hold on
hist(mean_baseline_all, 10);
ax_sp1 = gca;
set(gca,'color',[0 0 0]);
set(gcf,'inverthardcopy','off');
set(ax_sp1,'YColor','w', 'XColor','w')
set(gca,'TickDir','out');
% set(gca, 'Xcolor', 'k');
set(gcf,'inverthardcopy','off');
set(ax_sp1,'YColor','w', 'XColor','w')
set(gca,'TickDir','out');
ylabel('Number of Cells');
xlabel('Membrane Potential (mV)');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16);




h1 = findobj(gca,'Type','patch');
hold on
h1.FaceColor = [0.75 0.75 0.75];
h1.EdgeColor = 'w';

sp2 = subplot(1,2,2);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName','Times','fontsize',18);
hold on
hist(spike_amplitude_all, 10);
ax_sp2 = gca;
set(gca,'color',[0 0 0]);
set(gcf,'inverthardcopy','off');
set(ax_sp2,'YColor','w', 'XColor','w')
set(gca,'TickDir','out');
set(gca, 'Xcolor', 'k');
set(gcf,'inverthardcopy','off');
set(ax_sp2,'YColor','w', 'XColor','w')
set(gca,'TickDir','out');
ylabel('Number of Cells');
xlabel('Spike Amplitude (mV)');
b = get(gca,'XTickLabel');
set(gca,'XTickLabel',b,'FontName','Times','fontsize',16);


h2 = findobj(gca,'Type','patch');
h2.FaceColor = [0.75 0.75 0.75];
h2.EdgeColor = 'w';

saveas(gcf, 'Baseline_SpikeAmp.tif')
print(['Baseline_SpikeAmp.eps'], '-depsc2', '-tiff', '-r300', '-painters')
pause(1)
close gcf


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

fid = fopen('Analysis_Paramter_Preprocessing.txt','w');
fprintf(fid, '%6s %12s\r\n', 'input_checklength', num2str(input_checklength));
fprintf(fid, '%6s %12s\r\n', 'min_flight_duration_in_sec', num2str(min_flight_duration_in_sec));
fprintf(fid, '%6s %12s\r\n', 'min_flight_freq', num2str(min_flight_freq));
fprintf(fid, '%6s %12s\r\n', 'artifact_during_flight', num2str(artifact_during_flight));
fprintf(fid, '%6s %12s\r\n', 'min_interflight_interval', num2str(min_interflight_interval));
fprintf(fid, '%6s %12s\r\n', 'max_preflight_activity_in_sec', num2str(max_preflight_activity_in_sec));
fprintf(fid, '%6s %12s\r\n', 'analysis_window_pre_onset_in_sec', num2str(analysis_window_pre_onset_in_sec));
fprintf(fid, '%6s %12s\r\n', 'analysis_window_post_onset_in_sec', num2str(analysis_window_post_onset_in_sec));
fprintf(fid, '%6s %12s\r\n', 'analysis_window_pre_offset_in_sec', num2str(analysis_window_pre_offset_in_sec));
fprintf(fid, '%6s %12s\r\n', 'analysis_window_post_offset_in_sec', num2str(analysis_window_post_offset_in_sec));
fprintf(fid, '%6s %12s\r\n', 'scriptversion', scriptversion);

fclose(fid);
disp('FINISH')
