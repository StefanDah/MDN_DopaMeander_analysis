%% Load abf files and set parameters
filepath='path';
cd(filepath)
addpath(genpath(filepath));

%% List of channels that depend on the experiment
VMChanName = 'Vm_scaled';
ICurrentChanName = 'I_MTest 1';
CameraChanName = 'Camera1';
VelocXChanName = 'Veloc X';
VelocZChanName = 'Veloc Z';
VelocYChanName = 'Veloc Y';
WingBeatChanName = 'WingBeat1';
LEDChanName = 'DAC1';
LED2ChanName = 'DAC0';

smoothing_factor = 2;
nr_spikes_avg = 300;

%% Load data from flysheet 'Fly_ID'

[Date, FlyNumber, TrialNumber, CellNumber] = fcn_Excelimport('filepath_placeholder\Fly_ID.xlsx', 'Sheet1', [1,99]);

%Shorten These Variables in case empty values need to be excluded
Date(find(isnan(TrialNumber)==1)) = [];
FlyNumber(find(isnan(TrialNumber)==1)) = [];
CellNumber(find(isnan(TrialNumber)==1)) = [];
TrialNumber(find(isnan(TrialNumber)==1)) = [];


%% Load files and start analysis

files = dir([filepath '\*.abf']);  % Folder with .abf files


for fly = 1:length(files)

    file=num2str(files(fly).name(1:end-4)); % File name without '.abf' (-4 characters)

    [d1,h_1] = abfload_Sander(strcat(filepath, file, '.abf'));

    sampling_rate = 1/(h_1.si/1000000); %Sampling Rate in Hz from ABF FIle(stored as microseconds)

    %Find Channels

    ICurrentChan = strcmp(h_1.recChNames, ICurrentChanName);
    ICurrentID = find(ICurrentChan==1);

    VMChan = strcmp(h_1.recChNames, VMChanName);
    VMID = find(VMChan==1);

    WingBeatChan = strcmp(h_1.recChNames, WingBeatChanName);
    WingbeatID = find(WingBeatChan==1);

    LEDChan = strcmp(h_1.recChNames, LEDChanName);
    LEDID = find(LEDChan==1);

    LED2Chan = strcmp(h_1.recChNames, LED2ChanName);
    LED2ID = find(LED2Chan==1);

    LED = d1(:, LEDID);
    LED2 = d1(:, LED2ID);

    VelocXChan = strcmp(h_1.recChNames, VelocXChanName);
    VelocXID = find(VelocXChan==1);
    VelocYChan = strcmp(h_1.recChNames, VelocYChanName);
    VelocYID = find(VelocYChan==1);
    VelocZChan = strcmp(h_1.recChNames, VelocZChanName);
    VelocZID = find(VelocZChan==1);

    VM = d1(:, VMID);               %Load VM values from d1 array in VM variable
    Wingbeat = d1(:,WingbeatID);
    VelocX = d1(:,VelocXID);
    VelocY = d1(:,VelocYID);
    VelocZ = d1(:,VelocZID);

    x = (1:length(VM)) / sampling_rate;

%     %% Truncate data
%     raw_data_fig = figure;
%     plot(VM)
%     set(raw_data_fig, 'position', [1, 600, 1900, 450]);
% 
%     truncate_data = questdlg('Do you want to truncate data?','Truncate data?','Yes','No', 'No');
% 
%     if strcmpi (truncate_data, 'Yes')
%         cutofffig = figure;
%         set(cutofffig, 'position', [1, 600, 1900, 450]);
%         cutoff = ginput(1);
%         close (gcf);
% 
%         VM = VM(1:floor(cutoff(1)));
%         Wingbeat = Wingbeat(1:floor(cutoff(1)));
%         VelocX = VelocX(1:floor(cutoff(1)));
%         VelocY = VelocY(1:floor(cutoff(1)));
%         VelocZ = VelocZ(1:floor(cutoff(1)));
%         LED = LED(1:floor(cutoff(1)));
%     else
%         close(gcf);
%     end

    %% Smooth VM

    VM_smooth = smooth(VM,50, 'loess');

    %% Spike detection

    [spikes,spiketimes] = spikedetector(VM_smooth, sampling_rate);

    %% Plot mean spike shape
 if length(spiketimes)<nr_spikes_avg
        nr_spikes_avg = length(spiketimes);
        disp(['Number of avaraged spikes was reduced to ' num2str(nr_spikes_avg) ' because there are not enought spikes'])
    end
    %Remove first spike if it is at the ebegining of the rec
    %and would result in -intra_smooth values
    if spiketimes(1)-0.02*sampling_rate<VM_smooth(1)
        spiketimes(1) = [];
    end

    for i=2:nr_spikes_avg
        allspikes(:,i-1)=VM_smooth(spiketimes(i)-0.02*sampling_rate:spiketimes(i)+0.02*sampling_rate);
    end
    meanspike = mean(allspikes');
    figure
    plot(meanspike, 'r');


    %% Print out means

    mean_spikerate = mean(diff(spiketimes))/sampling_rate;
    spike_interval = diff(spiketimes)/sampling_rate;
    mean_spike_interval = mean(spike_interval);
    mean_spike_rate = 1/mean_spike_interval;
    mean_spike_amp = abs(meanspike(1)-max(meanspike));

    disp(['The mean spike interval is ' num2str(mean_spike_interval) 'Hz'])
    disp(['The mean spike rate is ' num2str(mean_spike_rate) 'Hz'])

    VM_medfilt = medfilt1(VM,1000);
    VM_median = median(VM_medfilt)

    %% Create analysis struct

    analysis(fly).VM = VM;
    analysis(fly).spikes = spikes;
    analysis(fly).medianVM = VM_median;
    analysis(fly).meanSpikeHz = mean_spike_rate;
    analysis(fly).meanSpikeIntervalHz = mean_spike_interval;
    analysis(fly).meanSpikeAmp = mean_spike_amp;
    analysis(fly).VelocX = VelocX;
    analysis(fly).VelocY = VelocY;
    analysis(fly).VelocZ = VelocZ;
    analysis(fly).file = file;
    analysis(fly).ncell = CellNumber(fly);
    analysis(fly).ntrial = TrialNumber(fly);
    analysis(fly).nfly = FlyNumber(fly);
    analysis(fly).ID = string(analysis(fly).file)+'_'+string(analysis(fly).nfly)...
        +'_'+string(analysis(fly).ncell)+'_'+string(analysis(fly).ntrial);
    analysis(fly).LED = LED;
    analysis(fly).Wingbeat = Wingbeat;

end

%%
%save('pre_processing','analysis','-v7.3');


