function [spikes,spiketimes] = spikedetector(VM_smooth, sampling_rate)

 checkbox_diff=0;                % Recommended value is '0'. If set to 0, threshold is determined on RAW data. If set to 1, threshold is determined on 'diff' voltage trace
    if checkbox_diff==0
        refractory=0.005;
    else
        refractory=0.005;
    end
    
    
    threshfig_raw=figure('name','RAW');
    set(threshfig_raw, 'position', [1, 600, 1900, 450]);
    plot(VM_smooth,'k');
    threshfig_diff=figure('name','DIFF');
    set(threshfig_diff, 'position', [0, 0, 1900, 450])
    plot(diff(VM_smooth),'k');
    
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
    
    
    threshfig=figure;
    set(threshfig, 'position', [1, 1, 1900, 1000])
    
    if checkbox_diff==1
        plot(diff(VM_smooth),'k');
    else
        plot(VM_smooth,'k');
    end
    set(threshfig, 'position', [1, 1, 2000, 1000])
    xlim([5*sampling_rate,150*sampling_rate])
    if checkbox_diff==1
        ylim([0, 2])%max(diff(VM_smooth))])y
    end
    
    thresh=ginput(1);
    thresh=thresh(1,2);
    
    close(threshfig);
    
    overthresh=zeros(size(VM_smooth));
    if checkbox_diff ==1
        overthresh(diff(VM_smooth)>thresh)=1;
        tricky=diff(overthresh);
    else
        overthresh(VM_smooth>thresh)=1;
        tricky=(overthresh);
    end
    
    
    spikes=nan(size(VM_smooth));
    spikes(tricky==1)=1;
    
    
    for s=3:length(spikes)
        if spikes(s-1)==1;
            spikes(s:s+refractory*sampling_rate)=nan;
        else
        end
    end
    
    %                 if checkbox_safespike ==1 %% Find and delete false positive spikes
    %                     for j=3:length(spikes)-2 %Does not take the first two and last two "spikes" into account
    %                         if spikes(j)==1
    %                             tmp_prespike=VM_smooth(j-1);
    %                             tmp_postspike=VM_smooth(j+1);
    %
    %                             if tmp_prespike>VM_smooth(j) && tmp_postspike<VM_smooth(j)
    %                                 spikes(j)=nan;
    %                                 idx_spike_missmatch = idx_spike_missmatch+1;
    %                             end
    %
    %                         end
    %                     end
    %                 end
    %
    %                 if isempty(idx_spike_missmatch)&& checkbox_safespike ==1
    %                     disp('No Spikes missmatches detected')
    %                 elseif checkbox_safespike ==1
    %                     disp([num2str(idx_spike_missmatch) ' Spikes missmatches detected'])
    %                 end
    
    spiketimes=find(spikes==1);