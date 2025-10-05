function analysis = binning_spikes(analysis, fieldName, binsize_factor, sampling_rate)
% BINVELOCITY Bins a specified velocity field (X, Y, or Z) in mm.
%
%   analysis = BINVELOCITY(analysis, fieldName, binsize_factor, sampling_rate)
%
%   Inputs:
%       analysis         - Struct array containing velocity fields
%       fieldName        - Name of field to bin (e.g. 'xveloc_in_mm')
%       binsize_factor   - Factor to multiply by sampling_rate to get bin size
%       sampling_rate    - Sampling rate of the data
%
%   Output:
%       analysis          - Updated struct with new binned field
%
%   Example:
%       analysis = binVelocity(analysis, 'xveloc_in_mm', 0.5, 1000);

    binsize = binsize_factor * sampling_rate;

    for k = 1:length(analysis)
        % Check that field exists
        if ~isfield(analysis(k), fieldName)
            warning('Field "%s" not found in analysis(%d). Skipping...', fieldName, k);
            continue;
        end

        data = analysis(k).(fieldName);
        binstarts = 1:binsize:length(data);
        binned = nan(1, length(binstarts)-1);


        for bin = 1:length(binstarts)-1
            binned(bin) = sum(data(binstarts(bin):binstarts(bin+1)), 'omitnan');
        end

        binned = binned *(1/binsize_factor);

        % Store result in a new field, e.g. xveloc_in_mm_binned
        newField = [fieldName '_binned'];
        analysis(k).(newField) = binned;
    end
end