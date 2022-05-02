function EMG_envelope_subs = getEMGEnvelop(emgData, fs_emg)
%Get the envelope of the EMG analog channels

    emgData(:, isnan(emgData(1, :))) = [];   % remove nan data, normally at the end

    % Get EMG Linear Envelope
    %% first high pass filter 25 Hz to remove wire artifacts
    [b, a] = butter(2, 25/(fs_emg/2), 'high'); 
    EMG_filter = filtfilt(b, a, emgData');

    %% then rectify the EMG signals to get all positive values
    EMG_rectify = abs(EMG_filter); 

    %% finally use the 6 Hz low pass filter to get envelop
    [b, a] = butter(2, 6/(fs_emg/2), 'low'); 
    EMG_envelope = filtfilt(b, a, EMG_rectify);
    
%     numEMG = size(EMG_envelope, 2);
%     y = zeros(1, numEMG);
%     
%     for e = 1:numEMG
%         fig1 = figure();
%         plot(EMG_envelope(:, e), 'linewidth', 0.5)
%         ylim([0, 0.1])
%         [~, y(e)] = ginput(1);
%         close(fig1)
%     end

    EMG_envelope_subs = EMG_envelope - min(EMG_envelope);
    
    % EMG_envelope is column data now
end

