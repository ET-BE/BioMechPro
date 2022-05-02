function EMG_Envelope_MovingWindow = getMovingAverage(EMG_Envelope,movingWindow)
%UNTITLED this function calculate the emg envelope with a moving window of
%50 ms

[row, col] = size(EMG_Envelope);

if row < col % make sure that the EMG_Envelope's rows are frames
    EMG_Envelope = EMG_Envelope';
end

movingWindow = round(movingWindow);

EMG_Envelope_MovingWindow = zeros(max(row, col) - movingWindow, min(row, col));

for i = 1:max(row, col) - movingWindow + 1
    EMG_Envelope_MovingWindow(i, :) = mean(EMG_Envelope(i:i+movingWindow-1, :), 1);
end

if row < col  % rotate back the orignal data structure
    EMG_Envelope_MovingWindow = EMG_Envelope_MovingWindow';
end

end

