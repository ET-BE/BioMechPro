function zeroOffset = getMaximumEMGInAveragedGait(hsMatrix, EMGdata)
% Get zero offset of the averaged EMG signals

    numStrike = size(hsMatrix, 1);
    zeroOffset = zeros(1, size(EMGdata, 2));
    
    for e = 1:size(EMGdata, 2)  
        for i = 1:numStrike-2
            zeroOffset(e) = max(zeroOffset(e), max(EMGdata(hsMatrix(i, 1):hsMatrix(i, 2), e)));
        end
    end
end