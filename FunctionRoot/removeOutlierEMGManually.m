function hsMatrixUpdate = removeOutlierEMGManually(hsMatrix, EMGdata, EMG_id)
%UNTITLED manually remove the outlier based on the EMG recordings.

    goodFlag = 0;

    while ~goodFlag

        numStrike = size(hsMatrix, 1);

        fig1 = figure();
        rowNum = ceil(size(EMGdata, 2)/2);
        for e = 1:size(EMGdata, 2)
            subplot(rowNum, 2, e)
            for i = 1:numStrike
                plot((e-1)*1000 + 1: (e-1)*1000 + 1 + hsMatrix(i, 2) - hsMatrix(i, 1),...
                    EMGdata(hsMatrix(i, 1):hsMatrix(i, 2), e), 'linewidth', 0.5)
                hold on
            end

            if e == 1
                plot(0, mean(EMGdata(:, 1)), 'go')
            end
            hold off

            title(EMG_id(e))
            if e == 9 || e == 10
                xlabel('Time (s)')
            end
            if e == 2
                legend('EMG')
            end
        end

        [x,y]=ginput(1);

        if sqrt(x^2 + (y-mean(EMGdata(:, 1)))^2) < 0.5
            goodFlag = 1;
        else
            iemg = ceil(x/1000);
            FyError = 10;
            for i = 1:size(hsMatrix, 1)
                 FyInterp = EMGdata(hsMatrix(i, 1) + floor(x)-(iemg-1)*1000 - 1, iemg)*(1-x+floor(x)) +...
                            EMGdata(hsMatrix(i, 1) + ceil(x)-(iemg-1)*1000 - 1, iemg)*(1+x-ceil(x));

                if abs(FyInterp - y) <= FyError
                    FyError = abs(FyInterp - y);
                    Error_trial = i;
                end
           end
           hsMatrix(Error_trial, :) = [];
           numStrike = numStrike - 1;
       end
       close(fig1)

    end
    
    hsMatrixUpdate = hsMatrix;
    
end

