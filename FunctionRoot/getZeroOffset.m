function zeroOffset = getZeroOffset(hsMatrix, EMGdata, EMG_id)
% Get zero offset of the averaged EMG signals

    numStrike = size(hsMatrix, 1);
    rowNum = ceil(size(EMGdata, 2)/2);
    
    zeroOffset = [];

    fig1 = figure();
    for e = 1:size(EMGdata, 2)
        subplot(rowNum, 2, e)
        for i = 1:numStrike
            plot((e-1)*1000 + 1: (e-1)*1000 + 1 + hsMatrix(i, 2) - hsMatrix(i, 1),...
                EMGdata(hsMatrix(i, 1):hsMatrix(i, 2), e), 'linewidth', 0.5)
            hold on
        end

        plot((e-1)*1000, mean(EMGdata(:, e)), 'go')
        hold on

        title(EMG_id(e))

        ylabel('% MVC')
        if e == 9 || e == 10
            xlabel('Time (s)')
        end
        if e == 2
            legend('EMG')
        end
    

        [x,y]=ginput(1);

        if sqrt((x - (e-1)*1000)^2 + (y-mean(EMGdata(:, e)))^2) < 0.5
            zeroOffset = [zeroOffset, 0];
        else
            zeroOffset = [zeroOffset, y];
        end
    
    end
close(fig1) 
    
end