function [Datastr] = S17_detHS(Datastr)
% gBMPDynUI
%
% plot the joint motion data of two systems, after the sychorinization
% INPUT)
% - Datastr, the data structure with at least the fields:
%          .Resample.Sych.IKAngData
%          .Resample.Sych.IDTrqData
%          .Resample.Sych.IMUAngData
%          .Resample.Sych.IDAngData_Portable
%          .Resample.Sych.ForcePlateGRFData
%          .Resample.Sych.InsoleGRFData
%
% - trial_name, the trial names that will excute this average function
%% check variables
if ~ isfield(Datastr.Resample.Sych, 'ForcePlateGRFData')
    warning('No force plate force data for averaging')
end
    
%% average the gait cycle
if strcmpi(Datastr.Info.Trial(1:4), 'walk') || strcmpi(Datastr.Info.Trial(1:3), 'run')
   % in walking and running trials, using Fy to detect gait cycles
   hsMatrix_l = getHS(Datastr.Resample.Sych.ForcePlateGRFData(:, 3), Datastr.Info.Trial(1:4));
   hsMatrix_r = getHS(Datastr.Resample.Sych.ForcePlateGRFData(:, 9), Datastr.Info.Trial(1:4));

   Datastr.Resample.Sych.Average.hsMatrix_left = hsMatrix_l;
   Datastr.Resample.Sych.Average.hsMatrix_right = hsMatrix_r;

else
    fig1 = figure();
        plot(Datastr.Resample.Sych.ForcePlateGRFData(:, 9))
        hold on
        plot(Datastr.Resample.Sych.ForcePlateGRFData(:, 3))
        hold on
        legend('Right Fy', 'Left Fy')
        ylabel('N')
        xlabel('Frames')
        
        hsMatrix_r = zeros(8, 2);
        
    % select 8 movement cycles for jump, land, lunge, and squat 
    for points = 1:16
        [x, y] = ginput(1);
        if rem(points, 2)
            plot(x, y, 'ro', 'linewidth', 2)
            hold on
            hsMatrix_r(ceil(points/2), 1) = round(x);
        else
            plot(x, y, 'bo', 'linewidth', 2)
            hold on
            hsMatrix_r(ceil(points/2), 2) = round(x);
        end
    end
    close(fig1)
Datastr.Resample.Sych.Average.hsMatrix_right = hsMatrix_r;
end

   
end



