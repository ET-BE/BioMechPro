function [Datastr] = S12_Sychronize(Datastr, MotionSychSign, ForceSychSign, EMGSychDelay)
% gBMPDynUI MotionSychSign=1; ForceSychSign=1; EMGSychDelay=1;
%
% Sychronize the portable and glod-standard (optical) system through the least
% square optimization.
% INPUT)
% - Datastr, the data structure with at least the fields:
%          .Resample.IMU
%          .Resample.IKID
%          .Resample.Force
%          .Resample.Insole
%
% - MotionSychSign, the sign of doing motion sychronization or not
% - ForceSychSign, the sign of doing force sychronization or not


%% remove fields if exist
if isfield(Datastr.Resample, 'Sych')
    if isfield(Datastr.Resample.Sych, 'EMG')
        Datastr.Resample.Sych = rmfield(Datastr.Resample.Sych, 'EMG');
    end
end

%% check if the data are suffcient for the sychronization
if strcmpi(Datastr.Info.Trial, 'static_pose')
    
else
    if strcmp(MotionSychSign, 'True')
       if isfield(Datastr.Resample.Sych, 'IKAngData') && isfield(Datastr.Resample, 'IMU')

           % select joints for the sychronization calculation [L/R hip, knee, ankle]
           IK_joints = Datastr.Resample.Sych.IKAngData(:, [15, 18, 19, 8, 11, 12]);
           IMU_joints = Datastr.Resample.IMU(:, [17, 20, 23, 7, 10, 13]);

           motion_SychTime = getSychronizationTime(IK_joints, IMU_joints, Datastr.Info.Trial)/Datastr.Resample.FrameRate;
           
           t_ik = 0:1/Datastr.Resample.FrameRate:(length(Datastr.Resample.Sych.IKAngData(:,19))-1)/Datastr.Resample.FrameRate;
           t_imu = 0:1/Datastr.Resample.FrameRate:(length(Datastr.Resample.IMU(:,23))-1)/Datastr.Resample.FrameRate;
           
           figure()
           plot(t_ik, Datastr.Resample.Sych.IKAngData(:, 19), '-', 'linewidth', 2)
           hold on
           plot(t_imu - motion_SychTime, Datastr.Resample.IMU(:, 23), '-', 'linewidth', 2)
           hold off
           legend('Optical', 'Wearable')
           title('motion sychronization')

       else
           error('Not enough motion data to apply schronization')
       end
    end

    if strcmp(ForceSychSign, 'True')
       if isfield(Datastr.Resample.Sych, 'ForcePlateGRFData') && isfield(Datastr.Resample, 'Insole')

           % select parameters for the sychronization calculation [L/R Fy]

           Fy_FP = [Datastr.Resample.Sych.ForcePlateGRFData(:, 3),...
                    Datastr.Resample.Sych.ForcePlateGRFData(:, 9)];
           Fy_IS_left = Datastr.Resample.Insole.Left(:, 23);
           Fy_IS_right = Datastr.Resample.Insole.Right(:, 23);

           minR = min(size(Fy_IS_left, 1), size(Fy_IS_right, 1));

           Fy_IS = [Fy_IS_left(1:minR), Fy_IS_right(1:minR)];

           force_SychTime = getSychronizationTime(Fy_FP, Fy_IS, Datastr.Info.Trial)/Datastr.Resample.FrameRate;     
           
           t_fp = 0:1/Datastr.Resample.FrameRate:(length(Datastr.Resample.Sych.ForcePlateGRFData(:,3))-1)/Datastr.Resample.FrameRate;
           t_is = 0:1/Datastr.Resample.FrameRate:(length(Datastr.Resample.Insole.Left(:, 23))-1)/Datastr.Resample.FrameRate;
           
            figure()
            plot(t_fp, Datastr.Resample.Sych.ForcePlateGRFData(:,3), '-', 'linewidth', 2)
            hold on
            plot(t_is - force_SychTime, Datastr.Resample.Insole.Left(:, 23), '-', 'linewidth', 2)
            hold off
            legend('Optical', 'Wearable')
            title('force sychronization')

       else
           error('Not enough force data to apply schronization')
       end
    end

    Datastr.Resample.Sych.DeltaT.Motion = motion_SychTime;
    Datastr.Resample.Sych.DeltaT.Force = force_SychTime;
    Datastr.Resample.Sych.DeltaT.EMG = EMGSychDelay;
    
    % save sych EMG data based on given EMGSychDelay value;
    nFramesIKAng = size(Datastr.Resample.Sych.IKAngData, 1);
    nFramesEMG = size(Datastr.Resample.EMG, 1);
    Datastr.Resample.Sych.EMG = interp1((0:nFramesEMG-1)/Datastr.Resample.FrameRate, Datastr.Resample.EMG(1:nFramesEMG, :),...
          (0:nFramesIKAng-1)/Datastr.Resample.FrameRate + EMGSychDelay, 'linear', 'extrap');

end

end