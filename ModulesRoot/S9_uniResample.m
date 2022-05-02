function Datastr = S9_uniResample(Datastr, resampleRate, removeFrames)
% gBMPDynUI resampleRate=1; removeFrames = 1;
% 
% INPUT)
% - Datastr: structure, containing at least the fields:
% .Marker.FrameRate
% .Marker.MarkerData
% 
% .Force.FrameRate
% .Force.RightForceData
% .Force.LeftForceData
%
% .EMG.FrameRate
% .EMG.EMGLinEnv
%
% .Insole.FrameRate
% .Insole.FilteredLeft
% .Insole.FilteredRight
%
% .IMU.IMUFrameRate
% .IMU.IMUData
% 
% - resampleRate: the resampling rate for all signals
% - beginEndRemove: whether to remove 1 second at the beginning and end for
% inverse kinematics and inverse dynamics. This is to avoid the unexpected
% peaks.
% 
% OUTPUT)
% 
% NOTES)
%

% remove the original resampling data in case the data lengths are different
if isfield(Datastr, 'Resample')
    
    if isfield(Datastr.Resample, 'Marker')
        Datastr.Resample = rmfield(Datastr.Resample, 'Marker');
    end
    if isfield(Datastr.Resample, 'Force')
        Datastr.Resample = rmfield(Datastr.Resample, 'Force');
    end
    if isfield(Datastr.Resample, 'EMG')
        Datastr.Resample = rmfield(Datastr.Resample, 'EMG');
    end
    if isfield(Datastr.Resample, 'Insole')
        Datastr.Resample = rmfield(Datastr.Resample, 'Insole');
    end
    if isfield(Datastr.Resample, 'IMU')
        Datastr.Resample = rmfield(Datastr.Resample, 'IMU');
    end
    
end

%% check input signals

if isempty(resampleRate)
   resampleRate = min([Datastr.Marker.FrameRate, Datastr.Force.FrameRate,...
       Datastr.EMG.FrameRate, Datastr.Insole.FrameRate, Datastr.IMU.IMUFrameRate]); 
end

Datastr.Resample.FrameRate = resampleRate;

if removeFrames
    remove_samples = removeFrames;
else
    remove_samples = 0;
end



%% start resampling
for mk = 1:length(Datastr.Marker.MarkerData(:, 1, 1))
    
    MrData = squeeze(Datastr.Marker.MarkerData(mk, :, :));
    time_frame = linspace(0,...
        (length(MrData(1, :))-1)/Datastr.Marker.FrameRate,...
        length(MrData(1, :)));

    Datastr.Resample.Marker(:, :, mk) = ...
        resample(MrData', time_frame, resampleRate, 'spline');
    
    % ressign the NaNs to the resampled marker trajectories
    if ~isempty(find(isnan(MrData(1, :))==1))
        id_nan_ori = find(isnan(MrData(1,:))==1);
        id_nan_new = ceil(1 + resampleRate*(id_nan_ori-1)/Datastr.Marker.FrameRate);
        id_nan_new(id_nan_new > size(Datastr.Resample.Marker, 1)) = [];
        Datastr.Resample.Marker(id_nan_new, :, mk) = NaN;
    end
end

% remove 10 frames of data from the beginning and end
Datastr.Resample.Marker(1:remove_samples, :, :) = [];
Datastr.Resample.Marker(end-remove_samples:end, :, :) = [];

if isfield(Datastr, 'Force')

    Datastr.Resample.Force.Right = ...
        resample(Datastr.Force.RightForceData, linspace(0,...
        (length(Datastr.Force.RightForceData(:, 1))-1)/Datastr.Force.FrameRate,...
        length(Datastr.Force.RightForceData(:, 1))), resampleRate, 'spline');

    % remove 10 frames of data from the beginning and end
    Datastr.Resample.Force.Right(1:remove_samples, :) = [];
    Datastr.Resample.Force.Right(end-remove_samples:end, :) = [];

    Datastr.Resample.Force.Left = ...
        resample(Datastr.Force.LeftForceData, linspace(0,...
        (length(Datastr.Force.LeftForceData(:, 1))-1)/Datastr.Force.FrameRate,...
        length(Datastr.Force.LeftForceData(:, 1))), resampleRate, 'spline');

    % remove 10 frames of data from the beginning and end
    Datastr.Resample.Force.Left(1:remove_samples, :) = [];
    Datastr.Resample.Force.Left(end-remove_samples:end, :) = [];
    
end

if isfield(Datastr, 'EMG')

    if isfield(Datastr.EMG, 'EMGLinEnv')
        Datastr.Resample.EMG = ...
        resample(Datastr.EMG.EMGLinEnv, linspace(0,...
        (length(Datastr.EMG.EMGLinEnv(:, 1))-1)/Datastr.EMG.FrameRate,...
        length(Datastr.EMG.EMGLinEnv(:, 1))), resampleRate, 'spline');
    elseif isfield(Datastr.EMG, 'EMGLinEnvNor')
        Datastr.Resample.EMG = ...
        resample(Datastr.EMG.EMGLinEnvNor, linspace(0,...
        (length(Datastr.EMG.EMGLinEnvNor(:, 1))-1)/Datastr.EMG.FrameRate,...
        length(Datastr.EMG.EMGLinEnvNor(:, 1))), resampleRate, 'spline');
    end

    % remove 10 frames of data from the beginning and end
    Datastr.Resample.EMG(1:remove_samples, :) = [];
    Datastr.Resample.EMG(end-remove_samples:end, :) = [];
end

if isfield(Datastr, 'Insole')
    Datastr.Resample.Insole.Right = ...
        resample(Datastr.Insole.FilteredRight, linspace(0,...
        (length(Datastr.Insole.FilteredRight(:, 1))-1)/Datastr.Insole.FrameRate,...
        length(Datastr.Insole.FilteredRight(:, 1))), resampleRate, 'spline');

    % remove 10 frames of data from the beginning and end
    Datastr.Resample.Insole.Right(1:remove_samples, :) = [];
    Datastr.Resample.Insole.Right(end-remove_samples:end, :) = [];

    Datastr.Resample.Insole.Left = ...
        resample(Datastr.Insole.FilteredLeft, linspace(0,...
        (length(Datastr.Insole.FilteredLeft(:, 1))-1)/Datastr.Insole.FrameRate,...
        length(Datastr.Insole.FilteredLeft(:, 1))), resampleRate, 'spline');

    % remove 10 frames of data from the beginning and end
    Datastr.Resample.Insole.Left(1:remove_samples, :) = [];
    Datastr.Resample.Insole.Left(end-remove_samples:end, :) = [];
end

if isfield(Datastr, 'IMU')
    Datastr.Resample.IMU = ...
        resample(Datastr.IMU.IMUData, linspace(0,...
        (length(Datastr.IMU.IMUData(:, 1))-1)/Datastr.IMU.IMUFrameRate,...
        length(Datastr.IMU.IMUData(:, 1))), resampleRate, 'spline');

    % remove 10 frames of data from the beginning and end
    Datastr.Resample.IMU(1:remove_samples, :) = [];
    Datastr.Resample.IMU(end-remove_samples:end, :) = [];

    Datastr.Resample.CoM = ...
        resample(Datastr.IMU.CoM, linspace(0,...
        (length(Datastr.IMU.CoM(:, 1))-1)/Datastr.IMU.IMUFrameRate,...
        length(Datastr.IMU.CoM(:, 1))), resampleRate, 'spline');

    % remove 10 frames of data from the beginning and end
    Datastr.Resample.CoM(1:remove_samples, :) = [];
    Datastr.Resample.CoM(end-remove_samples:end, :) = [];

end

end