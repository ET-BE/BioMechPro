function [Datastr] = S8_procEMG(Datastr)
% gBMPDynUI
%EMGs_Input={'TAL','GML','RFL'};
% Detrend, rectify, and low pass filter EMG
%
% INPUT)
% Datastr: structure, with at least the fields
%     .EMG.channels
%     .EMG.EMGFrameRate
%
% mvcNormalization: whether use mvc to normalize EMG
% dynNormalization: whether use dynamic mvc to normalize EMG
%
% OUTPUT)
% Datastr: structure with no new added fields
%     .EMG.envelope
%

%% Check input

if isfield(Datastr,'EMG')
    fprintf(['Trial ' Datastr.Info.Trial ' processed (flagEMG=1).\n ']) % All trials were processsed during the first trial        
else
    warning(['No EMG data in file ' Datastr.Info.Trial '. Skipping.'])
end
    
    
    %% Do Stuff
    %Input BioMechPro
    
    % %% ------------------------------------------------------------------------
    % %                      START/STOP COMPUTATION
    % %--------------------------------------------------------------------------

    if isfield(Datastr,'EMG') %only .mat files with EMG (and skip flagEMG file)
        
        emgData = Datastr.EMG.Channels;
        fs_emg = Datastr.EMG.FrameRate;
        % get EMG envelope
        EMG_envelope = getEMGEnvelop(emgData', fs_emg);

        Datastr.EMG.EMGLinEnv = EMG_envelope;
        Datastr.EMG.mvcNorFlag = 0;
        Datastr.EMG.dynNorFlag = 0;

    else
            warning(['No EMG data in file ' Datastr.Info.Trial '. Skipping.']);
    end

end