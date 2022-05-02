function [Datastr] = S20_extractMVCAverageGaits(Datastr, MVCtrial,...
                dynNormalization, leftLegEMG, rightLegEMG, timeMovingWindow)
% gBMPDynUI MVCtrial=1;  dynNormalization = 1; leftLegEMG=1; rightLegEMG=1; timeMovingWindow = 1;
%
%
% INPUT)
% Datastr: structure, with at least the fields
%     .Resample.Sych.EMG
% MVCtrial: trial name of the MVC recording
% dynNormalization: whether apply normization using the MVC or the higest
% emg peaks among all trials
%
% NOTES)
% Function based on MOtion data elaboration TOolbox for
% NeuroMusculoSkeletal applications (MOtoNMS).
% Copyright (C) 2012-2014 Alice Mantoan, Monica Reggiani

% TODO
%This module could be more efficient, loading the trial to be normalized in
%each trial processing from BioMechPro. Nevbertheless, to keep MOtoNMS
%functions the same, all the processing is done when we loading the first
%trial in BioMechpro. All trials are processed in the first loading. The next loadings are skipped 
%(as they were already proccessed during the first trial. Read the code!)


%% Check input

EMG_id = Datastr.EMG.DataLabel;

if strcmp(dynNormalization, 'True')
    
    fprintf('extract dynamic MVC from all data trials...\n')
            
    subj_folder = Datastr.Info.SubjRoot;
        
    EMGdata = Datastr.Resample.Sych.EMG;
    
    maxEMGAve = zeros(1, size(EMGdata, 2));
    
    % calculate the frame number of the moving window
    frameMovingWindow_Gait = Datastr.Resample.FrameRate*0.1;
    
    if strcmp(rightLegEMG, 'True')  % if right leg EMG were recorded
        hsMatrix_r = Datastr.Resample.Sych.Average.hsMatrix_right;

        emg_index_r = [];
        for emg_i = 1:length(EMG_id)
            if contains(EMG_id(emg_i), '_r')
                emg_index_r = [emg_index_r, emg_i];
            end
        end
               
        % apply moving average of 50 ms windows
        % movingAveEMG_r = getMovingAverage(EMGdata(:, emg_index_r), frameMovingWindow_Gait);
        maxEMGAve(emg_index_r) = getMaximumEMGInAveragedGait(hsMatrix_r,...
                                 EMGdata(:, emg_index_r));
    end
    
    if strcmp(leftLegEMG, 'True')
        hsMatrix_l = Datastr.Resample.Sych.Average.hsMatrix_left;

        emg_index_l = [];
        for emg_i = 1:length(EMG_id)
            if contains(EMG_id(emg_i), '_l')
                emg_index_l = [emg_index_l, emg_i];
            end
        end

        % movingAveEMG_l = getMovingAverage(EMGdata(:, emg_index_l), frameMovingWindow_Gait);
        maxEMGAve(emg_index_l) = getMaximumEMGInAveragedGait(hsMatrix_l,...
                                  EMGdata(:, emg_index_l));
    end
  
    if ~isempty(MVCtrial)
        if isfile([subj_folder, '\MVCvalue.mat'])
            MVCvalue = importdata([subj_folder, '\MVCvalue.mat']);
        else
            extractMVCtrial([Datastr.Info.SubjRoot, MVCtrial], timeMovingWindow, EMG_id)
            MVCvalue = importdata([subj_folder, '\MVCvalue.mat']);
        end
    end
    
    if isfile([subj_folder, '\dynMVCvalue.mat'])
       dynMVCvalue = importdata(strcat(subj_folder, '\dynMVCvalue.mat'));
       MVCvalueMid = max(MVCvalue, dynMVCvalue);
       MVCvalueMax = max(MVCvalueMid, maxEMGAve);
       
       save(strcat(subj_folder, '\dynMVCvalue.mat'),  'MVCvalueMax')
    else
        MVCvalueMax = max(maxEMGAve, MVCvalue);
        save(strcat(subj_folder, '\dynMVCvalue.mat'),  'MVCvalueMax')
    end

elseif strcmp(dynNormalization, 'False')
    
    if isempty(MVCtrial)
        fprintf('No MVC trial provided, cannot extract...\n')
        
    else
        extractMVCtrial([Datastr.Info.SubjRoot, MVCtrial], timeMovingWindow, EMG_id);
    end
end

    function MVCvalue = extractMVCtrial(mvcFileName, timeMovingWindow, EMG_id)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract MVCvalue from MVC trial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the MVC data
        MVCdata = importdata(mvcFileName);

        % Get sample frequency
        fs_emg = MVCdata.Analog.Frequency;

        % Get EMG
        emgData = MVCdata.Analog.Data(17:25, :);

        % EMG envelope
        EMG_envelope = getEMGEnvelop(emgData, fs_emg);
        
        % calculate the frame number of the moving window
        frameMovingWindow = fs_emg*timeMovingWindow;
        
        % apply moving average of 50 ms windows
        EMG_envelope_movingaverage = getMovingAverage(EMG_envelope, frameMovingWindow);
        
        fig2 = figure();
        for ee = 1:9
            subplot(5, 2, ee)
            plot(EMG_envelope_movingaverage(:, ee), 'linewidth', 0.5)
            
            title(EMG_id(ee))
            ylabel('% MVC')
            if ee == 9 || ee == 10
                xlabel('Time (s)')
            end
            if ee == 2
                legend('EMG')
            end
        end
        
        [~, yMVC] = ginput(length(EMG_id));

        MVCvalue = min(yMVC', max(EMG_envelope_movingaverage));
        
        slashIndexes = strfind(mvcFileName, '\');
        
        save([mvcFileName(1:slashIndexes(end-1)) 'MVCvalue.mat'], 'MVCvalue')
        
        close(fig2)

    end

end