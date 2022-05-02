function [Datastr] = S19_removeEMGBottomOffset(Datastr, leftLegEMG, rightLegEMG, examEMGID)
% gBMPDynUI leftLegEMG = 1;  rightLegEMG = 1; examEMGID = 1;
%
%
% INPUT)
% Datastr: structure, with at least the fields
%     .Resample.Sych.EMG
%
% NOTES)
% Function based on MOtion data elaboration TOolbox for
% NeuroMusculoSkeletal applications (MOtoNMS).
% Copyright (C) 2012-2014 Alice Mantoan, Monica Reggiani

% TODO

%% Check input

    EMG_id = Datastr.EMG.DataLabel;

    EMGdata = Datastr.Resample.Sych.EMG;

    if strcmp(rightLegEMG, 'True')  % if right leg EMG were recorded
        hsMatrix_r = Datastr.Resample.Sych.Average.hsMatrix_right;

        emg_index_r = [];
        for emg_i = 1:length(examEMGID)
            if examEMGID(emg_i) == 1 && contains(EMG_id(emg_i), '_r')
                emg_index_r = [emg_index_r, emg_i];
            end
        end

        zeroOffset_r = getZeroOffset(hsMatrix_r,...
                                  EMGdata(:, emg_index_r), EMG_id(emg_index_r));

        Datastr.Resample.Sych.EMG(:, emg_index_r) = EMGdata(:, emg_index_r) - zeroOffset_r;
        
    end
    
    if strcmp(leftLegEMG, 'True')
        hsMatrix_l = Datastr.Resample.Sych.Average.hsMatrix_left;

        emg_index_l = [];
        for emg_i = 1:length(examEMGID)
            if examEMGID(emg_i) == 1 && contains(EMG_id(emg_i), '_l')
                emg_index_l = [emg_index_l, emg_i];
            end
        end

        zeroOffset_l = getZeroOffset(hsMatrix_l,...
                                  EMGdata(:, emg_index_l), EMG_id(emg_index_l));

        Datastr.Resample.Sych.EMG(:, emg_index_l) = EMGdata(:, emg_index_l) - zeroOffset_l;
        
    end

end