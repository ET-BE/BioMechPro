function [Datastr] = S18_removeOutliersFromEMG(Datastr, leftLegEMG, rightLegEMG)
% gBMPDynUI leftLegEMG = 1; rightLegEMG = 1;
% 
% INPUT)
% Datastr: structure, with at least the fields
%     .Resmaple.Sych.EMG
% 
% NOTES)
% Function based on MOtion data elaboration TOolbox for
% NeuroMusculoSkeletal applications (MOtoNMS).
% Copyright (C) 2012-2014 Alice Mantoan, Monica Reggiani

% TODO
% 

%% Check input

    EMGdata = Datastr.Resample.Sych.EMG;

    if strcmp(rightLegEMG, 'True')  % if right leg EMG were recorded
        hsMatrix_r = Datastr.Resample.Sych.Average.hsMatrix_right;

        emg_index_r = [];
        for emg_i = 1:length(EMG_id)
            if contains(EMG_id(emg_i), '_r')
                emg_index_r = [emg_index_r, emg_i];
            end
        end

        hsMatrixUpdate_r = removeOutlierEMGManually(hsMatrix_r,...
                                  EMGdata(:, emg_index_r), EMG_id(emg_index_r));

        Datastr.Resample.Sych.Average = rmfield(Datastr.Resample.Sych.Average, 'hsMatrix_right');
        Datastr.Resample.Sych.Average.hsMatrix_right = hsMatrixUpdate_r;
    end
    
    if strcmp(leftLegEMG, 'True')
        hsMatrix_l = Datastr.Resample.Sych.Average.hsMatrix_left;

        emg_index_l = [];
        for emg_i = 1:length(EMG_id)
            if contains(EMG_id(emg_i), '_l')
                emg_index_l = [emg_index_l, emg_i];
            end
        end

        hsMatrixUpdate_l = removeOutlierEMGManually(hsMatrix_l,...
                                  EMGdata(:, emg_index_l), EMG_id(emg_index_l));

        Datastr.Resample.Sych.Average = rmfield(Datastr.Resample.Sych.Average, 'hsMatrix_left');
        Datastr.Resample.Sych.Average.hsMatrix_right = hsMatrixUpdate_l;
    end
end