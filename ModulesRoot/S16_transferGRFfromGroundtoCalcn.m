function [Datastr] = S16_transferGRFfromGroundtoCalcn(Datastr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function transfert the GRF from the Ground coordinate to the Calcn
% coordinate, for better comparison between the wearable system and the
% golden standard optical system.
%
% By: Huawei Wang
% Date: 1/16/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the Opensim libraries
import org.opensim.modeling.*

% subject opensim model path
osMod = [Datastr.Info.SubjRoot, '\OS\', Datastr.Info.subjosmodfile];

% load the OpenSim model
OsimModel = Model(osMod);

% extract variable names
variableNames = string(Datastr.Resample.Sych.IKAngDataLabel(2:end));

% transfer the left and right GRF into the Calcn coordinate and saved it to
% the sychronized data layer

for irow =  1:length(Datastr.Resample.Sych.IKAngData(:, 1))

    % % extract variable data
    variableValues = Datastr.Resample.Sych.IKAngData(irow, 2:end);  % extract the state variable values

    % transfer the joint angles from degree to radian
    variableValues(1:3) = variableValues(1:3)*(pi/180);
    variableValues(7:end) = variableValues(7:end)*(pi/180);

    % get the transform matries
    [RotL, RotR] = getTransformMatrixBetweenGroundandCalcn(OsimModel, variableNames, variableValues);
    
    % prepare GRF vectors for the transformation
    lGRF = Datastr.Resample.Sych.ForcePlateGRFData(irow, 2:4);
    lCoP = Datastr.Resample.Sych.ForcePlateGRFData(irow, 5:7);
    
    rGRF = Datastr.Resample.Sych.ForcePlateGRFData(irow, 8:10);
    rCoP = Datastr.Resample.Sych.ForcePlateGRFData(irow, 11:13);
    
    if sum(lGRF) ~= 0
       lGRF_Calcn = RotL(1:3, 1:3)*lGRF';
       lCoP_Calcn = RotL*([lCoP, 1]');
       
       Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(irow, 2:4) = lGRF_Calcn;
       Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(irow, 5:7) = lCoP_Calcn;
    else
       Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(irow, 2:4) = [0, 0, 0];
       Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(irow, 5:7) = [0, 0, 0];
    end
    
    if sum(rGRF) ~= 0
       rGRF_Calcn = RotR(1:3, 1:3)*rGRF';
       rCoP_Calcn = RotR*([rCoP, 1]');
       
       Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(irow, 8:10) = rGRF_Calcn;
       Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(irow, 11:13) = rCoP_Calcn;
    else
       Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(irow, 8:10) = [0, 0, 0];
       Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(irow, 11:13) = [0, 0, 0];
    end
   
end

% remove the rest column of data if any
if length(Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(1, :)) > 13
    Datastr.Resample.Sych.ForcePlateGRFDataInCalcn(:, 14:end) = [];
end

end

