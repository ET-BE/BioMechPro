function [lCalcnOriginInGround, rCalcnOriginInGround] = getCalcnLocationInGroundFrame(osimModelFile, ikData, ikDataLabel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relocate the force at the ground to the foot frame
% By: Huawei Wang
% Date: Nov. 10, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import org.opensim.modeling.*

% load the OpenSim model
OsimModel = Model(osimModelFile);

% select left calcn frame
framLCalcn = getOsimFrame(OsimModel, 'calcn_l');

% select right calcn frame
framRCalcn = getOsimFrame(OsimModel, 'calcn_r');

% extract variable names
variableNames = ikDataLabel;

variableValues = zeros(size(ikData));
% transfer the joint angles from degree to radian
variableValues(:, 1:3) = ikData(:, 1:3)*(pi/180);
variableValues(:, 7:end) = ikData(:, 7:end)*(pi/180);

lCalcnOriginInGround = zeros(size(variableValues, 1), 3);
rCalcnOriginInGround = zeros(size(variableValues, 1), 3);


% initilize the simulation model
OsimModelSim = initSystem(OsimModel);

% get the Opensim model states
% OsimState = OsimModelSim.State;

for f = 1:size(variableValues, 1)
        
    % update the state variable (posture: joint angles)
    [OsimModel, OsimModelSim] = updOsimStateVariableValue(OsimModel, OsimModelSim, variableNames, variableValues(f, :));

    lCalcnOriginInGround(f, :) = osimVec3ToArray(framLCalcn.getPositionInGround(OsimModelSim));
    rCalcnOriginInGround(f, :) = osimVec3ToArray(framRCalcn.getPositionInGround(OsimModelSim));
end

end








