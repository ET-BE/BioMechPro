function [RotL, RotR] = getTransformMatrixBetweenGroundandCalcn(OsimModel, variableNames, variableValues)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the transform matrix between ground and calcn based on the opensim
% model and the inverse kinematics results.
% By: Huawei Wang
% Date: Nov. 10, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% select ground frame
framGround = getOsimFrame(OsimModel, 'ground');

% select left calcn frame
framLCalcn = getOsimFrame(OsimModel, 'calcn_l');

% select right calcn frame
framRCalcn = getOsimFrame(OsimModel, 'calcn_r');

% get the Opensim model states
OsimModelSim = initSystem(OsimModel);

% update coordinates
CoordinateSet = OsimModel.getCoordinateSet();
numCoordinate = OsimModel.getNumCoordinates();

for i = 1:numCoordinate
    
    coordinate = CoordinateSet.get(i-1);
    
    coorId = find(strcmp(variableNames, coordinate.getName()));
    
    coordinate.setValue(OsimModelSim, variableValues(coorId))
    
end

% extract the transform matrix between ground and the Calcn
transLCalcn = framGround.findTransformBetween(OsimModelSim, framLCalcn);
transRCalcn = framGround.findTransformBetween(OsimModelSim, framRCalcn);

% transfer the Opensim transform matrix to Matlab matrix
lR = osimMat33toMatrx(transLCalcn.R().asMat33());
lT = osimVec3ToArray(transLCalcn.T());

rR = osimMat33toMatrx(transRCalcn.R().asMat33());
rT = osimVec3ToArray(transRCalcn.T());

RotL = [lR, lT'];
RotR = [rR, rT'];

end











