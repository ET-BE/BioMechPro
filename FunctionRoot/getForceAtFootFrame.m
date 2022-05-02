
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test
% relocate the force at the ground to the foot frame
% By: Huawei Wang
% Date: Nov. 10, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

% load the Opensim libraries
import org.opensim.modeling.*

% load the OpenSim model
osimModelFile = 'C:\Users\WangH\Documents\PortableSystem_DataAnalysis\Processed_data\Subj04\OS\UTmodel\gait2392_simbody_subj04_scaled_1.osim';
OsimModel = Model(osimModelFile);

% get frame list
frameList = OsimModel.getFrameList();

% select ground frame
framGround = getOsimFrame(OsimModel, 'ground');

% select left calcn frame
framLCalcn = getOsimFrame(OsimModel, 'calcn_l');

% select right calcn frame
framRCalcn = getOsimFrame(OsimModel, 'calcn_r');

% initilize the simulation model
% OsimModel.buildSystem();

% load the inverse kinematics results
ikFile = 'C:\Users\WangH\Documents\PortableSystem_DataAnalysis\Processed_data\Subj04\OS\DataFiles\Subj04walk_36IK.mot';
ikData = importdata(ikFile);

% extract variable names
variableNames = string(ikData.colheaders(2:end));

% % extract variable data
variableValues = ikData.data(11, 2:end);  % extract the state variable values

% transfer the joint angles from degree to radian
variableValues(1:3) = variableValues(1:3)*(pi/180);
variableValues(7:end) = variableValues(7:end)*(pi/180);


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

RotL = [lR, lT']
RotR = [rR, rT']












