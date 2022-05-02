function [lCalcnOriginInGround,rCalcnOriginInGround] = getCalcnLocationInGroundFramebyFrame(osimModelFile, variableNames, variableValues)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% load the Opensim libraries
        import org.opensim.modeling.*

        % load the OpenSim model
        OsimModel = Model(osimModelFile);

        % select left calcn frame
        framLCalcn = getOsimFrame(OsimModel, 'calcn_l');

        % select right calcn frame
        framRCalcn = getOsimFrame(OsimModel, 'calcn_r');
        
        % initilize the simulation model
        OsimModelSim = initSystem(OsimModel);

        % get the Opensim model states
        OsimState = OsimModelSim.State;
        
        % update the state variable (posture: joint angles)
        [OsimModel, OsimState] = updOsimStateVariableValue(OsimModel, OsimState, variableNames, variableValues);

        lCalcnOriginInGround = osimVec3ToArray(framLCalcn.getPositionInGround(OsimModelSim));
        rCalcnOriginInGround = osimVec3ToArray(framRCalcn.getPositionInGround(OsimModelSim));
end

