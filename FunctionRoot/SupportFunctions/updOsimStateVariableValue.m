
function [OsimModel, OsimModelSim] = updOsimStateVariableValue(OsimModel, OsimModelSim, variableNames, variableValues)

    import org.opensim.modeling.*
    % initilize the simulation model
%     OsimModelSim = initSystem(OsimModel);

    % get the Opensim model states
%     OsimState = OsimModelSim.State();

    % get the Opensim model state names
    OsimStateNames = OsimModel.getStateVariableNames();

    for i = 0:OsimStateNames.getSize()-1  % run over all States
       stateName = string(OsimStateNames.get(i));  % get current state name
       bshIndex = strfind(stateName, '/');  % get the state type

       if extractBetween(stateName, bshIndex(end)+1, strlength(stateName)) == "value"  % if the state type is 'value'

           % find the location of the state name inside the variable list
           variableIndex = find(variableNames == extractBetween(stateName, bshIndex(end-1)+1, bshIndex(end)-1));

           if  ~isempty(variableIndex) % if the location is not empty
               % set the given value of this state
               OsimModel.setStateVariableValue(OsimModelSim.State, OsimStateNames.get(i), variableValues(variableIndex));
           end

       end

    end

end
