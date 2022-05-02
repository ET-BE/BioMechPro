function [Datastr] = M4_getHumanRJCI(Datastr)
% See getRJCI for info

%% Get info

subjsex = Datastr.Info.subjsex;
subjmass = Datastr.Info.subjmass;
trial = Datastr.Info.Trial;

if isempty(subjsex)||isempty(subjmass)||isnan(subjsex)||isnan(subjmass)
    warning('M4_getRJCI:missingInfo',['Missing subject sexe and/or mass in trial ' trial]);
    return;
end

%% Do stuff

% Get joints, COM positions and inertia tensors (only applicable to humans)
[Datastr] = getRJCI(Datastr,subjsex,subjmass);


end