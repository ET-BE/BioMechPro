function [Datastr] = Gen_rmAll(Datastr)
% Delete all data inside the structure, except for the most basic info

subjroot = Datastr.Info.SubjRoot;
trial = Datastr.Info.Trial;

Datastr = [];
Datastr.Info.SubjRoot = subjroot;
Datastr.Info.Trial = trial;

end