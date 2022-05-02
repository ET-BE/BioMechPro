function [Datastr] = TestFun(Datastr,inputA,inputB,inputC)
% gBMPDynUI inputA=1; inputB=7; inputC=3;
% 
% When creating a module function:
% - The first input must always be the data structure
% - There must always be one output: the data structure
% - You can pass an empty output. Then the data will not be saved.
% - No other outputs may be provided
% - You cannot use varargin or varargout
% 
% When integrating the module in the UI:
% - You can use a SINGLE line of text starting with gBMP DynUI* as above
% to  dynamically create UI fields for each input. You can specify the 
% number of lines to display for every input.
% - Every input in addition to Datastr must be specified, not more not less
% - The number specifies the number of lines of the edit field.
% - There MUST be a semicolon at the end of every specification.
% 
% You could try loading this function as a module in the UI, and edit
% the line on top to see how this affects the UI (requires
% closing and reloading the module).
% 
% * Here there's a space inserted because this gBMPDyn... term may only
% occur once in the document.

disp(['foo = ' num2str( inputA + inputB + inputC) ]);

end