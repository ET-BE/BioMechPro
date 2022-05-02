function [InsoleProc] =  procInsole(InsoleRaw, leftID, rightID)
% ProcInsole
% 
% Split the insole data and orgnize dimensions
% 
% INPUT)
% InsoleRaw : structure, the importted insole txt file data
% leftID: integer, the number ID of the left insole data. Default: 0
% rightID: integer, ther number ID of the right insole data. Default: 1
% 
% OUTPUT)
% InsoleProc : structure, containing the splitted and orgnized data, 
% as well as raw analog data
% 

%% Init

if nargin == 1
   leftID = 0;
   rightID = 1;
end

%% split left and right insole data

leftInsole = InsoleRaw.data(InsoleRaw.data(:, 2) == leftID, [1, 3:length(InsoleRaw.data(1, :))]);
leftInsole(:, 1) = leftInsole(:, 1) - leftInsole(1, 1);

rightInsole = InsoleRaw.data(InsoleRaw.data(:, 2) == rightID, [1, 3:length(InsoleRaw.data(1, :))]);
rightInsole(:, 1) = rightInsole(:, 1) - rightInsole(1, 1);

% save to InsoleProc strcture

InsoleProc.left = leftInsole;
InsoleProc.right = rightInsole;

end






