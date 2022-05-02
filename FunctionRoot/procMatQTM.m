function [DatastrProc] =  procMatQTM(DatastrRaw,MTrans)
% procMatQTM
% 
% Interpolate the data and swap dimensions
% 
% INPUT)
% DatastrRaw : structure, output from readC3D
% MTrans : matrix, 3x3 transformation matrix for marker data. 
% Default=eye(3)
% 
% OUTPUT)
% DatastrProc : structure, containing interpolated and rotated marker data, 
% as well as raw analog data
% 
% NOTES)

% TO DO : check velocities, remove all which are too high
% TO DO : remove 1 sample before / after each missing data

if nargin == 1
    MTrans = eye(3);
end

%% Process MarkerData

% Check analog field 
if ~isfield(DatastrRaw,'Analog')
    DatastrRaw.Analog.AnalogData = 0;
    DatastrRaw.Analog.AnalogFrameRate = 0;
end

% Get relevant data
MarkerData = DatastrRaw.Marker.MarkerData;
VideoFrameRate = DatastrRaw.Marker.VideoFrameRate;
AnalogData = DatastrRaw.Analog.AnalogData;
AnalogFrameRate = DatastrRaw.Analog.AnalogFrameRate;
ParameterGroup = DatastrRaw.ParameterGroup;

% % Browse the header for label information
units = [];
% for i = 1:size(ParameterGroup, 2)
%     if strcmp(char(ParameterGroup(i).name), 'POINT')
%         for j = 1:size(ParameterGroup(i).Parameter, 2)
%             if strcmp(char(ParameterGroup(i).Parameter(j).name), 'LABELS')
%                 MLabelCell = ParameterGroup(i).Parameter(j).data;
%             end
%             if strcmp(char(ParameterGroup(i).Parameter(j).name), 'UNITS')
%                 units = char(ParameterGroup(i).Parameter(j).data);
%             end
%         end
%     end
% end
MLabelCell = DatastrRaw.ParameterGroup.labels;
units=DatastrRaw.ParameterGroup.units;  % for now assuming that is always filled in

% Transform and Interpolate
MarkerProc = zeros(length(MarkerData),length(MLabelCell),3);
wrncheck = 0;
for iLabel=1:length(MLabelCell)
    if ~any( sum( abs(MarkerData(:,iLabel,:)) <= eps , 1) == size(MarkerData,1) ) %check for completely empty marker
        
        % Get specific marker data, convert to meters
        if strcmpi(units,'mm')
            dum=squeeze(MarkerData(:,iLabel,:))/1000;
        elseif strcmpi(units,'m')
            dum=squeeze(MarkerData(:,iLabel,:));
        else % Assume mm
            if wrncheck == 0
                warning('procC3D:units','Unknown units for marker data');
                wrncheck = 1;
            end
            dum=squeeze(MarkerData(:,iLabel,:))/1000;
        end
        
        % Transform to lab coordinate axis
        dum=[MTrans*dum']';
        
        % Interpolate
%         totInterpolation = sum(abs(dum(:,1))<=eps); %total interpolated samples
        try
            dum = interp1(find(abs(dum(:,1))>=eps),dum(abs(dum(:,1))>=eps,:),1:size(dum,1),'spline');
        catch
            disp(['Warning: Could not interpolate ' MLabelCell{iLabel}]);
%             totInterpolation = 0;
        end
        
%         if totInterpolation>10
%             disp(sprintf('Warning: In marker %s a total of %d samples is interpolated', MLabelCell{iLabel}, totInterpolation-1)); 
%         end
        
        MarkerProc(:,iLabel,:) = dum;
        
    end
end

%% Generate output

DatastrProc.Marker.MarkerData = MarkerProc;
DatastrProc.Marker.MarkerDataLabel = MLabelCell;
DatastrProc.Marker.MarkerFrameRate = VideoFrameRate;
DatastrProc.Marker = orderfields(DatastrProc.Marker);

if ~isempty(AnalogData)
    
    DatastrProc.Analog.AnalogData = AnalogData;
    DatastrProc.Analog.AnalogFrameRate = AnalogFrameRate;
    DatastrProc.Analog = orderfields(DatastrProc.Analog);
    
end


end
