function [C3DProc] =  procC3D_mod(C3DRaw,MTrans,MTrans2,rotax,rotang)
% ProcC3D
% 
% Interpolate the data and swap dimensions
% 
% INPUT)
% C3DRaw : structure, output from readC3D
% MTrans : matrix, 3x3 transformation matrix for marker data. 
% Default=eye(3)
% 
% OUTPUT)
% C3DProc : structure, containing interpolated and rotated marker data, 
% as well as raw analog data
% 
% NOTES)
%It allows rot in final coordinate frame from biomechpro. Only around x!.
% TO DO: include y and z rotation and translation
% TO DO : check velocities, remove all which are too high
% TO DO : remove 1 sample before / after each missing data

if nargin == 1
    MTrans = eye(3);
end

%% Process MarkerData

% Check analog field 
if ~isfield(C3DRaw,'Analog')
    C3DRaw.Analog.AnalogData = 0;
    C3DRaw.Analog.AnalogFrameRate = 0;
end

% Get relevant data
MarkerData = C3DRaw.Marker.MarkerData;
VideoFrameRate = C3DRaw.Marker.VideoFrameRate;
AnalogData = C3DRaw.Analog.AnalogData;
AnalogFrameRate = C3DRaw.Analog.AnalogFrameRate;
ParameterGroup = C3DRaw.ParameterGroup;

% Browse the header for label information
units = [];
for i = 1:size(ParameterGroup, 2)
    if strcmp(char(ParameterGroup(i).name), 'POINT')
        for j = 1:size(ParameterGroup(i).Parameter, 2)
            if strcmp(char(ParameterGroup(i).Parameter(j).name), 'LABELS')
                MLabelCell = ParameterGroup(i).Parameter(j).data;
            end
            if strcmp(char(ParameterGroup(i).Parameter(j).name), 'UNITS')
                units = char(ParameterGroup(i).Parameter(j).data);
            end
        end
    end
end

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
         for irot=1:length(rotang)
             ax=rotax{irot};
             ang=rotang(irot);
             switch ax
                 case 'X'
                     MRot=rotx(ang);
                 case 'Y'
                     MRot=roty(ang);
                 case 'Z'
                     MRot=rotz(ang);
             end
             dum=[MRot*dum']';
         end
        dum=[MTrans2*dum']';
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

C3DProc.Marker.MarkerData = MarkerProc;
C3DProc.Marker.MarkerDataLabel = MLabelCell;
C3DProc.Marker.MarkerFrameRate = VideoFrameRate;
C3DProc.Marker = orderfields(C3DProc.Marker);

if ~isempty(AnalogData)
    
    C3DProc.Analog.AnalogData = AnalogData;
    C3DProc.Analog.AnalogFrameRate = AnalogFrameRate;
    C3DProc.Analog = orderfields(C3DProc.Analog);
    
end


end
