function [clusterMarkerLocal,usedLEDs] = putProbe(C3Ddata,cName)
%% Manually insert probe point using data viewer
%
% INPUT)
% C3Ddata : C3Ddata structure. Must be in same coordinates and scale as all other data from which you obtain local marker information.
% clusterLabel : string, containing the cluster name to which the probe position belongs
% 
% OUTPUT)
% clusterMarkerLocal : matrix containing probe position in local marker frame, to be stored in Body structure for later use.
% usedLEDs : vector containing information about which LEDs of the cluster are used in the reconstruction.
% 
% NOTES)
% Requires soder.m to work
% Requires viewC3D.m to work
% Once you've set the probe point, close the GUI when you're satisfied.
% 
% No check is made on the availability of the cluster markers!
% It is assumed all physical markers are visible in the provided C3Ddata

%% Check cluster data 

% Extract markers of the relevant cluster (can be more than 3!)
clusterMarkerGlobal = C3Ddata.Marker.MarkerData(:, strncmpi(cName,C3Ddata.Marker.MarkerDataLabel,length(cName))==1 ,:);
if size(clusterMarkerGlobal,2) < 3
    error(['getJointsLocal: Insufficient markers with name ' cName],'insufficient markers');
end

% WARNING: NO CHECK ON AVAILABILITY OF THESE LEDS !
usedLEDs = (1:size(clusterMarkerGlobal,2))';

%% Specify object sizes (pixels)

% Pushbutton size
pbutWi = 150;
pbutHi = 25;

% Input field size 
inputWi = 200;
inputHi = 25;

%% Visualise data

% Open viewer
viewData(C3Ddata);

% Get figure handle
h_gui = figure(1337);

%% Insert UI objects

% Going to insert some fields above the main figure axes
axPos = get(gca,'position');

% Insert coordinate field
h_pcoord = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',[1 1 1],...
    'position',[axPos(1)+axPos(3)-pbutWi-inputWi axPos(2)+axPos(4) inputWi inputHi],...
    'style','edit',...
    'callback',{},...
    'string','x,y,z -> press setProbe');

% Insert pushbuttons
h_setprobe = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[axPos(1)+axPos(3)-pbutWi axPos(2)+axPos(4) pbutWi pbutHi],...
    'callback',{@viewC3Dfun_setprobe},...
    'style','pushbutton',...
    'backgroundcolor',[1 0.5 0],...
    'string','SET PROBE');

h_done = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[axPos(1)+axPos(3) axPos(2)+axPos(4) pbutWi pbutHi],...
    'callback',{},...
    'style','togglebutton',...
    'backgroundcolor',[1 0 0],...
    'string','CLICK IF DONE');

% Create probe point plot object
h_insert = plot3(NaN,NaN,NaN,'o','color',[1 0.5 0],'markerfacecolor',[1 0.5 0]);

%% Put handles in userdata

h_all = get(h_gui,'userdata');

h_all.h_setprobe = h_setprobe;
h_all.h_done = h_done;
h_all.h_pcoord = h_pcoord;
h_all.h_insert = h_insert;

set(h_gui,'userdata',h_all);

%% Pause

% Find the "Done" button we just created
h_child = get(h_gui,'children');
h_done = h_child(1);  % NOTE: WILL BE THE FIRST, AS ITS THE LAST OBJECT WE CREATED!

% Stop function execution until user presses the "Done" button
waitfor(h_done,'value',1);

xdat = get(h_all.h_insert,'xdata');
ydat = get(h_all.h_insert,'ydata');
zdat = get(h_all.h_insert,'zdata');

close(h_gui);

%% Get output data

if isnan(xdat) || isnan(ydat) || isnan(zdat)
    error('putProbe:probeValue','No probe point set');
end

% Probe tip in global coordinates
TipGlobal = [xdat ydat zdat];

% Get 10 evenly spaced indices
idx_vis10 = floor( linspace(size(clusterMarkerGlobal,1)/10,size(clusterMarkerGlobal,1) - size(clusterMarkerGlobal,1)/10 ,10) );  % 10 evenly spaced idx for median value

clusterMarkerLocal = zeros(size(clusterMarkerGlobal,2),3,length(idx_vis10));
for igettip = 1:length(idx_vis10)     

    % Determine local tip position with respect to relevant frame
    % Contains the probed position w.r.t. the markers in local coordinates
    clusterMarkerLocal(:,:,igettip) = squeeze(clusterMarkerGlobal(idx_vis10(igettip),:,:)) - repmat(TipGlobal,[size(clusterMarkerGlobal,2),1]);  % nmark,dim,igettip
end

clusterMarkerLocal = median(clusterMarkerLocal,3);

%% Callback functions

    % Inserted pushbutton
    function viewC3Dfun_setprobe(hObject,eventdata)
        
        % Get all handles
        h_all = get(1337,'userdata');
        
        % Get coordinates from input field
        instr = get(h_all.h_pcoord,'string');
        
        [xcoord,rem] = strtok(instr,';, ');
        [ycoord,rem] = strtok(rem,';, ');
        [zcoord,~] = strtok(rem,';, ');
        
        % Plot point with coordinates
        try
            eval(['set(h_all.h_insert,''xdata'',' xcoord ',''ydata'',' ycoord ',''zdata'',' zcoord ');']);
        catch
            set(h_all.h_pcoord,'string','NaN,NaN,NaN');
            set(h_all.h_insert,'xdata',NaN','ydata',NaN,'zdata',NaN);
            warning('putProbe:invalidinput','Invalid input coordinates');
        end
        
    end

end