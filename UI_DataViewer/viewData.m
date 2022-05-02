function [] = viewData(Datastr,varargin)
%% viewC3D
% Open GUI to view marker and/or force data in Datastr structure
% 
% INPUT)
% Datastr : structure format with at least the fields:
%     .Marker.MarkerData
%     .Marker.MarkerDataLabel
%     .Marker.MarkerFrameRate
% 
% connectTree : cell containing names of variables stored in ...Label fields 
% within the .Marker field. Datapoints with these names will be connected with
% a line, which allows you to make stick figures. For example:
% connectTree = {'a','b','c'};
% will connect the points with names a, b, and c
% 
% 
% OUTPUT)
% Data viewer GUI
% 
% 
% NOTES)
% If a Force field is available in the structure, it can be plotted.
% Force data is assumed to have the same coordinate system as the marker data
% 
% For single plate data, the assumed column order is Fx, Fy, Fz, Mx, My, Mz
% For double plate data, the assumed column order is FxL, FyL, FzL, MxL, MyL, MzL, 
% followed by the same order for the right plate.
%
% To Do : add an analog data viewer
% 
% Last edit(s):
% Mark Vlutters - October 2014 (finalized v1)
% Mark Vlutters - January 2015 (added "set" buttons for slider start/stop)

try
    close(1337)
end

%% Check inputs

if nargin == 2
    connectTree = varargin{1};
else
    connectTree = [];
end

% Check marker field
if ~isfield(Datastr,'Marker')
    error('viewC3D:No marker field available in structure'); 
end

% Check force field
if isfield(Datastr,'Force')
    isforce = true;
else
    isforce = false;
end

% Check event field
isgpm = false;
if isfield(Datastr,'Event');
    isgpm = isfield(Datastr.Event,'GaitPhaseM');
end

%% Collect information

% Marker list
labels = {};
labellen = zeros(4,1);
if isfield(Datastr.Marker,'MarkerDataLabel')
    labels = Datastr.Marker.MarkerDataLabel;
    labellen(1) = length(Datastr.Marker.MarkerDataLabel);
end
if isfield(Datastr.Marker,'ProbedDataLabel')
    labels(end+1:end+length(Datastr.Marker.ProbedDataLabel)) = Datastr.Marker.ProbedDataLabel;
    labellen(2) = length(Datastr.Marker.ProbedDataLabel);
end
if isfield(Datastr.Marker,'COMDataLabel')
    labels(end+1:end+length(Datastr.Marker.COMDataLabel)) = Datastr.Marker.COMDataLabel;
    labellen(3) = length(Datastr.Marker.COMDataLabel);
end
if isfield(Datastr.Marker,'JointDataLabel')
    labels(end+1:end+length(Datastr.Marker.JointDataLabel)) = Datastr.Marker.JointDataLabel;
    labellen(4) = length(Datastr.Marker.JointDataLabel);
end

% Force data pre-processing (optional)
if isforce
    if isfield(Datastr.Marker,'MarkerFrameRate') && isfield(Datastr.Force,'ForceFrameRate')
        % Get relevant sample frequencies and offsets
        fs_v = Datastr.Marker.MarkerFrameRate;
        fs_f = Datastr.Force.ForceFrameRate;
    
        % Factor for normalizing vector ( 1000 N <> 2 m )
        normfact = 500;
        
        % Indices for syncing
        if isfield(Datastr.Marker,'MarkerSyncIdx') && isfield(Datastr.Force,'ForceSyncIdx')
            idxSync_v = Datastr.Marker.MarkerSyncIdx;
            idxSync_f = Datastr.Force.ForceSyncIdx;
        else
            idxSync_v = 1;
            idxSync_f = 1;
            warning('viewC3D:nosync','Data might be presented out of sync.');
        end

        % Instantiate Force data container of same size as marker data
        fdata = zeros( size(Datastr.Marker.MarkerData,1), size(Datastr.Force.ForceData,2) );

        % Resample actual force data
        fdata_r = resample( Datastr.Force.ForceData(idxSync_f:end,:) , fs_v , fs_f );

        % Put the data in the container at the appropriate shift
        if length( idxSync_v:size(fdata,1) ) > length(fdata_r)
            fdata( idxSync_v:idxSync_v+size(fdata_r,1)-1 , :) = fdata_r;
        else % if length( idxSync_v:size(fdata,1) ) <= length(fdata_r)
            fdata( idxSync_v:end , :) = fdata_r( 1:length(idxSync_v:size(fdata,1)) , :);
        end

        % Find COPs and create force vectors
        % Start position is the COP, end is at vector length
        if size(Datastr.Force.ForceData,2) == 6
            COPdata = getCOP(fdata,25); % Fthresh = 25

            % TODO: make small adjustment to allow time-varying offsets
            % NOTE: offset addition for COPT only works if both have the same origin
            if isfield(Datastr.Force,'ForcePlateOffset') 
                forcePlateOffset = Datastr.Force.ForcePlateOffset; % Lxyz Rxyz
            else
                forcePlateOffset = zeros(1,6);
            end
            
            % Get force vectors                         % >> currently treated as global data, not passed to functions
            fvec_T(:,:,1) = COPdata.COPT + repmat(forcePlateOffset,[size(COPdata.COPL,1) 1]);
            fvec_T(:,:,2) = COPdata.COPT + fdata(:,1:3) ./ normfact + repmat(forcePlateOffset(1:3),[size(COPdata.COPL,1) 1]);
            fvec_L = zeros(size(fvec_T));
            fvec_R = zeros(size(fvec_T));
           
        elseif size(Datastr.Force.ForceData,2) == 12
            COPdata = getCOP(fdata,25); % Fthresh = 25
            
            % TODO: make small adjustment to allow time-varying offsets
            % NOTE: offset addition for COPT only works if both have the same origin
            if isfield(Datastr.Force,'ForcePlateOffset') 
                forcePlateOffset = Datastr.Force.ForcePlateOffset; % Lxyz Rxyz
            else
                forcePlateOffset = zeros(1,6);
            end
            
            fvec_T(:,:,1) = COPdata.COPT + repmat(forcePlateOffset(1:3),[size(COPdata.COPL,1) 1]); % Assumed left plate offset same as right
            fvec_T(:,:,2) = COPdata.COPT + ( fdata(:,1:3) + fdata(:,7:9) ) ./ normfact + repmat(forcePlateOffset(1:3),[size(COPdata.COPL,1) 1]);
            
            fvec_L(:,:,1) = COPdata.COPL + repmat(forcePlateOffset(1:3),[size(COPdata.COPL,1) 1]);
            fvec_L(:,:,2) = COPdata.COPL + fdata(:,1:3) ./ normfact + repmat(forcePlateOffset(1:3),[size(COPdata.COPL,1) 1]);
            
            fvec_R(:,:,1) = COPdata.COPR + repmat(forcePlateOffset(4:6),[size(COPdata.COPL,1) 1]);
            fvec_R(:,:,2) = COPdata.COPR + fdata(:,7:9) ./ normfact + repmat(forcePlateOffset(4:6),[size(COPdata.COPL,1) 1]);
            
        end
    
    else
        warning('viewC3D:framerate','No frame rate information available to show force data.')
        
        % Make zero force data
        fvec_T = zeros(size(Datastr.Marker.MarkerData,1),3,2);
        fvec_L = fvec_T;
        fvec_R = fvec_T;
        
    end

end

% Determine slider step size from data
datlen = size(Datastr.Marker.MarkerData,1);
if datlen <= 1e6
    minorstep = 1/datlen;
else
    minorstep = 1e-6; % minimum
    warning('viewC3D:nosync','Cannot go lower than slider minorstep size.');
end
if minorstep < 0.01
    majorstep = 0.01;
else
    majorstep = minorstep;
end

% Store info
h_all.labellen = labellen;
h_all.idx = 1;              % default sample index
h_all.idxhigh = 1;          % default highlighted marker number

% Movie object
h_all.h_mov = [];

%% Figure
scrsize = get(0,'screensize');
figwi = 950; 
fighi = 800;

h_gui = figure(1337);
set(h_gui,...
    'name','viewC3D',...
    'numbertitle','off',...
    'menubar','none',...
    'toolbar','figure',...
    'position',[(scrsize(3)-figwi)./2 (scrsize(4)-fighi)./2 figwi fighi],...
    'resize','off');

h_all.h_gui = h_gui;

%% Axis 
axwi = fighi - 100;
axhi = fighi - 100;
h_ax = axes;
set(h_ax,...
    'units','pixels',...
    'color',[1 1 1],...
    'fontsize',8,...
    'position',[(fighi-axwi)/2 (fighi-axhi)./1.5 axwi axhi],...
    'visible','on',...
    'xlim',[-1 1],'ylim',[-1 1],'zlim',[0 2],...
    'plotboxaspectratiomode','manual',... 
    'view',[0 0],...
    'xgrid','on','ygrid','on','zgrid','on');

h_all.h_ax = h_ax;


%% Data pre-plot
% Pre-plot data and set data source (to update later with refreshdata)
hold on;

% Basic LED data
if labellen(1) ~= 0
    data_mark_x = squeeze(Datastr.Marker.MarkerData(1,:,1));
    data_mark_y = squeeze(Datastr.Marker.MarkerData(1,:,2));
    data_mark_z = squeeze(Datastr.Marker.MarkerData(1,:,3));
    
    h_mark = plot3(data_mark_x,data_mark_y,data_mark_z,'ok','markerfacecolor','k');
    
    set(h_mark,'xdatasource','data_mark_x',...
               'ydatasource','data_mark_y',...
               'zdatasource','data_mark_z');
   
   h_all.h_mark = h_mark;

    % Set default index range 
    idxrange = [1 size(Datastr.Marker.MarkerData,1)];
    h_all.idxrange = idxrange;
end

% Probed data
if labellen(2) ~= 0
    data_prob_x = squeeze(Datastr.Marker.ProbedData(1,:,1));
    data_prob_y = squeeze(Datastr.Marker.ProbedData(1,:,2));
    data_prob_z = squeeze(Datastr.Marker.ProbedData(1,:,3));
    
    h_prob = plot3(data_prob_x,data_prob_y,data_prob_z,'or','markerfacecolor','r');
    
    set(h_prob,'xdatasource','data_prob_x',...
               'ydatasource','data_prob_y',...
               'zdatasource','data_prob_z',...
               'visible','off');  % Invisible by default
           
    h_all.h_prob = h_prob;
    
end

% Segment COM data
if labellen(3) ~= 0
    data_com_x = squeeze(Datastr.Marker.COMData(1,:,1));
    data_com_y = squeeze(Datastr.Marker.COMData(1,:,2));
    data_com_z = squeeze(Datastr.Marker.COMData(1,:,3));
    
    h_com = plot3(data_com_x,data_com_y,data_com_z,'ob','markerfacecolor','b');
    
    set(h_com,'xdatasource','data_com_x',...
              'ydatasource','data_com_y',...
              'zdatasource','data_com_z',...
              'visible','off');  % Invisible by default
          
    h_all.h_com = h_com;

end

% Joint data
if labellen(4) ~= 0
    data_joint_x = squeeze(Datastr.Marker.JointData(1,:,1));
    data_joint_y = squeeze(Datastr.Marker.JointData(1,:,2));
    data_joint_z = squeeze(Datastr.Marker.JointData(1,:,3));
    
    h_joint = plot3(data_joint_x,data_joint_y,data_joint_z,'om','markerfacecolor','m');
    
    set(h_joint,'xdatasource','data_joint_x',...
                'ydatasource','data_joint_y',...
                'zdatasource','data_joint_z',...
                'visible','off');  % Invisible by default
            
    h_all.h_joint = h_joint;

end

% Highlight marker data
if any(labellen ~= 0)
    
    data_high_x = 0;
    data_high_y = 0;
    data_high_z = 0;

    h_high = plot3(data_high_x,data_high_y,data_high_z,'og','linewidth',3,'markersize',10);
    
    set(h_high,'xdatasource','data_high_x',...
               'ydatasource','data_high_y',...
               'zdatasource','data_high_z',...
               'visible','off');

    h_all.h_high = h_high;
    
end

% Force data
if isforce
    data_fT_x = squeeze(fvec_T(1,1,:));
    data_fT_y = squeeze(fvec_T(1,2,:));
    data_fT_z = squeeze(fvec_T(1,3,:));

    data_fL_x = squeeze(fvec_L(1,1,:));
    data_fL_y = squeeze(fvec_L(1,2,:));
    data_fL_z = squeeze(fvec_L(1,3,:));

    data_fR_x = squeeze(fvec_R(1,1,:));
    data_fR_y = squeeze(fvec_R(1,2,:));
    data_fR_z = squeeze(fvec_R(1,3,:));

    h_fT = plot3(data_fT_x,data_fT_y,data_fT_z,'-','color',[0 0.8 0],'linewidth',2);
    h_fL = plot3(data_fL_x,data_fL_y,data_fL_z,'-','color',[0 0 0.8],'linewidth',2);
    h_fR = plot3(data_fR_x,data_fR_y,data_fR_z,'-','color',[0.8 0 0],'linewidth',2);

    set(h_fT,'xdatasource','data_fT_x',...
             'ydatasource','data_fT_y',...
             'zdatasource','data_fT_z',...
             'visible','off');  % Invisible by default
    
    set(h_fL,'xdatasource','data_fL_x',...
             'ydatasource','data_fL_y',...
             'zdatasource','data_fL_z',...
             'visible','off');  % Invisible by default
         
     set(h_fR,'xdatasource','data_fR_x',...
              'ydatasource','data_fR_y',...
              'zdatasource','data_fR_z',...
              'visible','off');  % Invisible by default
         
    h_all.h_fT = h_fT;
    h_all.h_fL = h_fL;
    h_all.h_fR = h_fR;
    
end

% Stick data
if ~isempty(connectTree)
    
    stick_x = 0;
    stick_y = 0;
    stick_z = 0;
    
    h_stick = plot3(stick_x,stick_y,stick_z,'-','color',[0 0 0],'linewidth',3);
    
    set(h_stick,'xdatasource','stick_x',...
        'ydatasource','stick_y',...
        'zdatasource','stick_z',...
        'visible','on');
    
    % Allocate
    stickData = zeros(size(Datastr.Marker.MarkerData,1),1,3);
    
    % Get relevant data points
    istore = 0;
    for inode = 1:length(connectTree)

        % Determine the set which the label is in
        idxmark = find(strcmpi(labels,connectTree{inode}));
        
        if ~isempty(idxmark)
            
            idxset = find(idxmark <= cumsum(labellen),1,'first');

            if ~isempty(idxset)
                istore = istore + 1;
            end
            
            if idxset == 1
                stickData(:,istore,:) = Datastr.Marker.MarkerData(:,idxmark,:);
            elseif idxset == 2
                stickData(:,istore,:) = Datastr.Marker.ProbedData(:,idxmark-labellen(1),:);
            elseif idxset == 3
                stickData(:,istore,:) = Datastr.Marker.COMData(:,idxmark-labellen(1)-labellen(2),:);
            elseif idxset == 4
                stickData(:,istore,:) = Datastr.Marker.JointData(:,idxmark-labellen(1)-labellen(2)-labellen(3),:);
            end
        end
    end
    
    % Store data
    h_all.h_stick = h_stick;
    
end


%% uicontrol object creation

% Time indicator
timhi = 25;
timwi = 100;
h_time = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[0 fighi-timhi timwi timhi],...
    'style','text',...
    'string','t = 0',...
    'fontsize',10,...
    'fontweight','bold');

h_all.h_time = h_time;

% Gait phase indicators, display the phasevector value if available
if isgpm
    h_gpm = uicontrol('parent',h_gui,...
        'units','pixels',...
        'backgroundcolor',get(h_gui,'color'),...
        'position',[timwi fighi-timhi timwi timhi],...
        'style','text',...
        'string','GPm = 0',...
        'fontsize',10,...
        'fontweight','bold');

    h_all.h_gpm = h_gpm;
end

% Slider bar properties
barwi = axwi;
barhi = 25;

% Play slider bar
h_bar1 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[(fighi-axwi)/2 ((fighi-axhi)/6) barwi barhi],...   %(fighi-axhi)/5
    'callback',{@viewC3Dfun_slider1,Datastr},...
    'style','slider',...
    'sliderstep',[minorstep,majorstep]);

h_all.h_bar1 = h_bar1;

% Slider bars for start and stop sample
h_bar2 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[(fighi-axwi)/2 (fighi-axhi)/6+barhi barwi barhi/3],...
    'callback',{@viewC3Dfun_slider2,Datastr},...
    'style','slider',...
    'sliderstep',[minorstep,majorstep],...
    'backgroundcolor',[0.8 0.8 1],...
    'value',0);

h_bar3 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[(fighi-axwi)/2 (fighi-axhi)/6-barhi/3 barwi barhi/3],...
    'callback',{@viewC3Dfun_slider3,Datastr},...
    'style','slider',...
    'sliderstep',[minorstep,majorstep],...
    'backgroundcolor',[0.8 0.8 1],...
    'value',1);

h_all.h_bar2 = h_bar2;
h_all.h_bar3 = h_bar3;

% Slider bars set button (for start / stop sample)
h_set2 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[(fighi-axwi)/2 - 25 (fighi-axhi)/6+0.5*barhi 25 (5/6).*barhi],...
    'callback',{@viewC3Dfun_set2,Datastr},...
    'style','pushbutton',...
    'backgroundcolor',[0.7 0.7 1],...
    'string','SET');
h_set3 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[(fighi-axwi)/2 - 25 (fighi-axhi)/6-(1/3)*barhi 25 (5/6).*barhi],...
    'callback',{@viewC3Dfun_set3},...
    'style','pushbutton',...
    'backgroundcolor',[0.7 0.7 1],...
    'string','SET');

h_all.h_set2 = h_set2;
h_all.h_set3 = h_set3;

% listbox for labelcells
boxwi = 150;
boxhi = axhi * 0.45;

uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + axhi-boxhi boxwi boxhi],...
    'style','text',...
    'string','-- Highlight Data --',...
    'fontsize',10,...
    'fontweight','bold');

h_lbl = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[axwi+75 (fighi-axhi)./1.5 + axhi-boxhi boxwi boxhi-25],...
    'backgroundcolor',[1 1 1],...
    'callback',{@viewC3Dfun_labellist,Datastr},...
    'string',labels,...
    'style','listbox',...
    'max', sum(labellen) ,...
    'fontsize',10);

% Plot button
pbutwi = 150;
pbuthi = barhi;

h_plot = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',[0.8 0.8 1],...
    'position',[axwi+75 (fighi-axhi)./1.5 + axhi-boxhi-pbuthi pbutwi pbuthi],...
    'callback',{@viewC3Dfun_plot,Datastr},...
    'style','pushbutton',...
    'string','PLOT',...
    'fontweight','bold');

% Movie button
h_mov = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',[0.8 1 0.8],...
    'position',[axwi+75 (fighi-axhi)./1.5 + axhi-boxhi-2.*pbuthi pbutwi pbuthi],...
    'callback',{@viewC3Dfun_mov,Datastr},...
    'style','togglebutton',...
    'string','RECORD',...
    'fontweight','bold');


% Checkboxes for data visibilty
chckwi = boxwi;
chckhi = 25;

enable1 = 'off'; enable2 = 'off'; 
enable3 = 'off'; enable4 = 'off';
enable5 = 'off';
if labellen(1) ~= 0
    enable1 = 'on';
end
if labellen(2) ~= 0
    enable2 = 'on';
end
if labellen(3) ~= 0
    enable3 = 'on';
end
if labellen(4) ~= 0
    enable4 = 'on';
end
if isforce
    enable5 = 'on';
end

uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 3.8*chckhi chckwi chckhi],...
    'style','text',...
    'string','-- Show Position Data --',...
    'fontsize',10,...
    'fontweight','bold');

h_chck1 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 4.5*chckhi chckwi chckhi],...
    'callback',{@viewC3Dfun_chck1},...
    'style','checkbox',...
    'value',1,...
    'string','Marker',...
    'fontsize',10,...
    'enable',enable1);

h_chck2 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+150 (fighi-axhi)./1.5 + (axhi-boxhi) - 4.5*chckhi chckwi chckhi],...
    'callback',{@viewC3Dfun_chck2},...
    'style','checkbox',...
    'value',0,...
    'string','Probe',...
    'fontsize',10,...
    'enable',enable2);

h_chck3 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 5.5*chckhi chckwi chckhi],...
    'callback',{@viewC3Dfun_chck3},...
    'style','checkbox',...
    'value',0,...
    'string','COM',...
    'fontsize',10,...
    'enable',enable3);

h_chck4 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+150 (fighi-axhi)./1.5 + (axhi-boxhi) - 5.5*chckhi chckwi chckhi],...
    'callback',{@viewC3Dfun_chck4},...
    'style','checkbox',...
    'value',0,...
    'string','Joint',...
    'fontsize',10,...
    'enable',enable4);

uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 6.8*chckhi chckwi chckhi],...
    'style','text',...
    'string','--- Show Force Data ---',...
    'fontsize',10,...
    'fontweight','bold');

h_chck5 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 7.5*chckhi chckwi/3 chckhi],...
    'callback',{@viewC3Dfun_chck5},...
    'style','checkbox',...
    'value',0,...
    'string','All',...
    'fontsize',10,...
    'enable',enable5);

h_chck6 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+125 (fighi-axhi)./1.5 + (axhi-boxhi) - 7.5*chckhi chckwi/3 chckhi],...
    'callback',{@viewC3Dfun_chck6},...
    'style','checkbox',...
    'value',0,...
    'string','Left',...
    'fontsize',10,...
    'enable',enable5);

h_chck7 = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+175 (fighi-axhi)./1.5 + (axhi-boxhi) - 7.5*chckhi chckwi/3 chckhi],...
    'callback',{@viewC3Dfun_chck7},...
    'style','checkbox',...
    'value',0,...
    'string','Right',...
    'fontsize',10,...
    'enable',enable5);

h_all.h_chck1 = h_chck1; h_all.h_chck2 = h_chck2; 
h_all.h_chck3 = h_chck3; h_all.h_chck4 = h_chck4; 
h_all.h_chck5 = h_chck5; h_all.h_chck6 = h_chck6; h_all.h_chck7 = h_chck7;

% Field of view fields
fovwi = chckwi*0.5;
fovhi = chckhi;
txtwi = chckwi*0.5;
txtoffset = 3; % pixels

uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 9.5*chckhi - txtoffset chckwi fovhi],...
    'style','text',...
    'string','-- Axis Options --',...
    'fontsize',10,...
    'fontweight','bold');

h_xlim = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',[1 1 1],...
    'position',[axwi+75+txtwi (fighi-axhi)./1.5 + (axhi-boxhi) - 10.5*chckhi fovwi fovhi],...
    'style','edit',...
    'callback',{@viewC3Dfun_xlim},...
    'string','-1 1');

uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 10.5*chckhi - txtoffset txtwi fovhi],...
    'style','text',...
    'string','x-lim       : ',...
    'fontsize',10);

h_ylim = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',[1 1 1],...
    'position',[axwi+75+txtwi (fighi-axhi)./1.5 + (axhi-boxhi) - 11.5*chckhi fovwi fovhi],...
    'style','edit',...
    'callback',{@viewC3Dfun_ylim},...
    'string','-1 1');

uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 11.5*chckhi - txtoffset txtwi fovhi],...
    'style','text',...
    'string','y-lim       : ',...
    'fontsize',10);

h_zlim = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',[1 1 1],...
    'position',[axwi+75+txtwi (fighi-axhi)./1.5 + (axhi-boxhi) - 12.5*chckhi fovwi fovhi],...
    'style','edit',...
    'callback',{@viewC3Dfun_zlim},...
    'string','0 2');

uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',get(h_gui,'color'),...
    'position',[axwi+75 (fighi-axhi)./1.5 + (axhi-boxhi) - 12.5*chckhi - txtoffset txtwi fovhi],...
    'style','text',...
    'string','z-lim       : ',...
    'fontsize',10);

% Viewpoint buttons
vbuthi = 50;
vbutwi = 50;

h_viewxy = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[axwi+75 (fighi-axhi)./1.5 vbuthi vbutwi],...
    'callback',{@viewC3Dfun_viewxy},...
    'style','pushbutton',...
    'string','XY');

h_viewxz = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[axwi+75+vbutwi (fighi-axhi)./1.5 vbuthi vbutwi],...
    'callback',{@viewC3Dfun_viewxz},...
    'style','pushbutton',...
    'string','XZ');

h_viewyz = uicontrol('parent',h_gui,...
    'units','pixels',...
    'position',[axwi+75+2*vbutwi (fighi-axhi)./1.5 vbuthi vbutwi],...
    'callback',{@viewC3Dfun_viewyz},...
    'style','pushbutton',...
    'string','YZ');

% Play / pause button
pbutwi = 150;
pbuthi = barhi;

h_play = uicontrol('parent',h_gui,...
    'units','pixels',...
    'backgroundcolor',[1 0 0],...
    'position',[axwi+75 (fighi-axhi)./5 pbutwi pbuthi],...
    'callback',{@viewC3Dfun_play,Datastr},...
    'style','togglebutton',...
    'string','PLAY',...
    'fontweight','bold');


%% Store h_all in root figure
% Because you shouldn't try accessing output variables of callback functions
% I find this easier than having to look for the handles in each callback 

set(h_gui,'userdata',h_all);  % h_gui = 1337


%% uicontrol callback functions


    % Main Slider function 
    function viewC3Dfun_slider1(hObject,eventdata,C3Ddata)
        h_all = get(1337,'userdata');
                
        % Get marker index from slider
        idx = round( get(hObject,'value') * size(C3Ddata.Marker.MarkerData,1) );
        if idx == 0 % Because sliderbar goes from 0 to 100%
            idx = 1;
        end
        
        % Update time indicator
        set(h_all.h_time,'string',['t = ' num2str( idx ./ C3Ddata.Marker.MarkerFrameRate )] );
        
        % Update gait phase indicator (if possible)
        if isfield(h_all,'h_gpm')
            set(h_all.h_gpm,'string',['GPm = ' num2str(C3Ddata.Event.GaitPhaseM(idx,end)) ]);
        end
        
        % Update data in viewer
        if (h_all.labellen(1) ~= 0) && get(h_all.h_chck1,'value')
            data_mark_x = squeeze(C3Ddata.Marker.MarkerData(idx,:,1));
            data_mark_y = squeeze(C3Ddata.Marker.MarkerData(idx,:,2));
            data_mark_z = squeeze(C3Ddata.Marker.MarkerData(idx,:,3));
            refreshdata(h_all.h_mark,'caller');
        end
        if (h_all.labellen(2) ~= 0) && get(h_all.h_chck2,'value')
            data_prob_x = squeeze(C3Ddata.Marker.ProbedData(idx,:,1));
            data_prob_y = squeeze(C3Ddata.Marker.ProbedData(idx,:,2));
            data_prob_z = squeeze(C3Ddata.Marker.ProbedData(idx,:,3));
            refreshdata(h_all.h_prob,'caller');
        end
        if (h_all.labellen(3) ~= 0) && get(h_all.h_chck3,'value')
            data_com_x = squeeze(C3Ddata.Marker.COMData(idx,:,1));
            data_com_y = squeeze(C3Ddata.Marker.COMData(idx,:,2));
            data_com_z = squeeze(C3Ddata.Marker.COMData(idx,:,3));
            refreshdata(h_all.h_com,'caller');
        end
        if (h_all.labellen(4) ~= 0) && get(h_all.h_chck4,'value')
            data_joint_x = squeeze(C3Ddata.Marker.JointData(idx,:,1));
            data_joint_y = squeeze(C3Ddata.Marker.JointData(idx,:,2));
            data_joint_z = squeeze(C3Ddata.Marker.JointData(idx,:,3));
            refreshdata(h_all.h_joint,'caller');
        end
        
        % Force data
        if isfield(h_all,'h_fT')
            if get(h_all.h_chck5,'value')
                data_fT_x = squeeze(fvec_T(idx,1,:));
                data_fT_y = squeeze(fvec_T(idx,2,:));
                data_fT_z = squeeze(fvec_T(idx,3,:));
                refreshdata(h_all.h_fT,'caller');
            end
            if get(h_all.h_chck6,'value')
                data_fL_x = squeeze(fvec_L(idx,1,:));
                data_fL_y = squeeze(fvec_L(idx,2,:));
                data_fL_z = squeeze(fvec_L(idx,3,:));
                refreshdata(h_all.h_fL,'caller');
            end
            if get(h_all.h_chck7,'value')
                data_fR_x = squeeze(fvec_R(idx,1,:));
                data_fR_y = squeeze(fvec_R(idx,2,:));
                data_fR_z = squeeze(fvec_R(idx,3,:));
                refreshdata(h_all.h_fR,'caller');
            end
        end
        
        % Stickdata
        if isfield(h_all,'h_stick')
            stick_x = squeeze(stickData(idx,:,1));
            stick_y = squeeze(stickData(idx,:,2));
            stick_z = squeeze(stickData(idx,:,3));
            refreshdata(h_all.h_stick,'caller');
        end
        
        % Highlighted marker data
        idxhigh_all = h_all.idxhigh; % Note : these are ALL highlighted markers
        
        % Clear data_high (because it's filled up in the loop)
        data_high_x = zeros(1,length(idxhigh_all)); 
        data_high_y = zeros(1,length(idxhigh_all)); 
        data_high_z = zeros(1,length(idxhigh_all)); 
        
        % Update highlighted marker
        for idxhigh = idxhigh_all
            n = idxhigh == idxhigh_all; % indexing variable
            
            idxset = find(idxhigh <= cumsum(h_all.labellen),1,'first');
            if idxset == 1
                data_high_x(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,1));
                data_high_y(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,2));
                data_high_z(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,3));
            elseif idxset == 2
                data_high_x(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),1)); % Correction for the set which the marker belongs to
                data_high_y(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),2));
                data_high_z(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),3));
            elseif idxset == 3
                data_high_x(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1)+h_all.labellen(2) ),1));
                data_high_y(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1)+h_all.labellen(2) ),2));
                data_high_z(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1)+h_all.labellen(2) ),3));
            elseif idxset == 4
                data_high_x(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1)+h_all.labellen(2)+h_all.labellen(3) ),1));
                data_high_y(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1)+h_all.labellen(2)+h_all.labellen(3) ),2));
                data_high_z(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1)+h_all.labellen(2)+h_all.labellen(3) ),3));
            end
        end
        refreshdata(h_all.h_high,'caller');
    %         refreshdata([h_all.h_mark h_all.h_prob h_all.h_high],'caller'); % This doesn't work for some reasons
    
        h_all.idx = idx;
        set(1337,'userdata',h_all);
    end

    function viewC3Dfun_slider2(hObject,eventdata,C3Ddata)
        % Note: this is the top slider, bar2
        
        h_all = get(1337,'userdata');
        
        % Get idx values
        idxA = round( get(hObject,'value') * size(C3Ddata.Marker.MarkerData,1) );
        idxB = round( get(h_all.h_bar3,'value') * size(C3Ddata.Marker.MarkerData,1) );
        if idxA == 0
            idxA = 1;
        end
        if idxB == 0
            idxB = 1;
        end
        
        % Sort values and set range
        idxrange = sort([idxA idxB],'ascend');
        
        % Update handle
        h_all.idxrange = idxrange;
        set(1337,'userdata',h_all);
    end

    function viewC3Dfun_slider3(hObject,eventdata,C3Ddata)
        % Note: this is the bottom slider, bar3
        
        h_all = get(1337,'userdata');
        
        % Get idx values
        idxA = round( get(hObject,'value') * size(C3Ddata.Marker.MarkerData,1) );
        idxB = round( get(h_all.h_bar2,'value') * size(C3Ddata.Marker.MarkerData,1) );
        if idxA == 0
            idxA = 1;
        end
        if idxB == 0
            idxB = 1;
        end
        
        % Sort values and set range
        idxrange = sort([idxA idxB],'ascend');
        
        % Update handle
        h_all.idxrange = idxrange;
        set(1337,'userdata',h_all);
    end

    % Slider set pushbutton callbacks
    function viewC3Dfun_set2(hObject,eventdata,C3Ddata)
        h_all = get(1337,'userdata');
        
        % Set start slider to current main slider value
        set(h_all.h_bar2,'value',get(h_all.h_bar1,'value') );
        
        % Update idxrange by calling slider callback
        viewC3Dfun_slider2(h_all.h_bar2,eventdata,C3Ddata);
    end
    function viewC3Dfun_set3(hObject,eventdata)
        h_all = get(1337,'userdata');
        
        % Set start slider to current main slider value
        set(h_all.h_bar3,'value',get(h_all.h_bar1,'value') );
        
        % Update idxrange by calling slider callback
        viewC3Dfun_slider3(h_all.h_bar3,eventdata,Datastr);
    end

    % Listbox callback
    function viewC3Dfun_labellist(hObject,eventdata,C3Ddata)
        h_all = get(1337,'userdata');
        
        idx = h_all.idx;
        idxhigh_all = get(hObject,'value');
        
        % Switch on / off when clicking the same
        if length(idxhigh_all) == 1
            if strcmpi(get(h_all.h_high,'visible'),'on') && any( idxhigh_all == h_all.idxhigh ) && ( length(h_all.idxhigh) == 1 )
                set(h_all.h_high,'visible','off');
            elseif strcmpi(get(h_all.h_high,'visible'),'off')
                set(h_all.h_high,'visible','on');
            end
        else % Multiple selected!
            set(h_all.h_high,'visible','on');
        end
        
        % update data
        h_all.idxhigh = idxhigh_all;
        set(1337,'userdata',h_all);
        
        % Clear data_high (because it's filled up in the loop)
        data_high_x = zeros(1,length(idxhigh_all)); 
        data_high_y = zeros(1,length(idxhigh_all)); 
        data_high_z = zeros(1,length(idxhigh_all)); 
        
        % Update highlighted marker
        for idxhigh = idxhigh_all
            n = idxhigh == idxhigh_all; % indexing variable
            
            idxset = find(idxhigh <= cumsum(h_all.labellen),1,'first');
            if idxset == 1
                data_high_x(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,1));
                data_high_y(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,2));
                data_high_z(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,3));
            elseif idxset == 2
                data_high_x(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),1)); % Correction for the set which the marker belongs to
                data_high_y(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),2));
                data_high_z(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),3));
            elseif idxset == 3
                data_high_x(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2) ),1));
                data_high_y(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2) ),2));
                data_high_z(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2) ),3));
            elseif idxset == 4
                data_high_x(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2)+ h_all.labellen(3) ),1));
                data_high_y(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2)+ h_all.labellen(3) ),2));
                data_high_z(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2)+ h_all.labellen(3) ),3));
            end
        end
        refreshdata(h_all.h_high,'caller');
        
        % Stickdata
        if isfield(h_all,'h_stick')
            stick_x = squeeze(stickData(idx,:,1));
            stick_y = squeeze(stickData(idx,:,2));
            stick_z = squeeze(stickData(idx,:,3));
            refreshdata(h_all.h_stick,'caller');
        end
        
    end

    % Checkbox callbacks
    function viewC3Dfun_chck1(hObject,eventdata)
        h_all = get(1337,'userdata');
        checkval = get(hObject,'value');
        
        if checkval == 1
            set(h_all.h_mark,'visible','on');
        elseif checkval == 0
            set(h_all.h_mark,'visible','off');
        end
        
        set(1337,'userdata',h_all);
    end
    function viewC3Dfun_chck2(hObject,eventdata)
        h_all = get(1337,'userdata');
        checkval = get(hObject,'value');
        
        if checkval == 1
            set(h_all.h_prob,'visible','on');
        elseif checkval == 0
            set(h_all.h_prob,'visible','off');
        end
        
        set(1337,'userdata',h_all);
    end
    function viewC3Dfun_chck3(hObject,eventdata)
        h_all = get(1337,'userdata');
        checkval = get(hObject,'value');
        
        if checkval == 1
            set(h_all.h_com,'visible','on');
        elseif checkval == 0
            set(h_all.h_com,'visible','off');
        end
        
        set(1337,'userdata',h_all);
    end
    function viewC3Dfun_chck4(hObject,eventdata)
        h_all = get(1337,'userdata');
        checkval = get(hObject,'value');
        
        if checkval == 1
            set(h_all.h_joint,'visible','on');
        elseif checkval == 0
            set(h_all.h_joint,'visible','off');
        end
        
        set(1337,'userdata',h_all);
    end
    function viewC3Dfun_chck5(hObject,eventdata)
        h_all = get(1337,'userdata');
        checkval = get(hObject,'value');
        
        if checkval == 1
            set(h_all.h_fT,'visible','on');
        elseif checkval == 0
            set(h_all.h_fT,'visible','off');
        end
        
        set(1337,'userdata',h_all);
    end
    function viewC3Dfun_chck6(hObject,eventdata)
        h_all = get(1337,'userdata');
        checkval = get(hObject,'value');
        
        if checkval == 1
            set(h_all.h_fL,'visible','on');
        elseif checkval == 0
            set(h_all.h_fL,'visible','off');
        end
        
        set(1337,'userdata',h_all);
    end
    function viewC3Dfun_chck7(hObject,eventdata)
        h_all = get(1337,'userdata');
        checkval = get(hObject,'value');
        
        if checkval == 1
            set(h_all.h_fR,'visible','on');
        elseif checkval == 0
            set(h_all.h_fR,'visible','off');
        end
        
        set(1337,'userdata',h_all);
    end


    % Axis limits callbacks
    function viewC3Dfun_xlim(hObject,eventdata)
        h_all = get(1337,'userdata');
        
        strval = get(hObject,'string');
        
        try
            eval(['set(h_all.h_ax,''xlim'',[' strval '])']);
        end
    end
    function viewC3Dfun_ylim(hObject,eventdata)
        h_all = get(1337,'userdata');
        
        strval = get(hObject,'string');
        
        try
            eval(['set(h_all.h_ax,''ylim'',[' strval '])']);
        end
    end
    function viewC3Dfun_zlim(hObject,eventdata)
        h_all = get(1337,'userdata');
        
        strval = get(hObject,'string');
        
        try
            eval(['set(h_all.h_ax,''zlim'',[' strval '])']);
        end
    end
    
    % View pushbutton functions
    function viewC3Dfun_viewxy(hObject,eventdata)
        h_all = get(1337,'userdata');
        set(h_all.h_ax,'view',[0 90]);
    end
    function viewC3Dfun_viewxz(hObject,eventdata)
        h_all = get(1337,'userdata');
        set(h_all.h_ax,'view',[0 0]);
    end
    function viewC3Dfun_viewyz(hObject,eventdata)
        h_all = get(1337,'userdata');
        set(h_all.h_ax,'view',[90 0]);
    end


    % Play function
    function viewC3Dfun_play(hObject,eventdata,C3Ddata)
        h_all = get(1337,'userdata');
        idx = h_all.idx;
        
        butval = get(hObject,'value');
        
        % change button color
        if butval == 1
            set(hObject,'backgroundcolor',[0 1 0]);
        else
            set(hObject,'backgroundcolor',[1 0 0]);
        end
        
        idxrange = h_all.idxrange;
        while get(hObject,'value') 

            % Update idx
            if any(idx == idxrange(1):idxrange(end)-1)
                idx = idx + 1;
            else
                idx = idxrange(1);  % Loop forever
            end

            % Update time indicator
            set(h_all.h_time,'string',['t = ' num2str( idx ./ C3Ddata.Marker.MarkerFrameRate )] );
        
            % Update gait phase indicator (if possible)
            if isfield(h_all,'h_gpm')
                set(h_all.h_gpm,'string',['GPm = ' num2str(C3Ddata.Event.GaitPhaseM(idx,end)) ]);
            end
            
            % Update data in viewer
            if (h_all.labellen(1) ~= 0) && get(h_all.h_chck1,'value')
                data_mark_x = squeeze(C3Ddata.Marker.MarkerData(idx,:,1));
                data_mark_y = squeeze(C3Ddata.Marker.MarkerData(idx,:,2));
                data_mark_z = squeeze(C3Ddata.Marker.MarkerData(idx,:,3));
                refreshdata(h_all.h_mark,'caller');
            end
            if (h_all.labellen(2) ~= 0) && get(h_all.h_chck2,'value')
                data_prob_x = squeeze(C3Ddata.Marker.ProbedData(idx,:,1));
                data_prob_y = squeeze(C3Ddata.Marker.ProbedData(idx,:,2));
                data_prob_z = squeeze(C3Ddata.Marker.ProbedData(idx,:,3));
                refreshdata(h_all.h_prob,'caller');
            end
            if (h_all.labellen(3) ~= 0) && get(h_all.h_chck3,'value')
                data_com_x = squeeze(C3Ddata.Marker.COMData(idx,:,1));
                data_com_y = squeeze(C3Ddata.Marker.COMData(idx,:,2));
                data_com_z = squeeze(C3Ddata.Marker.COMData(idx,:,3));
                refreshdata(h_all.h_com,'caller');
            end
            if (h_all.labellen(4) ~= 0) && get(h_all.h_chck4,'value')
                data_joint_x = squeeze(C3Ddata.Marker.JointData(idx,:,1));
                data_joint_y = squeeze(C3Ddata.Marker.JointData(idx,:,2));
                data_joint_z = squeeze(C3Ddata.Marker.JointData(idx,:,3));
                refreshdata(h_all.h_joint,'caller');
            end
            
            % Force data
            if isfield(h_all,'h_fT')
                if get(h_all.h_chck5,'value')
                    data_fT_x = squeeze(fvec_T(idx,1,:));
                    data_fT_y = squeeze(fvec_T(idx,2,:));
                    data_fT_z = squeeze(fvec_T(idx,3,:));
                    refreshdata(h_all.h_fT,'caller');
                end
                if get(h_all.h_chck6,'value')
                    data_fL_x = squeeze(fvec_L(idx,1,:));
                    data_fL_y = squeeze(fvec_L(idx,2,:));
                    data_fL_z = squeeze(fvec_L(idx,3,:));
                    refreshdata(h_all.h_fL,'caller');
                end
                if get(h_all.h_chck7,'value')
                    data_fR_x = squeeze(fvec_R(idx,1,:));
                    data_fR_y = squeeze(fvec_R(idx,2,:));
                    data_fR_z = squeeze(fvec_R(idx,3,:));
                    refreshdata(h_all.h_fR,'caller');
                end
            end
            
            % Stickdata
            if isfield(h_all,'h_stick')
                stick_x = squeeze(stickData(idx,:,1));
                stick_y = squeeze(stickData(idx,:,2));
                stick_z = squeeze(stickData(idx,:,3));
                refreshdata(h_all.h_stick,'caller');
            end
            
            % Highlighted marker
            idxhigh_all = h_all.idxhigh;
            
            % Clear data_high (because it's filled up in the loop)
            data_high_x = zeros(1,length(idxhigh_all)); 
            data_high_y = zeros(1,length(idxhigh_all)); 
            data_high_z = zeros(1,length(idxhigh_all)); 
            
            % Update highlighted markers
            for idxhigh = idxhigh_all
                n = idxhigh == idxhigh_all; % indexing variable
                
                idxset = find(idxhigh <= cumsum(h_all.labellen),1,'first');

                if idxset == 1
                    data_high_x(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,1));
                    data_high_y(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,2));
                    data_high_z(1,n) = squeeze(C3Ddata.Marker.MarkerData(idx,idxhigh,3));
                elseif idxset == 2
                    data_high_x(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),1));
                    data_high_y(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),2));
                    data_high_z(1,n) = squeeze(C3Ddata.Marker.ProbedData(idx,idxhigh - h_all.labellen(1),3));
                elseif idxset == 3
                    data_high_x(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1) + h_all.labellen(2) ),1));
                    data_high_y(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1) + h_all.labellen(2) ),2));
                    data_high_z(1,n) = squeeze(C3Ddata.Marker.COMData(idx,idxhigh - ( h_all.labellen(1) + h_all.labellen(2) ),3));
                elseif idxset == 4
                    data_high_x(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1) + h_all.labellen(2)+ h_all.labellen(3) ),1));
                    data_high_y(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1) + h_all.labellen(2)+ h_all.labellen(3) ),2));
                    data_high_z(1,n) = squeeze(C3Ddata.Marker.JointData(idx,idxhigh - ( h_all.labellen(1) + h_all.labellen(2)+ h_all.labellen(3) ),3));
                end
            end
            refreshdata(h_all.h_high,'caller');
    %         refreshdata([h_all.h_mark h_all.h_prob h_all.h_high],'caller'); % This doesn't work for some reasons

            % update slider bar
            set(h_all.h_bar1,'value', idx / size(C3Ddata.Marker.MarkerData,1) )   % TO DO : use something else as the size here, because Marker.Markerdata might not exist.
    
            % update index
            h_all.idx = idx;
            
            % Capture movie frame if required
            % TODO: it might be better not to write the video in the loop
            if ~isempty(h_all.h_mov)
                writeVideo(h_all.h_mov,getframe(h_all.h_ax))
            end
            
            % Set to userdata
            set(1337,'userdata',h_all);
            
            % Pause a bit
            pause(1/C3Ddata.Marker.MarkerFrameRate);
            
        end
        
        
    end


    % Function to generate a plot from highlighted markers over selected range
    function viewC3Dfun_plot(hObject,eventdata,C3Ddata)

        % Get highlighted markers
        h_all = get(1337,'userdata');
        idxhigh_all = h_all.idxhigh;

        if ~isempty(idxhigh_all);

            % Get idxrange
            idxrange = h_all.idxrange;

            % Create new figure;
            hfig = 1337 + 1;
            while ishandle(hfig)
                hfig = hfig + 1;
            end
            figure(hfig);

            % Creat time axis
            tax = (0:size(C3Ddata.Marker.MarkerData,1)-1) ./ C3Ddata.Marker.MarkerFrameRate;

            % Plot highlighted marker in new subplot
            nplot = numel(idxhigh_all); 
            n = 0; % Plot number variable
                for idxhigh = idxhigh_all
                    n = n + 1;

                idxset = find(idxhigh <= cumsum(h_all.labellen),1,'first');	

                subplot(nplot,1,n); hold on;
                    if idxset == 1
                        plot( tax , squeeze(C3Ddata.Marker.MarkerData(:,idxhigh,1)) , '-b' , 'linewidth', 2);
                        plot( tax , squeeze(C3Ddata.Marker.MarkerData(:,idxhigh,2)) , '-r' , 'linewidth', 2);
                        plot( tax , squeeze(C3Ddata.Marker.MarkerData(:,idxhigh,3)) , '-g' , 'linewidth', 2);
                        ylabel( C3Ddata.Marker.MarkerDataLabel(idxhigh) );
                    elseif idxset == 2
                        plot( tax , squeeze(C3Ddata.Marker.ProbedData(:,idxhigh - h_all.labellen(1),1)) , '-b' , 'linewidth', 2); % Correction for the set which the marker belongs to
                        plot( tax , squeeze(C3Ddata.Marker.ProbedData(:,idxhigh - h_all.labellen(1),2)) , '-r' , 'linewidth', 2);
                        plot( tax , squeeze(C3Ddata.Marker.ProbedData(:,idxhigh - h_all.labellen(1),3)) , '-g' , 'linewidth', 2);
                        ylabel( C3Ddata.Marker.ProbedDataLabel(idxhigh  - h_all.labellen(1)) );
                    elseif idxset == 3
                        plot( tax , squeeze(C3Ddata.Marker.COMData(:,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2) ),1)) , '-b' , 'linewidth', 2);
                        plot( tax , squeeze(C3Ddata.Marker.COMData(:,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2) ),2)) , '-r' , 'linewidth', 2);
                        plot( tax , squeeze(C3Ddata.Marker.COMData(:,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2) ),3)) , '-g' , 'linewidth', 2);
                        ylabel( C3Ddata.Marker.COMDataLabel(idxhigh - ( h_all.labellen(1)+ h_all.labellen(2) ) ) );
                    elseif idxset == 4
                        plot( tax , squeeze(C3Ddata.Marker.JointData(:,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2)+ h_all.labellen(3) ),1)) , '-b' , 'linewidth', 2);
                        plot( tax , squeeze(C3Ddata.Marker.JointData(:,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2)+ h_all.labellen(3) ),2)) , '-r' , 'linewidth', 2);
                        plot( tax , squeeze(C3Ddata.Marker.JointData(:,idxhigh - ( h_all.labellen(1)+ h_all.labellen(2)+ h_all.labellen(3) ),3)) , '-g' , 'linewidth', 2);
                        ylabel( C3Ddata.Marker.JointDataLabel(idxhigh - ( h_all.labellen(1)+ h_all.labellen(2)+ h_all.labellen(3) ) ) );
                    end

                end

            % Time label on last plot
            xlabel('Time [s]');

            % Link x-axes of all subplots
            linkaxes([findall(hfig,'type','axes')],'x');

            % Set x-axes
            set([findall(hfig,'type','axes')],'xlim',[idxrange(1) idxrange(end)]./ C3Ddata.Marker.MarkerFrameRate);

            % Put in a legend
            legend('x-direction','y-direction','z-direction');

        end

    end

    % Movie generation function
    function viewC3Dfun_mov(hObject,eventdata,C3Ddata)

        % Get highlighted markers
        h_all = get(1337,'userdata');
        
        butval = get(hObject,'value');
        
        % change button color
        if butval == 1
            set(hObject,'backgroundcolor',[1 0.8 0.8],'string','STOP RECORD');
        else
            set(hObject,'backgroundcolor',[0.8 1 0.8],'string','RECORD');
        end

        % Initiate movie
        if butval == 1

            % uiputfile
            [savename,savepath] = uiputfile('.avi');
            
            if savename == 0
                % Toggle back
                set(hObject,'backgroundcolor',[0.8 1 0.8],'string','RECORD','value',0);
                butval = 0;
            else
                % Instantiate movie
                h_all.h_mov = VideoWriter([savepath savename],'Motion JPEG AVI');
%                 h_all.h_mov.FrameRate = 30;  % Frames / sec
                h_all.h_mov.FrameRate = C3Ddata.Marker.MarkerFrameRate;  % Frames / sec ; not advised to use the marker frame rate here, you'll get a 100fps movie
                open(h_all.h_mov);
                
                % Disable the toolbar (because rotating can change the video frame size)
                set(h_all.h_gui,'toolbar','none');
                
            end

        end
        
        % Write to movie file
        % SEE WHILE LOOP UNDER PLAY BUTTON
        
        % Save movie
        if (butval == 0) && ~isempty(h_all.h_mov)
            close(h_all.h_mov);
            h_all.h_mov = [];
            
            % Enable the toolbar
            set(h_all.h_gui,'toolbar','figure');
        end
        
        % Set to userdata
        set(1337,'userdata',h_all);
        
    end


end