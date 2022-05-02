function [Datastr] = getOSMAinM(Datastr,subjosmasto,varargin)
%% getOSMAinM
%
% Import OpenSim muscle analysis data from a .sto file into
% the corresponding data structure.
%
% INPUT)
% Datastr : Datastr structure
% subjosmasto : string - filename of the .sto file containing the MA variables.
%
% OUTPUT)
%  Datastr : Datastr structure with added MA variables
%
% NOTES)


%% Check input

folders = split(subjosmasto,'\');
mc = folders{end-1}; %mocap folder
%% Load data
MA_vars={'Length', 'MomentArm_knee_angle_r', 'MomentArm_ankle_angle_r'};
for i=1:length(MA_vars)
    % Read file header and data header
    fid = fopen([subjosmasto '\_MuscleAnalysis_' MA_vars{i} '.sto'],'r');
    if fid == -1
        warning(['Failed to open \_MuscleAnalysis_' MA_vars{i} '.sto in'  subjosmasto '. Skipping.']);
        return;
    end
    
    tline = '';
    ntimeout = 0;
    while isempty(strfind(tline,'time')) && ntimeout < 20
        
        tline = fgetl(fid);
        
        if ~isempty(strfind(tline,'nRows='))
            nRows = str2double( tline(strfind(tline,'nRows=') + 6 : end) );
        elseif ~isempty(strfind(tline,'nColumns='))
            nCols = str2double( tline(strfind(tline,'nColumns=') + 9 : end) );
        end
        
        ntimeout = ntimeout + 1;
    end
    fclose(fid);
    
    if ntimeout == 20
        error('getOSMAinM:DataHeader','Unable to find data column header in file.');
    else
        stoData = dlmread([subjosmasto '\_MuscleAnalysis_' MA_vars{i} '.sto'],'\t',[ntimeout,0,ntimeout+nRows-1,nCols-1]);
        % NOTE: assumed here the data is directly below the column header
    end
    
    
    %% Sort data based on header
    %
    % NOTE: It is assumed the column header is delimited with tabs
    
    headercell = strsplit(tline,'\t');
    
    % Column and layer index for data matrix
    % stoidx = zeros(1,length(headercell));
        mtuDataLabel = {'tfl_r','rect_fem_r','vas_int_r','vas_lat_r','vas_med_r','bifemlh_r','bifemsh_r','semimem_r','semiten_r','per_brev_r','per_long_r','per_tert_r','tib_ant_r','lat_gas_r','med_gas_r','soleus_r'};

    % Create indices and labelcells for data storage based on header names
  stoidx = zeros(1,length(mtuDataLabel));

% Create indices and labelcells for data storage based on header names
for ihd = 1:length(headercell)
    switch headercell{ihd}
         %Quad
        case 'tfl_r'
            stoidx(1) = ihd;
        case 'rect_fem_r'
            stoidx(2) = ihd;
        case 'vas_int_r'
            stoidx(3) = ihd;
        case 'vas_lat_r'
            stoidx(4) = ihd;
        case 'vas_med_r'
            stoidx(5) = ihd;
        %Hamstrings
        case 'bifemlh_r'
            stoidx(6) = ihd;
        case 'bifemsh_r'
            stoidx(7) = ihd;
        case 'semimem_r'
            stoidx(8) = ihd;
        case 'semiten_r'
            stoidx(9) = ihd;
         %Peroneals & Tib   
        case 'per_brev_r'
            stoidx(10) = ihd;
        case 'per_long_r'
            stoidx(11) = ihd;   
        case 'per_tert_r'
            stoidx(12) = ihd;
        case 'tib_ant_r'
            stoidx(13) = ihd;                          
        %Calfs    
        case 'lat_gas_r'
            stoidx(14) = ihd;
        case 'med_gas_r'
            stoidx(15) = ihd;         
        case 'soleus_r'
            stoidx(16) = ihd;
            
        otherwise % time / other MTU
            % Do nothing
    end
end

    
    % Create corresponding labelcell
    % NOTE: you should expand this cell if you want to add other vars.
    
    
    %% Assign data
    
    % Pre-alloc new data container
    mtuData = zeros(nRows,length(mtuDataLabel));
    
    % Put the data in a new matrix
      for ihd = 1:length(stoidx)
        mtuData(:,ihd) = stoData(:,stoidx(ihd));
    end
    
    
    switch MA_vars{i}
        case 'Length'
            if strcmp(mc,'OS')
                Datastr.MTU.lengthosData= mtuData;
            elseif strcmp(mc,'IMU')
                Datastr.MTU.lengthimuData= mtuData;
            end
        case 'MomentArm_knee_angle_r'
            if strcmp(mc,'OS')
                Datastr.MTU.MAkneeosData= mtuData;
            elseif strcmp(mc,'IMU')
                Datastr.MTU.MAkneeimuData= mtuData;
            end
        case 'MomentArm_ankle_angle_r'
            if strcmp(mc,'OS')
                Datastr.MTU.MAankleosData= mtuData;
            elseif strcmp(mc,'IMU')
                Datastr.MTU.MAankleimuData= mtuData;
            end
    end
end
Datastr.MTU.MTUDataLabel=mtuDataLabel;
Datastr.MTU.MTUFrameRate=1/mean(diff(stoData(:,1))); %time column from sto file
Datastr.MTU = orderfields(Datastr.MTU);


end