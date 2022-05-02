function [C3DRaw] = readC3D(FullFileName)
% GetC3D:	Getting 3D coordinate/analog data from a C3D file
%
% Input:	FullFileName - file (including path and extension) to be read
%
% Output:
% MarkerData         3D-marker data [Nmarkers x NvideoFrames x Ndim(=3)]
% VideoFrameRate     Frames/sec
% ANALOGdat          Analog signals [Nsignals x NanalogSamples ]
% AnalogFrameRate    Samples/sec
% Event              Event(Nevents).time ..value  ..name
% ParameterGroup     ParameterGroup(Ngroups).Parameters(Nparameters).data ..etc.
% CameraInfo         MarkerRelated CameraInfo [Nmarkers x NvideoFrames]
% ResidualError      MarkerRelated ErrorInfo  [Nmarkers x NvideoFrames]

% AUTHOR(S) AND VERSION-HISTORY
% Ver. 1.0 Creation (Alan Morris, Toronto, October 1998) [originally named "getc3d.m"]
% Ver. 2.0 Revision (Jaap Harlaar, Amsterdam, april 2002)

% Modified by May Liu, Dec 2004. Added the HeaderGroup, timeVector, and
% changed some of the type parameters (e.g., int8, float32, etc)

% Modified by Michael Schwartz, Dec 2004. Changed data reading; eliminated
% loops (read in blocks with skip) - dramatically faster (150x - 500x).

% Modified by Francois Cauwe, Leuven, Aug 2007. Changed data reading: check
% if measurement is valid to ensure all markers all visible during
% registration

% Modified by Mark Vlutters, Enschede, June 2011.
% Machine format 'd' for DEC-VAX is no longer supported by Matlab 2007 and
% higher. Use freadVAXD instead. This unfortunately slows down the process
% for files with DEC-VAX. Can now read files containing 2^16 frames instead
% of 2^15 (int16->uint16 in header). AnalogData is scaled with the genScale.
% Changed output to structure format.

% ################################################
%                                               ##
LoadAnalog = 1; % Set 0 to skip analog data ; ##
%                                               ##
% ################################################

MarkerData = [];
VideoFrameRate = 0;
ANALOGdat = [];
AnalogFrameRate = 0;
Event = [];
ParameterGroup = [];
CameraInfo = [];
ResidualError = [];
HeaderGroup = [];

% ###############################################
% ##                                           ##
% ##    open the file and get general info     ##
% ##                                           ##
% ###############################################

%FullFileName = 'g:\Matlab\c3dreader\900008 21.c3d'
ind = findstr(FullFileName,'\');

if ind > 0
    FileName = FullFileName(ind(length(ind))+1:length(FullFileName));
else
    FileName=FullFileName;
end

fid = fopen(FullFileName,'r','ieee-le');                  % 'l' = little endian - works for both DEC and INTEL/PC

if fid == -1
    h = errordlg(['File: ',FileName,' could not be opened'],'application error');
    uiwait(h)
    return
end

NrecordFirstParameterblock = fread(fid,1,'int8');    % Reading record number of parameter section
key = fread(fid,1,'int8');                           % key = 80;

if key ~= 80,
    h = errordlg(['File: ',FileName,' does not comply to the C3D format'],'application error');
    uiwait(h)
    fclose(fid)
    return
end


fseek(fid,512*(NrecordFirstParameterblock-1)+3,'bof');      % jump to processortype - field
proctype = fread(fid,1,'int8')-83;                          % proctype: 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)

% if proctype==-83 % (M)
%    warning('ReadC3D_mhs:UnknownProcType','Processor Type not found in C3D file. Will attempt to read anyway.');
%    % If processor type not found and format is DEC-VAX data is read incorrect.
% end

% NOTE : THIS IS NOT SUPPORTED IN MATLAB VERSIONS 2008 AND HIGHER!
% USE freadVAXD INSTEAD !
% if proctype == 2,
%     fclose(fid);
%     fid = fopen(FullFileName,'r','d');                      % DEC VAX D floating point and VAX ordering
% end

% ###############################################
% ###############################################
% ##                                           ##
% ##             READ HEADER                   ##
% ##                                           ##
% ###############################################
% ###############################################

fseek(fid,2,'bof');                                     % set pointer just before word 2
Nmarkers = fread(fid,1,'uint16');                        % word 2: number of markers
NanalogSamplesPerVideoFrame = fread(fid,1,'uint16');     % word 3: number of analog measurements = chann x #anl frames per video frame
StartFrame = fread(fid,1,'uint16');                      % word 4: # of first video frame
EndFrame = fread(fid,1,'uint16');                        % word 5: # of last video frame
MaxInterpolationGap = fread(fid,1,'uint16');             % word 6: maximum interpolation gap allowed (in frame)
if proctype==2
    Scale = freadVAXD(fid,1,'float32');                 % word 7-8: floating-point scale factor to convert 3D-integers to ref system units
else
    Scale = fread(fid,1,'float32');
end
NrecordDataBlock = fread(fid,1,'uint16');                % word 9: starting record number for 3D point and analog data
NanalogFramesPerVideoFrame = fread(fid,1,'uint16');      % word 10: number of analog samples per 3d frame
if proctype==2
    VideoFrameRate = freadVAXD(fid,1,'float32');        % word 11-12: 3D frame rate
else
    VideoFrameRate = fread(fid,1,'float32');
end

if NanalogFramesPerVideoFrame > 0,
    NanalogChannels = NanalogSamplesPerVideoFrame/NanalogFramesPerVideoFrame;
else
    NanalogChannels = 0;
end

AnalogFrameRate = VideoFrameRate*NanalogFramesPerVideoFrame;


% ###############################################
% ###############################################
% ##                                           ##
% ##           READ EVENTS                     ##
% ##                                           ##
% ###############################################
% ###############################################

fseek(fid,298,'bof');                       % place pointer before 150th word (bytes 299 and 300)
EventIndicator = fread(fid,1,'int16');      % word 150: key value (12345 decimal) indicates 4 char event labels

if EventIndicator == 12345
    
    Nevents = fread(fid,1,'uint16');         % word 151: number of events
    fseek(fid,2,'cof');                     % skip one word (2 bytes)
    
    % read in the event times, values, and names
    if Nevents>0
        for i = 1:Nevents
            if proctype==2
                Event(i).time = freadVAXD(fid,1,'float');                   % read in event times
            else
                Event(i).time = fread(fid,1,'float');
            end
        end
        fseek(fid,188*2,'bof');
        for i = 1:Nevents
            Event(i).value = fread(fid,1,'int8');                   % read in event values: 0x00 = ON, 0x01 = OFF
        end
        fseek(fid,198*2,'bof');
        for i = 1:Nevents
            Event(i).name = cellstr(char(fread(fid,4,'char')'));    % read in event names
        end
    end
    
end

% ###############################################
% ###############################################
% ##                                           ##
% ##        READ 1st PARAMETER BLOCK           ##
% ##                                           ##
% ###############################################
% ###############################################

fseek(fid,512*(NrecordFirstParameterblock-1),'bof');

dat1 = fread(fid,1,'int8');
key2 = fread(fid,1,'int8');                     % key = 80
NparameterRecords = fread(fid,1,'int8');       % number of parameter blocks to follow
proctype = fread(fid,1,'int8')-83;              % proctype: 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)


% This is the initial read of Nchar... and Group... Subsequently, these are
% read from within the while loop (below).

Ncharacters = fread(fid,1,'int8');   	% characters in group/parameter name
GroupNumber = fread(fid,1,'int8');      % id number -ve = group / +ve = parameter

while Ncharacters > 0                   % While loop to read in parameter section...
    % The end of the parameter record is indicated by <0 characters for group/parameter name
    
    if GroupNumber < 0                  % Group data
        
        GroupNumber = abs(GroupNumber);
        GroupName = fread(fid,[1,Ncharacters],'char');
        ParameterGroup(GroupNumber).name = cellstr(char(GroupName));                    % group name
        
        filepos = ftell(fid);                           % present file position
        offset = fread(fid,1,'int16');                  % offset in bytes
        nextrec = filepos+offset;                       % position of beginning of next record
        
        deschars = fread(fid,1,'int8');                                                 % description characters
        GroupDescription = fread(fid,[1,deschars],'char');
        ParameterGroup(GroupNumber).description = cellstr(char(GroupDescription));      % group description
        
        ParameterNumberIndex(GroupNumber)=0;
        
        %fseek(fid,offset-3-deschars,'cof');        % pointer to next group
        fseek(fid,nextrec,'bof');                   % pointer to next group
        
    else                                % parameter data
        
        clear dimension;
        
        ParameterNumberIndex(GroupNumber) = ParameterNumberIndex(GroupNumber)+1;
        ParameterNumber = ParameterNumberIndex(GroupNumber);                            % index all parameters within a group
        
        ParameterName=fread(fid,[1,Ncharacters],'char');                                % name of parameter
        
        % read parameter name
        if size(ParameterName) > 0
            ParameterGroup(GroupNumber).Parameter(ParameterNumber).name...
                = cellstr(char(ParameterName));                                         %save parameter name
        end
        
        % read offset
        filepos = ftell(fid);                           % present file position
        offset = fread(fid,1,'int16');                  % offset of parameters in bytes
        nextrec = filepos + offset;                     % position of beginning of next record
        
        % read type
        type = fread(fid,1,'int8');                     % type of data: -1=char/1=byte/2=integer*2/4=real*4
        ParameterGroup(GroupNumber).Parameter(ParameterNumber).datatype = type;
        
        
        % read number of dimensions
        dimnum = fread(fid,1,'int8');
        
        if dimnum == 0
            datalength = abs(type);                                 %length of data record
        else
            mult = 1;
            for j = 1:dimnum
                dimension(j) = fread(fid,1,'uint8');                %changed to "uint8" (was previously int8)
                mult = mult*dimension(j);
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).dim(j)...
                    = dimension(j);                                 %save parameter dimension data
            end
            datalength = abs(type)*mult;							%length of data record for multi-dimensional array
        end
        
        % Read in the data
        %==================================================================
        %                       CHARACTER
        %==================================================================
        if type == -1                                               %datatype == 'char'
            
            wordlength = dimension(1);                              %length of character word
            
            if dimnum == 2 & datalength>0                           %& parameter(idnumber,index,2).dim>0
                
                for j = 1:dimension(2)
                    data = fread(fid,[1,wordlength],'char');        %character word data record for 2-D array
                    ParameterGroup(GroupNumber).Parameter(ParameterNumber).data(j)...
                        = cellstr(char(data));
                end
                
            elseif dimnum==1 & datalength>0
                
                data=fread(fid,[1,wordlength],'char');              %numerical data record of 1-D array
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data...
                    = cellstr(char(data));
            end
            
            %==================================================================
            %                       BOOLEAN
            %==================================================================
        elseif type == 1    %1-byte for boolean
            
            Nparameters = datalength/abs(type);
            data = fread(fid,Nparameters,'int8');
            ParameterGroup(GroupNumber).Parameter(ParameterNumber).data = data;
            
            %==================================================================
            %                       INTEGER
            %==================================================================
        elseif type == 2 & datalength > 0                       %integer
            
            Nparameters = datalength/abs(type);
            data = fread(fid,Nparameters,'int16');
            
            if dimnum > 1
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data...
                    = reshape(data,dimension);
            else
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data...
                    = data;
            end
            
            %==================================================================
            %                       FLOATING POINT
            %==================================================================
        elseif type == 4 & datalength > 0
            
            Nparameters = datalength/abs(type);
            if proctype==2
                data = freadVAXD(fid,Nparameters,'float');
            else
                data = fread(fid,Nparameters,'float');
            end
            if dimnum > 1
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data...
                    = reshape(data,dimension);
            else
                ParameterGroup(GroupNumber).Parameter(ParameterNumber).data...
                    = data;
            end
            
        end
        
        deschars = fread(fid,1,'int8');							%description characters
        
        if deschars > 0
            description = fread(fid,[1,deschars],'char');
            ParameterGroup(GroupNumber).Parameter(ParameterNumber).description...
                = cellstr(char(description));
        end
        
        fseek(fid,nextrec,'bof');       %moving ahead to next record
        
    end
    
    % check group/parameter characters and idnumber to see if more records present
    Ncharacters=fread(fid,1,'uint8');                % characters in next group/parameter name
    GroupNumber=fread(fid,1,'int8');				% id number -ve=group / +ve=parameter
    
    
end                     % end of while loop reading parameter section

% ###############################################
% ###############################################
% ##                                           ##
% ##                READ DATA                  ##
% ##                                           ##
% ###############################################
% ###############################################

NvideoFrames = EndFrame - StartFrame + 1;

% Mike's code for reading in the data follows
% The code assumes Scale < 0, this can be easily modified

% renamed variables, only for readability
NVF = NvideoFrames;                     % number of video frames
NAPVF = NanalogFramesPerVideoFrame;     % number of analog frames per video frame
NM = Nmarkers;                          % number of markers
NC = NanalogChannels;                   % number of analog channels

% Read in Markers as repeated blocks
fseek(fid,(NrecordDataBlock-1)*512,'bof');

if proctype==2
    % NOTE: freadVAXD does not work in blocks ;                 (M)
    % need to shift all upper and lower bits individually
    
    h=waitbar(0,[FileName ' is loading ; freadVAXD']);
    tmpMKR=zeros(NVF,NM,3);
    tmpANL=zeros(NVF*NAPVF,NC);
    for i=1:NVF
        for j=1:NM
            tmpMKR(i,j,1:3) = freadVAXD(fid,3,'float');
            a=fix(freadVAXD(fid,1,'float32'));
        end
        if LoadAnalog==0
            fseek(fid,NAPVF*NC*4,'cof');
        else
            for j=1:NAPVF
                tmpANL(j+NAPVF*(i-1),1:NC) = freadVAXD(fid,NC,'float')';
            end
        end
        waitbar(i/NVF);
    end
    close(h);
    %     MarkerData = permute(tmpMKR,[3,2,1]); % Optional. Also see line 407 ; switch on to get [dim,marker,data] instead of [data,marker,dim];
    MarkerData=tmpMKR; % Switch this off if you switch the line above on.
    ANALOGdat = tmpANL;
    
else % Read in blocks
    reps = [int2str(4*NM),'*float32'];
    tmpMKR = fread(fid,4*NM*NVF,reps,4*NAPVF*NC);   % note four bytes per analog channel data skipped
    
    % Reshape/resize/reorder data
    % First do the markers
    tmpMKR = reshape(tmpMKR,4,NM,NVF);              % go from one long column to a 4xNMxNVF matrix
    MarkerData = tmpMKR(1:3,:,:);                      % trim off the residual/camera contribution stuff%
    MarkerData = permute(MarkerData,[3,2,1]);  % Optional. Also see line 396 ; switch off to get [dim,marker,data] instead of [data,marker,dim];
    a = squeeze(fix(tmpMKR(4,:,:)));
end

%% modification Francois Cauwe, Leuven, Aug 2007: Check if the measurement is valid
%% See note for programmers: http://www.c3d.org/HTML/notesforprogrammers4.htm
%% If a point is invalid (i.e. marker not visible) then the fourth word is -1 and the XYZ are zero

check=0;
if check~=0
    checkWord4 = (~(tmpMKR(4,:,:)==intmax('uint16'))*1);
    
    %checkIsZero = sum(abs(tmpMKR(1:3,:,:)),1);
    checkNaN = checkWord4./checkWord4; % Becomes NaN when 0/0, otherwise it's 1 (1/1)
    MarkerData(1,:,:) = MarkerData(1,:,:).* checkNaN;
    MarkerData(2,:,:) = MarkerData(2,:,:).* checkNaN;
    MarkerData(3,:,:) = MarkerData(3,:,:).* checkNaN;
end


if (NC > 0)&&(proctype~=2)&&(LoadAnalog~=0)
    % Read in analog data as repeated blocks
    fseek(fid,(NrecordDataBlock-1)*512+4*4*NM,'bof');
    reps = [int2str(NAPVF*NC),'*float32'];
    
    tmpANL = fread(fid,NAPVF*NC*NVF,reps,4*4*NM);   % note: four bytes per 3d data skipped
    
    % Reshape the analog data
    ANALOGdat = reshape(tmpANL,NC,NAPVF*NVF);
    ANALOGdat = permute(ANALOGdat,[2 1]);
elseif (proctype==2)&&(LoadAnalog~=0)
    % Do nothing ; ANALOGdat loaded above;
else
    ANALOGdat = 0;
end

if ANALOGdat ~= 0;
    
    
    % Find scales
    for iGroup=1:length(ParameterGroup)
        if ~isempty(ParameterGroup(iGroup).name)       % added SS Fricke, only do this if not empty
            if strcmp(ParameterGroup(iGroup).name{1},'ANALOG')
                for iSubgroup=1:length(ParameterGroup(iGroup).Parameter)
                    if strcmp(ParameterGroup(iGroup).Parameter(iSubgroup).name{1},'OFFSET')
                        offsets=ParameterGroup(iGroup).Parameter(iSubgroup).data;
                    elseif strcmp(ParameterGroup(iGroup).Parameter(iSubgroup).name{1},'GEN_SCALE')
                        genScale=ParameterGroup(iGroup).Parameter(iSubgroup).data;
                    elseif strcmp(ParameterGroup(iGroup).Parameter(iSubgroup).name{1},'SCALE')
                        scale=ParameterGroup(iGroup).Parameter(iSubgroup).data;
                    end
                end
            end
        end
    end
    
    
    
    if proctype ~= -83
        offsets=offsets+(offsets<0)*2^16;
        ANALOGdat=((ANALOGdat+(ANALOGdat<0)*2^16)-repmat(offsets',size(ANALOGdat,1),1))*genScale.*repmat(scale',size(ANALOGdat,1),1)*20/2^16;
    else
        ANALOGdat=ANALOGdat.*genScale;
    end
    
end




% Get camera contribution/residual
highbyte=fix(a/256);
lowbyte=a-highbyte*256;
CameraInfo=highbyte;
ResidualError=lowbyte*abs(Scale);

% close the c3d file
fclose(fid);

HeaderGroup(1).name = 'nMarkers';
HeaderGroup(1).data = Nmarkers;
HeaderGroup(2).name = 'nAnalogChannels';
HeaderGroup(2).data = NanalogChannels;
HeaderGroup(3).name = 'startFrame';
HeaderGroup(3).data = StartFrame;
HeaderGroup(4).name = 'endFrame';
HeaderGroup(4).data = EndFrame;
HeaderGroup(5).name = 'videoSampleRate';
HeaderGroup(5).data = VideoFrameRate';
HeaderGroup(6).name = 'analogSampleRate';
HeaderGroup(6).data = AnalogFrameRate;
HeaderGroup(7).name = 'startRecord';
HeaderGroup(7).data = NrecordDataBlock;
HeaderGroup(8).name = 'maxInterpolationGap';
HeaderGroup(8).data = MaxInterpolationGap;

% Generate output structure
C3DRaw.Marker.MarkerData = MarkerData;
C3DRaw.Marker.MarkerDataLabel = ParameterGroup(1).Parameter(9).data; % This is generally where the labelcell is located
C3DRaw.Marker.VideoFrameRate = VideoFrameRate;
C3DRaw.Analog.AnalogData = ANALOGdat;
C3DRaw.Analog.AnalogFrameRate = AnalogFrameRate;
C3DRaw.Event = Event;
C3DRaw.ParameterGroup = ParameterGroup;
C3DRaw.CameraInfo = CameraInfo;
C3DRaw.ResidualError = ResidualError;

return
