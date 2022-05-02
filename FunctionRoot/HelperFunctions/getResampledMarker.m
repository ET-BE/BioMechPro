function C3Ddata = getResampledMarker(C3Ddata,sync_c3d,sync_xpc)
%% getResampledMarker
% Remove samples from the marker data and the sync signal
% This might be required if the sample frequencies are not exactly equal
% (e.g. if the target runs at 1000 Hz and the NI box to collect the sync at
% the VZ runs at 1000.01 Hz)
% 
% INPUT) 
% C3Ddata : data structure
% sync_c3d : random continuous sync signal recorded on the visualeyez
% sync_xpc : random continuous sync signal recorded (and generated by) the xpc
% 
% OUTPUT)
% C3Ddata : data structure, in which the Marker field and the Analog field
% are resampled
%
% NOTES)
% Data is ONLY resampled during the duration of the xpc sync.
% Before and after that the data is left unedited.
% 
% Samples are removed as equally spaced as possible

% Mark Vlutters - University of Twente - June 2015

%% Do stuff

% #####################################################
% factor = 10; % Between marker and VZ analog
% factor = 20; % Between marker and VZ analog
factor = C3Ddata.Analog.AnalogFrameRate ./ C3Ddata.Marker.MarkerFrameRate;
% #####################################################

datdelay = finddelay(sync_xpc,sync_c3d);

nSyncSmpl = 2500;
if datdelay > 0 % VZ started first
    
    % Find where the sync starts and ends
    idxStart_c3da = find( abs(sync_c3d) > 0.1 , 1 , 'first'); % Analog
    idxEnd_c3da = find( abs(sync_c3d) > 0.1 , 1 , 'last');
    
    idxStart_c3dm = round(idxStart_c3da ./ factor); % Marker
    idxEnd_c3dm = round(idxEnd_c3da ./ factor);

    idxStart_xpc = find( abs(sync_xpc) > 0.1 , 1 , 'first'); % XPC

else % Target started first
    
    % Find where the sync starts and ends
    idxStart_c3da = find( abs(sync_c3d) > 0.1 , 1 , 'first'); % Analog
    idxEnd_c3da = find( abs(sync_c3d) > 0.1 , 1 , 'last');

    idxStart_c3dm = round(idxStart_c3da ./ factor); % Marker
    idxEnd_c3dm = round(idxEnd_c3da ./ factor);
    
    % Try to skip the sync part that you missed with the VZ
    idxStart_xpc = finddelay(sync_c3d(idxStart_c3da:idxStart_c3da+nSyncSmpl) , sync_xpc);
    
end
% TO DO : CHECK THE BELOW LINE OF CODE. DOESN'T ALWAYS SEEM TO WORK.
if idxEnd_c3da(end) == length(sync_c3d) % VZ stopped before target   WARNING: THIS ONLY WORKS if the sync signal doesn't go below 0.1 at the end.
disp(' check' )
    idxEnd_xpc = finddelay(sync_c3d(idxEnd_c3da-nSyncSmpl+1:idxEnd_c3da) , sync_xpc ) + nSyncSmpl;

else % Target stopped before VZ
    idxEnd_xpc = find( abs(sync_xpc) > 0.1 , 1 , 'last');
end


% Get data lengths
nsmpl_c3d = idxEnd_c3da - idxStart_c3da + 1;
nsmpl_xpc = idxEnd_xpc - idxStart_xpc + 1;

% Number of samples to be removed / inserted into the c3d data
nsmpl_outa = nsmpl_c3d - nsmpl_xpc;
nsmpl_outm = abs(round(nsmpl_outa./factor));

if nsmpl_outa > 0 % c3d longer than xpc
    
    nsmpl_outa = abs(nsmpl_outa);
    
    % Indices of samples to be removed
    idx_outa = round( linspace(idxStart_c3da+nsmpl_c3d/nsmpl_outa,idxEnd_c3da, nsmpl_outa) );
    idx_outm = round( linspace(idxStart_c3dm+(nsmpl_c3d/nsmpl_outa)./factor, idxEnd_c3dm, nsmpl_outm) );
    
    names_m = fieldnames(C3Ddata.Marker);
    names_a = fieldnames(C3Ddata.Analog);
    
    % Loop through field names and remove samples
    % Remove analog samples
    for iname = 1:length(names_a)
        
        if isnumeric(C3Ddata.Analog.(names_a{iname})) 
            
            datsiz = size(C3Ddata.Analog.(names_a{iname}));
            
            if (length(datsiz) == 2) && (datsiz(1) >= idx_outa(end))  % Don't operate on scalar
                C3Ddata.Analog.(names_a{iname})(idx_outa) = [];
            elseif (length(datsiz) == 3) && (datsiz(1) >= idx_outa(end))
                C3Ddata.Analog.(names_a{iname})(idx_outa,:,:) = [];
            elseif (length(datsiz) == 4) && (datsiz(1) >= idx_outa(end))
                C3Ddata.Analog.(names_a{iname})(idx_outa,:,:,:) = [];
            end
        end
    end
    
    % Remove marker based samples
    for iname = 1:length(names_m)
        if isnumeric(C3Ddata.Marker.(names_m{iname}))
            
            datsiz = size(C3Ddata.Marker.(names_m{iname}));
            
            if (length(datsiz) ==2) && (datsiz(1) >= idx_outm(end))
                C3Ddata.Marker.(names_m{iname})(idx_outm) = [];
            elseif (length(datsiz) == 3) && (datsiz(1) >= idx_outm(end))
                C3Ddata.Marker.(names_m{iname})(idx_outm,:,:) = [];
            elseif (length(datsiz) == 4) && (datsiz(1) >= idx_outm(end))
                C3Ddata.Marker.(names_m{iname})(idx_outm,:,:,:) = [];
            end
        end
    end
    
else % c3d shorter than xpc
    % TO DO : insert samples
    disp('Samples should be inserted. Currently not available')
end


%% Generate output

C3Ddata.Marker.wasRe = true;

end