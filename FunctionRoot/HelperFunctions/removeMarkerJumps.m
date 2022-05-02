function [C3Ddata] = removeMarkerJumps(C3Ddata)
%% RemoveMarkerJumps
%
% Sometimes the visualeyez does really weird things with markers, 
% resulting in large jumps of the data to a certain value.
% 
% Preferably, don't use this function

	markerData = C3Ddata.Marker.MarkerData;
	
    % Detect if data is out of visibility for a long time 
    % (If more than half of the data is missing, don't bother interpolating).
    nMissing = sum( diff(markerData,1,1) == 0 , 1);
    outMaskMiss = repmat( nMissing < 0.55*size(markerData,1) , [size(markerData,1) 1 1]);

%     % Deletion based on median value
%     markerData_md = nanmedian(abs(markerData),1);
% 	outMask = abs(markerData) > outFact.*repmat(markerData_md,[size(markerData,1) 1 1]);

%     % Deletion based on constant values, plus leading sample
%     outMask = false(size(markerData));
%     foo = outMask; bar = outMask;
%     outMask(2:end,:,:) = diff(markerData,1,1) == 0;
%     foo(1:end-1,:,:) = diff(outMask,1,1) > 0;
%     outMask = outMask | foo;
    
    % Deletion based on mode value (if marker jumps to some same but non-zero value)
    outMask = (markerData - repmat(mode(markerData,1),[size(markerData,1) 1 1])) == 0;
    
	% Make outMask equal over all dimensions
% 	outMask = repmat(sum(outMask,3) > 0, [1 1 3]);

    % Merge masks
    outMask = outMask & outMaskMiss;
    
	% Zero data will be interpolated by procC3D
	markerData(outMask) = 0;

    % Store
	C3Ddata.Marker.MarkerData = markerData;

end