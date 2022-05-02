function [Datastr] = M7_medianfilterMarkers(Datastr,n)
% gBMPDynUI order=1;
%
% INPUT)
% - Datastr: the data structure with any number of fields
% - n: order of the median filter
% Function to filter markers with a median filter, specially useful where
% flickering of the LEDs is observed.


%% Find the fields to remove
markerdata=Datastr.Marker.MarkerData;
% Find
% figure
% plot(markerdata(:,24,2)); hold on
filt_markerdata(:,:,1)=medfilt1(markerdata(:,:,1),n);
filt_markerdata(:,:,2)=medfilt1(markerdata(:,:,2),n);
filt_markerdata(:,:,3)=medfilt1(markerdata(:,:,3),n);
% plot(filt_markerdata(:,24,2));
Datastr.Marker.MarkerData=filt_markerdata;
Datastr.Marker.MedianFilterOrder=n;



end