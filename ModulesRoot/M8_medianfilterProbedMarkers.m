function [Datastr] = M8_medianfilterProbedMarkers(Datastr,n)
% gBMPDynUI order=1;
%
% INPUT)
% - Datastr: the data structure with any number of fields
% - n: order of the median filter
% Function to filter markers with a median filter, specially useful where
% flickering of the LEDs is observed.


%% Find the fields to remove
markerdata=Datastr.Marker.ProbedData;
% Find
figure
ch=[1:3];
subplot(2,1,1)
plot(markerdata(:,ch,2)); hold on
ch2=[9:10];
subplot(2,1,2)
plot(markerdata(:,ch2,2)); hold on

filt_markerdata(:,:,1)=medfilt1(markerdata(:,:,1),n);
filt_markerdata(:,:,2)=medfilt1(markerdata(:,:,2),n);
filt_markerdata(:,:,3)=medfilt1(markerdata(:,:,3),n);
subplot(2,1,1)
plot(filt_markerdata(:,ch,2));
subplot(2,1,2)
plot(filt_markerdata(:,ch2,2));

Datastr.Marker.ProbedData=filt_markerdata;
Datastr.Marker.MedianFilterOrder=n;



end