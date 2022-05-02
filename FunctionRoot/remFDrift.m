function [fdata] = remFDrift(fdata,prctL,prctR)
%% remFDrift
% Try to get rid of drift in force channels
% 
% INPUT)
% - fdata: matrix, Nx12, containing dual plate force data.
% 
% - prctL / prctR, scalar between 1 and 100 (see doc prctile).
% Used to select 'no contact instances' from the vertical forces.
% What is a suitable prctile depends your data distribution, which in turn
% depends on the gait frequency. Suggestion: something between 5 and 25.
% 
% OUTPUT)
% - fdata: matrix, Nx12, with detrended data
% 
% NOTES)
% Input is assumed an Nx12 matrix with dual plate force data, and the
% vertical GRF data being in channels 3 (left) and 9 (right).
% 
% You could opt for a non-linear fit as well by swapping some code lines

% Mark Vlutters - 29-01-2016 - University of Twente

%% Do stuff

% Dummy time axis used for lsq fits
tax(:,2) = 0:size(fdata,1)-1;
tax(:,1) = 1;
% tax(:,3) = (0:size(fdata,1)-1).^2;

% lsq lin fit
coeffL = tax(:,1:2) \ fdata(:,3);
coeffR = tax(:,1:2) \ fdata(:,9);

% Remove slope but not intercept
fdata(:,3) = fdata(:,3) - tax(:,2)*coeffL(2);
fdata(:,9) = fdata(:,9) - tax(:,2)*coeffR(2);

% We now hope that the 'lower' values in Fz are those without floor contact
idxselL = fdata(:,3) < prctile(fdata(:,3),prctL);
idxselR = fdata(:,9) < prctile(fdata(:,9),prctR);

% For every channel, make a lsq fit to the data at the selected 
% indices and subtract it, making the data (near) zero at those points
for ichan = 1:size(fdata,2)
    
    if ichan <= 6
        coeff = tax(idxselL,:) \ fdata(idxselL,ichan);
        fdata(:,ichan) = fdata(:,ichan) - (coeff(1) + coeff(2)*tax(:,2));
%         fdata(:,ichan) = fdata(:,ichan) - (coeff(1) + coeff(2)*tax(:,2) + coeff(3)*tax(:,3));
    else
        coeff = tax(idxselR,:) \ fdata(idxselR,ichan);
        fdata(:,ichan) = fdata(:,ichan) - (coeff(1) + coeff(2)*tax(:,2));
%         fdata(:,ichan) = fdata(:,ichan) - (coeff(1) + coeff(2)*tax(:,2) + coeff(3)*tax(:,3));
    end
    
end

end
