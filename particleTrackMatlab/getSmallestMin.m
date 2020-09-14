function [smallestMinIdx, boolMinRed] = getSmallestMin(contactRegion,profileRed,profileGreen)
% Author: Eru K.
% Date: 8-May-2016

% objective: to identify and return the background fluorescence in
% profileRed. Green and red images are colocalized, therefore a trough
% in profileGreen(contactRegion(1):contactRegion(2)) corresponds to
% excluded background signal in profileRed

% Always, always, always! use profileGreen to spatially locate any region
% because green fluorescence is localized around the edges of the
% particles. 
% 
% In this case, you're locating the excluded region around the
% trough in profileGreen

addpath('~/poincareProgs/particleTrackMatlab/');

contactRegionStart = contactRegion(1);
contactRegionStop = contactRegion(2);

if(isequal(contactRegionStart,contactRegionStop)) %single peak, contactRegion has vanished
    smallestMinIdx = nan;
    boolMinRed = 0; %smallest local min. is NOT found using profileRed
end

if(contactRegionStart < contactRegionStop) %multiple peaks; contactRegion has a width
    %want to measure excluded fluorescence so start by finding the deepest trough in profileGreen 
    tempIdxGreen = double(getLocalMinima(profileGreen(contactRegionStart:contactRegionStop))) + double(contactRegionStart-1);
    tempIdxRed = double(getLocalMinima(profileRed(contactRegionStart:contactRegionStop))) + double(contactRegionStart-1); %indices of local min. inside contactRegion

    if(numel(tempIdxGreen) == 0) %no trough/local-minima in profileGreen
        smallestMinIdx = nan;
        boolMinRed = 0; %defualt. Smallest local min. is NOT found using profileRed

        if(numel(tempIdxRed) > 0) %check for local minima in profileRed as a last resort
            tempValRed = profileRed(tempIdxRed);
            smallest = find(tempValRed == nanmin(tempValRed));
            smallestMinIdx = nanmin(tempIdxRed(smallest));
            boolMinRed = 1;%change defualt setting. Yes, smallest local min. found using profileRed instead of profileGreen, altough profileGreen more accurately locates the smallestMin
        end            
    end

    if(numel(tempIdxGreen) > 0) %yes, there is a trough and local-minima in profileGreen
        tempValGreen = profileGreen(tempIdxGreen);
        smallest = find(tempValGreen == nanmin(tempValGreen));
        smallestMinIdx = nanmin(tempIdxGreen(smallest));
        boolMinRed = 0; %smallest local min. is NOT found using profileRed
    end
end
end