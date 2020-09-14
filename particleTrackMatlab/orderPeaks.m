function [distance, peakOrder] = orderPeaks(peaks)
% Author: Eru K.
% Date: 29Oct2013
% filename: orderPeaks.m

% Objective: return all the peaks/pts neighboring a single peak/pt in
% order, with the closest being first and the most distant being last


dist = zeros([size(peaks,1) size(peaks,1)]);
idx = dist;

for num_peak = 1:size(peaks,1)
    dx = abs(peaks(num_peak,1)-peaks(:,1)).^2;
    dy = abs(peaks(num_peak,2)-peaks(:,2)).^2;
    dz = abs(peaks(num_peak,3)-peaks(:,3)).^2;

    dist(:,num_peak) = sqrt(dx+dy+dz);

    [distOrder, idx(:,num_peak)] = sort(dist(:,num_peak));
end

distance = dist;
peakOrder = idx;
end

% Example to illustrate how orderPeaks.m works:
% (taken from data of hist23oct2013.m)
% 
% [dist peaksOrder] = orderPeaks(peaksXYZ);
% 
% peaksOrder(1:4,1:4) =
%      1     2     3     4
%     52    43    99    31
%     45    81    24    93
%     51    85    13    48

% Looking at the the first column, 
%   peaksXYZ(2,:) is the 52nd peak away from peaksXYZ(1,:)
%   peaksXYZ(3,:) is the 45th peak away from peaksXYZ(1,:)
%   peaksXYZ(4,:) is the 51st peak away from peaksXYZ(1,:)


    