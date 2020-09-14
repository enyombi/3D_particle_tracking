function r = peakPlacementA(center,ImAxes,dim)
%Author: Eru Kyeyune-Nyombi
% filename: peakPlacementA
% adopted from peakPlacementB on 5-Feb-2015

% Objective: To return the SPHERICAL region of radial values used to
% calculate the ideal mask particle image

% centers = particle positions 
% dim = image dimensions, dim(1)pixels X dim(2)pixels X dim(3)pixels

%                   ImAxes(3,2)
%                  /
%           y     /
%           _ ImAxes(1,2)
%           |   /
%           |  /
% ImAxes(2,1) /        ImAxes(2,2)
% |_________|__________|
%          /|
%         / |
%        /  |
%       /   |
%      /     - ImAxes(1,1)
%     /
% ImAxes(3,1) 
rangeR = cell([1,3]);
for jj = 1:3
    rangeR{jj} = linspace(ImAxes(jj,1),ImAxes(jj,2),dim(jj));
    if(dim(jj) <= 0)
        rangeR{jj} = 0;
    end
end

[row, col, slice] = ndgrid(rangeR{1},rangeR{2},rangeR{3}); %row == y, col == x, slice == z
r = sqrt((row - center(1)).^2 + (col - center(2)).^2 + (slice - center(3)).^2);       
end