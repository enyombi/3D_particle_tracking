function [pearsonCoeff, mandersCoeff, m1, m2, k1, k2] = getColocalizationCoeffs(ROIred, ROIgreen)
% Author: Eru K.
% Date: 15-Feb-2015

% Objective: to return colocalization coeffs.
% ref: Vadim Zinchuk and Olga Zinchuk, "Quantitative Colocalization
% Analysis of Confocal Fluorescence Microscopy Images", Current Protocols
% in Cell Biology: Supplement 39, John Wiley and Sons, Inc. (2008)

diffRed = ROIred - nanmean(ROIred(:));
diffGreen = ROIgreen - nanmean(ROIgreen(:));

pearsonCoeff = nansum(diffRed(:).*diffGreen(:))/sqrt(nansum((diffRed(:).^2))*nansum((diffGreen(:).^2)));
mandersCoeff = nansum(ROIred(:).*ROIgreen(:))/sqrt(nansum(ROIred(:).^2).*nansum(ROIgreen(:).^2));
end