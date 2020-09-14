function [slope,yint,correlation] = regressA(x,y)

N = numel(x);

xy = x.*y;
xsquared = x.^2;

numerator = (N*nansum(xy(:))) - (nansum(x(:))*nansum(y(:)));
denominator = (N*nansum(xsquared(:)))-(nansum(x(:))^2);

slope = numerator/denominator;
yint = (nansum(y(:))-slope*nansum(x(:)))/N;

% error wrt best fit line:
residual = y-(yint+slope.*x);
residual2 = residual.^2;
sumResidual = nansum(residual2(:));

% error wrt mean-value:
error = nanmean(y(:))-y;
error2 = (nanmean(y(:))-y).^2;
sumError = nansum(error2(:));
correlation = (sumError - sumResidual)/sumError; %coefficient of deterimation (same as R^2 in Microsoft Excel)
% coefficient of determination = correlation coeff^2 (aka R^2)
% R^2 == 1 means good fit
% R^2 == 0 meas fit is bad
end
