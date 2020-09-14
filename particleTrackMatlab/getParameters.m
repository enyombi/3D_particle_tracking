function [paramBest, RsquaredBest] = getParameters(avgProfile, Dactual, wActual, numTry)
% Author: Eru K.
% objective: use nonlinear regression on 'avgProfile' to compute parameters
% that best fit the function 'ssf'

% numTry = # of attempts made to compute the parameters 
% avgProfile = average linear profile of real image of a particle 

% initial guesses for parameters of ssf().m are:
% K = max(avgProfile(:))-min(avgProfile(:))
% D = frac(npp)*Dactual
% w = wActual
% offset = ((min1+min2)/2)
%   where,
%     min1 = minimum in the 1st half of the linear profile
%     min2 = minimum in the 2nd half of the linear profile

b = ceil((length(avgProfile)-1)/2);
r = linspace(-b,b,numel(avgProfile));
min1 = min(avgProfile(1:b));
min2 = min(avgProfile(b:numel(avgProfile)));

numerator = (-1*(numTry-ceil(numTry/2))):numTry;
frac = (numerator/numTry) + 1;

% tries nonlinear regression from ~Dactual/2 to 2*Dactual
%   e.g1., 
%       suppose numTry = 10 then (-1*(numTry-ceil(numTry/2))) == -5
%       so, 
%       (nummerator/numTry)+1 = (-5/10)+1, (-4/10)+1, ..., (0/10)+1, (1/10)+1, ..., (10/10)+1
%                             = 0.5000, 0.6000, ..., 1.0000, ..., 1.1000, 2.0000
% 
% e.g2., 
%       suppose numTry = 7 then (-1*(numTry-ceil(numTry/2))) == -3
%       so, 
%       (nummerator/numTry)+1 = (-3/7)+1, (-2/7)+1, (-1/7)+1, (0/7)+1, (1/7)+1, ..., (7/7)+1
%                             = 0.5714, 0.7143, 0.8571, 1.0000, ..., 2.0000

paramTemp = zeros([numel(frac),4],'single');
for npp = 1:numel(frac)
    try 
        paramTemp(npp,:) = nlinfit(r,avgProfile,'ssf',[max(avgProfile(:))-min(avgProfile(:)),frac(npp)*Dactual,wActual,((min1+min2)/2)]); 
%         disp(sprintf('\t attempt %d of %d using Dactaul*%.04f, SUCCESSFUL',npp,numTry,frac(npp)));
    catch err
        paramTemp(npp,:) = NaN;
%         disp(sprintf('\t attempt %d of %d using Dactual*%.04f, FAILED',npp,numTry,frac(npp)));
    end
end

nppGood = find(~isnan(paramTemp(:,1)));

if(isempty(nppGood))
    paramBest = NaN(1,4);
    RsquaredBest = NaN;
end

% compare Rsquared for all successful attempts of nonlinear regression and
% retain parameters with smallest Rsquared value
if(~isempty(nppGood))
    for num = 1:numel(nppGood)
        diff2 = (avgProfile - ssf(paramTemp(num,:),r)).^2;
        Rsquared(num) = sum(diff2(:));
    end
    idxBest = min(find(Rsquared == min(Rsquared)));
    RsquaredBest = Rsquared(idxBest);
    paramBest = paramTemp(idxBest,:);
end
end

% paramTemp(8,:)
% 
% ans =
% 
%     0.2072   71.7021    9.4753    0.6227
% 
% param(5,:)
% 
% ans =
% 
%     0.2072   71.7021    9.4753    0.6227