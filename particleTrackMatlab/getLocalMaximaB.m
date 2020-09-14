function maximaIdxFinal = getLocalMaximaB(profileGreen)
% Author: Eru K.
% Date: 28-May-2016 (Adopted from getLocalMaxima.m)
% objective: to get local maxima of discrete function, in this case the
% profile in Igreen between particles

%Find local maxima of profileGreen by doing the following...
%   step1:
%       find local minima of the ddprofileGreen (remember: 2nd
%       derivative > 0 means concave up (aka local minima) AND 2nd
%       deriv < 0 means concave down (aka local maxima))...  
%
%   step2:
%       of the pts found in step1, locate inflection pts in
%       dprofileGreen
%   
%   step3: 
%       verify pts vetted in step2 are:
%           1) near the absolute max. in profileGreen
%           2) have a local minimum between them

dprofileGreen = vertcat(nan(1),diff(profileGreen)); %derivative of profileGreen, diff() == fwd difference numerical differentiation
ddprofileGreen = gradient(gradient(profileGreen)); %2nd derivative of profileGreen, gradient() == center difference numerical differentiation 

idxProfile = 1:numel(profileGreen);

% step1: find concave up using 2nd deriv. (i.e., ddprofileGreen < 0) 
cutoffD2green = 0.3;
ptsD2green = idxProfile(ddprofileGreen <= cutoffD2green*nanmin(ddprofileGreen));%locates all pts. of the 2nd deriv. that are w/in 70percent of min(2nd deriv.)

% step2: find local maxima. local maxima in profileGreen appear as
% inflection pts. in dprofileGreen
ptsDgreen = cell(numel(ptsD2green),1); %inflection pts. found in dprofileGreen
spread = 5; %because there's noise and dprofileGreen isn't smooth, add 5 pts above and below each ptsD2green for computation of local max. 
for pNum = 1:numel(ptsD2green)   
    idxStart = ptsD2green(pNum) - spread;
    idxStop = ptsD2green(pNum) + spread;

    idxA = idxStart:idxStop;
    idxA = idxA((idxA >=1) & (idxA<=numel(profileGreen))); %only retain indices btwn 1 and numel(profileGreen) 

    negIdx = idxA(dprofileGreen(idxA) <0);%negative values of dprofileGreen that coincide with idxA 
    inflectionIdx = negIdx - 1;%values located just before negative values in dprofileGreen
    inflectionIdx(inflectionIdx<=0) = [];%corrects for zero or negative values resulting from subtracting 1
    ptsDgreen{pNum} = inflectionIdx(dprofileGreen(inflectionIdx)>0);%By definition, a true inflection point is a positive value of dprofileGreen followed by a negatvie value at dprofileGreen(negIdx)
end

maximaIdx2 = ptsDgreen{1};
for np = 2:numel(ptsDgreen)
    maximaIdx2 = union(maximaIdx2,ptsDgreen{np});
end

% thresholding...retain only those maximaIdx2 that have values >=
% cuttoffGreen
cutoffGreen = 0.5*nanmax(profileGreen);%defualt
if(nanmax(profileGreen)>=0.55)
    cutoffGreen = 0.45*(nanmax(profileGreen)-nanmin(profileGreen)); %relative absolute max.
end
maximaIdx3 = intersect(maximaIdx2,idxProfile(profileGreen >= cutoffGreen));

% locate localMinima between maximaIdx3
localMinimaIdx = nan([numel(maximaIdx3)-1,1]);
for kk = 1:(numel(maximaIdx3)-1)
    localMinimaIdx(kk) = nanmin(idxProfile(profileGreen == nanmin(profileGreen(maximaIdx3(kk):maximaIdx3(kk+1)))));%localMinimaIdx(kk) is the local minimum between maximaIdx2(kk) and maximaIdx2(kk-1), nanmin(...) is used to return only ONE value in the event there are to local minimina, most likely adjacent pts in a trough of profileGreen
end

maximaIdx4 = nan(size(maximaIdx3));
for kk = 1:numel(maximaIdx3)
    if(profileGreen(maximaIdx3(kk)) >= 0.80*nanmax(profileGreen)) %if maximaIdx3(kk) is anywhere near the absolute maximum, then it is valid local maximum
        maximaIdx4(kk) = maximaIdx3(kk);
    end
    if(profileGreen(maximaIdx3(kk)) < 0.80*nanmax(profileGreen))%...if maximaIdx3(kk) isn't close to nanmax(profileGreen) then consider adjacent localMinima
        if(isempty(localMinimaIdx))
            maximaIdx4 = maximaIdx3;
        end
        
        if(~isempty(localMinimaIdx))
            if(kk>1)%there's localMinimaIdx(kk-1) to the left of maximaIdx3(kk) and localMinimaIdx(kk) to the right
                if(((profileGreen(maximaIdx3(kk))) >= 1.13*profileGreen(localMinimaIdx(kk-1))) && ((profileGreen(maximaIdx3(kk))) >= 1.2*profileGreen(localMinimaIdx(kk-1))))%maximaIdx3(kk) must be 20percent larger than adjacent localMinimaIdx(kk) and localMinimaIdx(kk-1) to qualify as a valid local maximum
                    maximaIdx4(kk) = maximaIdx3(kk);
                end
            end

            if(kk == 1)%there's only localMinimaIdx(kk) to the right of localMaximaIdx(kk)
                if((profileGreen(maximaIdx3(kk))) >= 1.13*profileGreen(localMinimaIdx(kk))) 
                    maximaIdx4(kk) = maximaIdx3(kk);
                end
            end

            if(kk == numel(maximaIdx3))%there's only localMinimaIdx(kk-1) to the left of localMaximaIdx(kk)
                if((profileGreen(maximaIdx3(kk))) >= 1.13*profileGreen(localMinimaIdx(kk-1)))
                    maximaIdx4(kk) = maximaIdx3(kk);
                end
            end
        end
    end
end
maximaIdx4(isnan(maximaIdx4))= [];

maximaIdxFinal = maximaIdx4;
end