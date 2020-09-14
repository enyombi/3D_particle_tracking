function maximaIdxFinal = getLocalMaxima(profileGreen)
% Author: Eru K
% Date: 22-Sept-2014, modified 27-Oct-2014 (After writing getLocalMinima.m)
% objective: to get local maxima of discrete function, in this case the
% profile in Igreen between particles

dprofileGreen = gradient(profileGreen); %derivative of profileGreen, gradient() == center difference numerical differentiation
ddprofileGreen = gradient(dprofileGreen); %2nd derivative of profileGreen, gradient() == center difference numerical differentiation 

%Find local maxima of profileGreen by doing the following...
%   step1:
%       find local minima of the ddprofileGreen (remember: 2nd
%       derivative > 0 means concave up (aka local minima) AND 2nd
%       deriv < 0 means concave down (aka local maxima))...  
%
%   step2:
%       of the pts found in step1, locate corresponding pts in
%       dprofileGreen with values close to 0
%   
%   step3: 
%       of the pts in step2, locate corresponding pts in
%       profileGreen that are close to max(profileGreen). These pts
%       are the local maxima


% step1: find concave up using 2nd deriv. (i.e., ddprofileGreen < 0) 
cutoffD2green = 0.3;
ptsD2green = find(ddprofileGreen <= cutoffD2green*min(ddprofileGreen(:))); %locates all pts. of the 2nd deriv. that are w/in 70percent of min(2nd deriv.)

localMaxIdxFinal = cell(numel(ptsD2green),1);

%step2: find maximum by identifying which of the concave up pts has a 1st
%deriv. nearly equal to zero. (i.e., dprofileGreen ~ 0)
spread = 5; %because there's noise and dprofileGreen isn't smooth, add 5 pts above and below each ptsD2green for computation of local max. 
for pNum = 1:numel(ptsD2green)   
    idxStart = ptsD2green(pNum) - spread;
    idxStop = ptsD2green(pNum) + spread;

    if(idxStart <= 0)
        idxStart = 1;
    end

    if(idxStop > length(profileGreen))
        idxStop = length(profileGreen);
    end

    %local maxima in profileGreen appear as inflection pts. in
    %dprofileGreen (i.e., dprofileGreen == 0). 'idxGreen' locates these
    %inflection pts. 
    minDistFromZero = 0.50*min(abs(max(dprofileGreen)),abs(min(dprofileGreen))); 
    idxDgreen = find(dprofileGreen(idxStart:idxStop) <= minDistFromZero & dprofileGreen(idxStart:idxStop) >= (-1*minDistFromZero));
    ptsDgreen = idxStart + idxDgreen - 1; 
        %   ptsD2green = 13
        %   idxStart = 13-5 = 8
        %   idxStop = 13+5 = 13 
        %
        %idxDgreen = find(dprofileGreen(idxStart:idxStop) <= minDistFromZero & dprofileGreen(idxStart:idxStop) >= (-1*minDistFromZero));
        %          = 6,7
        %
        %ptsDgreen = idxDgreen(6,7) = 13,14  
        %   where dprofileGreen(ptsDgreen) = 0.0280,-0.0229

    %step3: find local maximum (i.e., max(profileGreen(idxStart:idxStop))
    if(numel(ptsDgreen) > 0)
        ptsMaxGreen = zeros([numel(ptsDgreen),1],'single'); %one local max. is retained for each ptsDgreen

        for ppNum = 1:numel(ptsDgreen)
            idxStart = ptsDgreen(ppNum) - spread;
            idxStop = ptsDgreen(ppNum) + spread;

            if(idxStart <= 0)
                idxStart = 1;
            end

            if(idxStop > length(profileGreen))
                idxStop = length(profileGreen);
            end 

            idxGreen = find(profileGreen(idxStart:idxStop) == max(profileGreen(idxStart:idxStop)));
            ptsMaxGreen(ppNum) = idxStart + max(idxGreen) - 1; %local max. in profileGreen; used max(idxGreen) to ensure that there's only one max.
        end
        
        %remove redundant values of local maxima...
        [value,idx] = sort(ptsMaxGreen);
        localMaxIdx = intersect(value,min(value):max(value));
        localMaxIdxFinal{pNum} = localMaxIdx;
    end
end

maxima = localMaxIdxFinal{1};
for np = 2:numel(ptsD2green) %size(localMaxIdxFinal)
    maxima = union(maxima,localMaxIdxFinal{np});
end

%remove redundant values of local maxima...
[value,idx] = sort(maxima);
maximaIdx = intersect(value,min(value):max(value));

% finally, locate and retain only those maxima whose adjacent pts are less
% than maximaIdx
maximaIdx2 = zeros(size(maximaIdx),'single');
for num = 1:numel(maximaIdx)
    numBefore = maximaIdx(num) - 1;
    numAfter = maximaIdx(num) + 1;
    if(numBefore < 1)
        numBefore = 1;
    end
    
    if(numAfter > numel(profileGreen))
        numAfter = numel(profileGreen);
    end
    
    %   numBefore < maximaIdx(num) & numAfter < maximaIdx(num)
    %
    %               maximaIdx(num)
    %                  /    \                
    %          numBefore      numAfter        
    %
    %                   -OR-
    %
    %   numBefore == maximaIdx(num) & numAfter < maximaIdx(num)
    %
    %           numBefore--maximaIdx(num)
    %                             \           
    %                             numAfter     
    %   
    %                   -OR-
    % 
    %   numBefore == maximaIdx(num) & numAfter == maximaIdx(num)
    %
    %           numBefore--maximaIdx(num)--numAfter
    %
    if(profileGreen(numBefore) < profileGreen(maximaIdx(num)))
        if(profileGreen(numAfter) <= profileGreen(maximaIdx(num)))
            maximaIdx2(num) = maximaIdx(num);
        end
    end
    
    if(profileGreen(numBefore) <= profileGreen(maximaIdx(num)))
        if(profileGreen(numAfter) < profileGreen(maximaIdx(num)))
            maximaIdx2(num) = maximaIdx(num);
        end
    end
end

maximaIdx2(maximaIdx2(:) == 0) = [];

% Retain only those local maxima which are true maxima (i.e.,close to the
% value of max(profileGreen)) 
spread = 5;
cutoffGreen = 0.5;
maximaIdx3 = zeros(size(maximaIdx2),'single');
for num = 1:numel(maximaIdx2)
    maximaIdx3(num) = 0;
    if(profileGreen(maximaIdx2(num)) >= cutoffGreen*max(profileGreen)) 
        maximaIdx3(num) = maximaIdx2(num);
    end
end
maximaIdx3(maximaIdx3(:) == 0) = [];
maximaIdxFinal = maximaIdx3;
end