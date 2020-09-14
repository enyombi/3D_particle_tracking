function [minimaIdxFinal, value] = getLocalMinima(profile)

% objective: find ALL local minima of 'profile'

% algorithm:
%   By definintion, a minimum occurs between 2 maxima so locate all maxima
%   and then find the minima between them.


% step1: threshold liberally to find local maxima....
cutoff = 0.5; 

% divide profile() into 3 sections and use the max. of each section as your
% threshold value
sect = round(linspace(1,numel(profile),5));
sect(numel(sect)) = numel(profile);

overlap = 3;
maximaIdxPre = cell([(numel(sect)-1),1]);
for j = 1:(numel(sect)-1)
    idxStart = sect(j) - overlap;
    idxStop = sect(j+1);
    
    if(idxStart < 1)
        idxStart = 1;
    end
    
    if(idxStop > numel(profile))
        idxStop = numel(profile);
    end
    idx = idxStart:idxStop;
    
    maximaIdxPre{j} = idx(find(profile(idx) >= cutoff*max(profile(idx))));
end

maximaIdx = maximaIdxPre{1};
for j = 2:(numel(sect)-1)
    maximaIdx = union(maximaIdx,maximaIdxPre{j});
end

% locate and retain only those maxima whose adjacent pts are less
% than maximaIdx
maximaIdx2 = zeros(size(maximaIdx),'single');
for num = 1:numel(maximaIdx)
    numBefore = maximaIdx(num) - 1;
    numAfter = maximaIdx(num) + 1;
    if(numBefore < 1)
        numBefore = 1;
    end
    
    if(numAfter > numel(profile))
        numAfter = numel(profile);
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
    if(profile(numBefore) < profile(maximaIdx(num)))
        if(profile(numAfter) <= profile(maximaIdx(num)))
            maximaIdx2(num) = maximaIdx(num);
        end
    end
    
    if(profile(numBefore) <= profile(maximaIdx(num)))
        if(profile(numAfter) < profile(maximaIdx(num)))
            maximaIdx2(num) = maximaIdx(num);
        end
    end
end

maximaIdx2(maximaIdx2(:) == 0) = [];

% step2: Find minima between each maxima
minimaIdx = zeros([(numel(maximaIdx2)-1),1],'int8');
if(numel(maximaIdx2) > 1)
    for j = 1:(numel(maximaIdx2)-1)
        idx = (maximaIdx2(j)+1):(maximaIdx2(j+1)-1); %exclude the local maxima
        if(numel(idx) >= 1) %you need at least 3 pts (i.e., 2 maxima and minimum in between), if there are only 3pts then numel(idx) = 1
            minimaIdx(j) = min(idx(find(min(profile(idx)) == profile(idx)))); 
        end
    end
end
minimaIdx(minimaIdx(:) == 0) = [];

minimaIdxFinal = minimaIdx;
value = profile(minimaIdxFinal);
end