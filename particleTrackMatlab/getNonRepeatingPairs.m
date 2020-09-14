function pairsNonRepeating = getNonRepeatingPairs(pairs)
% Author: Eru K.
% Date: 9-Nov-2014
% objective: to return all the non-repeating pairs listed in 'pairs'

% example:
% 
% pairs = [
% 
%   24      1010
%   36      37
%   43      969
%   10      11
%   37      36
%   501     512
%   1010    24];
% 
% pairBool = 
%           (24)    (36)    (43)    (10)    (37)    (501)   (1010)
%   (1010)  0       0       0       0       0       0       -1
%   (37)    0       0       0       0       -2      0       0
%   (969)   0       0
%   (11)    0       0
%   (36)    0       2
%   (512)   0       0
%   (24)    1       0
% 
% note: max(pairBool(:)) = # of duplicate/palindrome pair eliminated
% 
% pairNonRepeating = 
% 
%   24      1010
%   36      37
%   43      969
%   10      11
%   501     512

duplicates = intersect(pairs(:,1),pairs(:,2));

pairBool = zeros([size(pairs,1), size(pairs,1)],'int8');

for numD = 1:numel(duplicates)
    idx1 = find(pairs(:,1) == duplicates(numD));
    idx2 = find(pairs(:,2) == duplicates(numD));
    
    %check if pair(idx1,:) == pair(idx2,:)
    %     i.e., 
    %       pair(idx1,:) = 24,1010
    %       pair(idx2,:) = 1010,24
    if(pairBool(idx2,idx1) == 0) % pair(idx2,idx1) == 0 means that the [pair(idx1,1), pair(idx2,1)] has been processed
        if(pairs(idx1,1) == pairs(idx2,2))
            if(pairs(idx1,2) == pairs(idx2,1))
                pairBool(idx2,idx1) = numD;
                pairBool(idx1,idx2) = -1*numD;
            end
        end
    end
end
[row,col] = ind2sub(size(pairs),find(pairBool < 0));
pairs(col,:) = [];
pairsNonRepeating = pairs;
end