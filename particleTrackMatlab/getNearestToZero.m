function closestToZero = getNearestToZero(boolContactRegion)
% Author: Eru K.
% Date: 07-April-2017
% 
% boolContactRegion(:,1) = position
% boolContactRegion(:,1) = boolean value (0 or 1)
% 
% Objective: to return values closest to zero
% 
% boolContactRegion =
%
%    -10     0
%     -9     0
%     -8     0
%     -7     0
%     -6     0
%     -5     1  <--the last successive value closest to zero in the negative direction
%     -4     1
%     -3     1
%     -2     1
%     -1     1
%      0     1
%      1     1
%      2     1  <--the last successive value closest to zero in the positive direction
%      3     0
%      4     0
%      5     0
%      6     0
%      7     0
%      8     0
%      9     0
%     10     0

closestToZero = nan([1,2]);

% step1: find the position of all values == 1
verified = boolContactRegion(boolContactRegion(:,2) == 1,1); %positions of all values in boolContactRegion(:,2) that equal 1

% step2: locate the last successive values closest to zero in the negative
% AND positive direction 

% step2a (for negative direction)
verifiedNeg = verified(verified <= 0);
if(isempty(verifiedNeg) == 0)
    boolSuccessiveNeg = zeros([numel(verifiedNeg),1]);

    for jj = 1:numel(verifiedNeg)
        tmpSegment = boolContactRegion((boolContactRegion(:,1) <= 0) & (boolContactRegion(:,1) >= verifiedNeg(jj)),:);
        if(tmpSegment(:,2) == 1) %every boolean value in tmpSegment == 1
            boolSuccessiveNeg(jj) = 1;
        end
    end
    closestToZero(1) = nanmin(verifiedNeg(boolSuccessiveNeg==1)); %the last successive value closest to zero in the positive direction
end        

% step2b (for positive direction)
verifiedPos = verified(verified >= 0);
if(isempty(verifiedPos) == 0)
    boolSuccessivePos = zeros([numel(verifiedPos),1]);
    
    for jj = 1:numel(verifiedPos)
        tmpSegment = boolContactRegion((boolContactRegion(:,1) >= 0) & (boolContactRegion(:,1) <= verifiedPos(jj)),:);
        if(tmpSegment(:,2) == 1) %every boolean value in tmpSegment == 1
            boolSuccessivePos(jj) = 1;
        end
    end
    closestToZero(2) = nanmax(verifiedPos(boolSuccessivePos==1)); %the last successive value closest to zero in the positive direction
end
end