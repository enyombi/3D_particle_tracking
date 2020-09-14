function [rows, cols, slice] = getCutoutSectIdx(center1, center2, sizImage)
%Author: Eru K.
% adopted from getLineProfile.m on 15-Feb-2015

% objective: To return the rows colums and slice of a region-of-interest
% (ROI) in an image from center1 to center2

%         ______  
%        /     /|
%      1/_____/ |     
%       |     | |
%       |     | /2
%       |_____|/     
%             
r1 = uint16(center1(1));
r2 = uint16(center2(1));

c1 = uint16(center1(2));
c2 = uint16(center2(2));

s1 = uint16(center1(3));
s2 = uint16(center2(3));

% --------------------------(correction)---------------------------
if(r1 <= 0)
    r1 = 1;
end

if(r1 > sizImage(1))
    r1 = sizImage(1);
end

if(r2 <= 0)
    r2 = 1;
end

if(r2 > sizImage(1))
    r2 = sizImage(1);
end

if(c1 <= 0)
    c1 = 1;
end

if(c1 > sizImage(2))
    c1 = sizImage(2);
end

if(c2 <= 0)
    c2 = 1;
end

if(c2 > sizImage(2))
    c2 = sizImage(2);
end

if(s1 <= 0)
    s1 = 1;
end

if(s1 > sizImage(3))
    s1 = sizImage(3);
end

if(s2 <= 0)
    s2 = 1;
end

if(s2 > sizImage(3))
    s2 = sizImage(3);
end

% uint() used for consistency with getLineProfile.m
% now must change r1, r2, c1, etc. to doubles in order to use linspace()
r1 = double(r1);
r2 = double(r2);
c1 = double(c1);
c2 = double(c2);
s1 = double(s1);
s2 = double(s2);

rows = linspace(r1,r2,abs(r2-r1)+1); %plus 1 to account for matlab's indexing
cols = linspace(c1,c2,abs(c2-c1)+1);
slice = linspace(s1,s2,abs(s2-s1)+1);
end