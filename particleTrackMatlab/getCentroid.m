function centroid = getCentroid(centers)
% objective: to compute the centorid (aka 'center of mass') of a
% dimer/dumbell, trimer, ... any multi-mer

% centroid is eqivalent to the mean/avg of x-, y- and z- coords
%   
%   <row> = nanmean(centers(:,1))
%   <col> = nanmean(centers(:,2))
%   <slice> = nanmean(centers(:,3))

  centroid(1) = nanmean(centers(:,1));
  centroid(2) = nanmean(centers(:,2));
  centroid(3) = nanmean(centers(:,3));

end

% dimer/dumbell:
% 
%       1 
%       |
%       x
%       |
%       2
% 
% centroid = (center(1,:) + center(2,:))/2 
%
%   -or in terms of linear algebra-
% 
%       1 
%     ' |
%    '  x
%   '   |
%  0... 2
% 
% x == sum of vectors 01 and 02, where 0 == reference pt./origin (0,0,0)

% trimer:
% 
% 2
% |\
% | \
% |x \
% |___\
% 1    3

% "4-mer":
% 1_____2
% |     |
% |  x  |
% |_____|
% 4     3
