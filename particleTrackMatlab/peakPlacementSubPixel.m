function [DeltaR,pos,roundOff] = peakPlacementSubPixel(centerDp,stepSize) 
% Author: Eru K.
% Data: 22-July-2014
% Objective: to return the radial displacement or distance (DeltaR) from 
% 'centerDp'

% centerDp = sub-pixel displacement of a particle's center
%       so, -1 < centerDp(1) < 1 %displacement in x/row-direction
%           -1 < centerDp(2) < 1 %displacement in y/col-direction
%           -1 < centerDp(3) < 1 %displacement in z/slice-direction
%       1 == 1 pixel
% 
% rangeColormap = max(I) - min(I), where I == original/unedited 3D image

% Note: pixel values in the original image are indices refrencing one color
% from a colormap composed of 256 shades of color. A normalization of the
% image (I)...
% 
%    Inorm = (max(I) - I)/(max(I)-min(I))
% 
% ...This means that pixels in the normalized image (Inorm) have values
% that are some fraction of the range max(I)-min(I). 
%       i.e., 0, 1/(max(I)-min(I)) , 2/(max(I)-min(I)) , ..., 1
%
% Therefore, 1/(max(I)-min(I)) is the smallest of pixel value in Inorm and
% hence the maximum resolution of the image

% stepSize = single(1/rangeColormap);

% A visual...
%                                   _______________________    
%                    .             / |                    /|
%                   .             /--|-------------------/ | 
%            (-1,-1,-1+stepSize) /:  |      :           /  |
%                               /______________________/   |
%                    (-1,-1,-1) | :  |      :          |   |   
%                               | :  |      :          |   |
%          (-1+stepSize,-1,-1)  | :  |      :          |   |
%                    .          | :  |      x----------|-- |
%                    .          | :  |______:__________|___|
%                    .          | : /       :          |   /(1,1,1)
%                               |  /--------:          |  /
%                               | /        /           | /
%                               |_________/____________|/
% 
% note: (-1,-1,-1) and (1,1,1) are just position in the cartesian
% coordinate system overlaid on the voxel.  This does NOT change! ONLY the
% location of 'centerDp' changes. DeltaR is the radial from wherever
% 'centerDp' is located inside the voxel
% 
% centerDp = sub-pixel position given by the displacement from the center
%            of a voxel at pt x (0,0,0)
%       note: -1pixel < centerDp < 1pixel
% 
% stepSize = spacing between the Cartesian coord system overlaying the
%            voxel 

spacing = -1:stepSize:1;
[row, col, slice] = ndgrid(spacing,spacing,spacing); %3.1919seconds

pos = zeros([1,3],'single');

%find the position/index in 'spacing' that locate the displacement values
%in 'centerDp'
%   i.e., 
%       pos(1) = find(spacing == centerDp(1)) 
%       pos(2) = find(spacing == centerDp(2)) 
%       pos(3) = find(spacing == centerDp(3)) 
% 
% roundoff error prevents the displacement values in 'centerDp' from ever
% EXACTLY equaling those in spacing. 

 %So, in order to find the "correctIndex" where spacing == centerDp search
 %for the smallest difference between 'spacing' and centerDp (i.e., the
 %difference closest to zero)
 %     diff = spacing - centerDp(dim);          
 %     pos(dim) = find(diff(:) == min(diff(:)));
for dim = 1:3 %row, column, and slice position
    diff = abs(spacing-centerDp(dim));
    pos(dim) = find(diff(:) == min(diff(:))); 
end
DeltaR = sqrt(((row-spacing(pos(1))).^2)+(col-spacing(pos(2))).^2+(slice-spacing(pos(3))).^2); %Important Note: the bounds of the sub-pixel voxel (-1,-1,-1) to (1,1,1) remain the SAME. Only the position of the center (i.e., 'centerDp' or more accurately spacing(pos(:)) changes
roundOff = centerDp-[spacing(pos(1)),spacing(pos(2)),spacing(pos(3))]; %note: centerDp and [row(pos(1),1,1),row(pos(2),1,1),row(pos(3),1,1)] ought to be equal but there's a roundoff error. The for-loop above circumvents this 
end