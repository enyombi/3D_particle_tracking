function [voronoiVol, r] = peakPlacement(centersIndex,dim,D,w,dp)
%Author: Eru Kyeyune-Nyombi
% filename: peakPlacement

% Objective: To return the SPHERICAL region of radial values used to
% calculate the ideal mask particle image

% adopted from peak_placement.m. This is a trimmed down version of
% peak_placement.m that uses the actual particle size, i.e. ((D/2)+(2*w)),
% to define spherical region occupied by the particle

% centers = particle positions 
% dim = dimensions of the 3D image dim(1) = # of rows, dim(2) = # of cols,
% etc.
% D = particle diameter
% w = transition region at edges of the particle
% dp = displacement (used for locating the center at sub-pixel resolution
% by adding/subtracting <1 pixel to row, col, slice)
%   r = sqrt((row - center(1) + dp(1)).^2 + (col - center(2) + dp(2)).^2 + (slice - center(3) + dp(3)).^2);}

% [row, col, slice] = ind2sub(dim,centersIndex); %converts linear index to subscripted position(i.e. row,col,slice)
% centers(:,1) = row;
% centers(:,2) = col;
% centers(:,3) = slice;

[centers(:,1), centers(:,2), centers(:,3)] = ind2sub(dim,centersIndex); %converts linear index to subscripted position(i.e. row,col,slice)
centers = single(centers);

[row, col, slice] = ndgrid(1:dim(1),1:dim(2),1:dim(3));

if ~exist('dp','var')
    row = single(row);
    col = single(col);
    slice = single(slice);
end

if exist('dp','var')
    row = single(row)+dp(1);
    col = single(col)+dp(2);
    slice = single(slice)+dp(3);
end

% bkgrdVal = max(dim);%arbitrary value used so that when ipf([D w], r) is called the particle is distinct from the bkgrd (note: bkgrdVal should be outside the range of values used in constructing the mask particle image) 
r = max(dim)*ones(dim,'single'); %radial values for calculated image
voronoiVol = zeros(dim,'single'); %eqivalent to 'over' in Shattuck's 'pgrid.m'

for np = 1:size(centers,1)
%     rowTemp = row - centers(np,1); %y-displacement from center(np, )
%     colTemp = col - centers(np,2); %x-displacement from center(np, )
%     sliceTemp = slice - centers(np,3); %z-displacment from center(np, )

    %Eqn of Sphere:
    %     (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = Radius^2
    %     where center = (x0,y0,z0)
    rTemp = sqrt((row - centers(np,1)).^2 + (col - centers(np,2)).^2 + (slice - centers(np,3)).^2); 
             
    RadiusApparent = ((D(np)/2)+(2*w(np)));
    sphericalRegion = find(rTemp(:) <= RadiusApparent); 
    
    r(sphericalRegion) = rTemp(sphericalRegion);
    
    voronoiVol(sphericalRegion) = np;
%     disp(sprintf('In peakPlacement.m...finished locating peak %d of %d\n',np,size(centers,1)))
end
end