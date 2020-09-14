function neighbors = getNeighbors(coorParticles, minSepDist)

% posIdx = linear index identifying particle location
% dim = size(imageData);
% minSepDist = min. separation distance between particles for them to 
%              considered touching

%0 == pair of particles are NOT separated by 'minSepDist' and therefore not
%neighbors
%1 == pair of particles are separated by 'minSepDist' and therefore are
%neighbors

% Geometrical illustration:
% minSepDist = ___
% 
%       *p4
%  *p1
% *p3      *p5       *p2

% boolNeighbors =
% 
%       p1      p2      p3      p4      p5
% p1    0       0       1       0       0           
% p2    0       0       0       0       0
% p3    1       0       0       0       0
% p4    0       0       0       0       0
% p5    0       0       0       0       0

% particles 1 and 3 are separated by a distance that is less than or equal
% to 'minSepDist'.  And therefore are neighbors

% but particles 2, 4, and 5 are separated by a distance greater than
% 'minSepDist'

boolNeighbors = zeros(size(coorParticles,1),size(coorParticles,1));  

for np = 1:size(coorParticles,1) 
    for npNeighbor = 1:size(coorParticles,1)
        if(np ~= npNeighbor)
            dx = coorParticles(np,2) - coorParticles(npNeighbor,2); %x-direction == columns
            dx2 = (dx)^2;

            dy = coorParticles(np,1) - coorParticles(npNeighbor,1);%y-direction == rows
            dy2 = (dy)^2;

            dz = coorParticles(np,3) - coorParticles(npNeighbor,3);%z-direction == slices
            dz2 = (dz)^2;        

            actualSepDist = sqrt(dx2+dy2+dz2); %actual separation distance between particles

            if(actualSepDist <= minSepDist)
                boolNeighbors(np,npNeighbor) = 1;
            end
        end
    end
end
neighbors = boolNeighbors;
end

