function  [DcorrFinal, dimerBoolFinal] = getDimerBool(pksBest,D,w,scaleFactor,Dtrue,Dstd)
% Modified 19-Feb-2015
%   reinstated scalefactor as an input parameter
%   redefined correction factor, f, using both the manufactuer's mean
%   particle size and disperisty. (before, you incorrectly only used
%   disperisty)

% getDimerBool.m does the following:
%   randomly adjusts the sizes of all particles in a packing such that
%   all the particles do NOT overlap.  All particles that overlap are
%   identified in boolean matrix dimerBoolFinal. 
% 
%   Particle size is computed by  
%   combining 2 parameters of the mask image (or sphere spread function,
%   ssf) using a correction factor(f) that is based on the manufactuer's
%   specification for particle size (i.e., 5.06um+/-0.44um): 
%       Dcorr(np) = Dssf(np) + f*w_ssf(np) 
%       where,
%       f = [fLower, fUpper] >= 0
%       fUpper = ((5.06+0.44um*scaleFactor) - mean(Dssf(:))/mean(w_ssf(:))
%       fLower = ((5.06-0.44um*scaleFactor) - mean(Dssf(:))/mean(w_ssf(:))

% D = diameter of inner spherical region of mask image(Dssf), units: pixels
% w = thickness of outer bright edges of mask image (w_ssf), units: pixels
% scaleFactor = um/pixels

% Dtrue = manufacture's reported particle size
% Dstd = manufacture's reported stand. dev.

% ---------------------------(Start: Derivation)---------------------------
% 
% recall, peakPlacementB.m,
%   RadiusApparent = D/2 + 2w
%                  = the diameter of the spherical region in which the mask
%                  image of a particle is drawn.
% 
% RadiusApparent is a gross overestimate of the actual size of the particle
% in the mask image and by default an overestimate of the size of the real
% image of the particle that the mask image represents:
% 
%       R_real < RadiusApparent
%              < D/2 + 2w 
%       D_real < D + 4w
% 
% IMPORTANT: all units are in pixels
% 
% The diameter of the inner region of the mask image (D) is also a gross
% underestimate of particle size (D_real):
%       
%       D/2 < R_real
%       D   < D_real
% 
% so,
%   D < D_real < D + 4*w
% 
% You know that real particle size (D_real) relates to the size of the mask
% particle:
%   D_real = f(2*RadiusApparent)
%          = f(D,w)
% 
% let,
%   D_real = D + f*w
% or
%   D_real(np) = D(np) + f*w(np)
%      where, f = correction factor
% 
% Solve for f using what info. provided by manufacturer regarding the
% real particle size.
% 
% According to manufacturer,
%   mean(D_real) = (5.06um +/- 10%)/scaleFactor %units: pixels
%   i.e.,
%   0.1*mean(D_real) <= mean(D_real) <= 1.1*mean(D_real)
%   4.6200um/scaleFactor <= mean(D_real) <= 5.5000um/scaleFactor
% 
% and so,
%   mean(D_real) = mean(D) + f*mean(w)
%   where, 
%       f = [fLower,fUpper]
% 
% solving for Lowerbound of f:
%   0.8*mean(D_real) = mean(D) + fLower*mean(w)
%   4.048um/scaleFactor = mean(D) + fLower*mean(w)
% 
%   fLower = (0.8*mean(D_real) - mean(D))/mean(w)
% 
% solving for upperbound of f:
%   1.2*mean(D_real) = mean(D) + fLower*mean(w)
%   6.072um/scaleFactor = mean(D) + fUpper*mean(w)
% 
%   fUpper = (1.2*mean(D_real) - mean(D))/mean(w)
% 
% so finally,
%   D_real(np) = D(np) + f*w(np)
%   where, 
%       f = [fLower,fUpper]
% 
% Remember: D(np) < D_real(np), therefore fLower must be >=0 
% 
% ---------------------------(Stop: Derivation)---------------------------

% note, 8-oct-2015: 'Dstd' is really just wiggle-room for adjusting
% particle sizes so that no two particles overlap
% 
% Ideally, Dstd = stand dev. reported by manufacturer
% But feel free to adjust Dstd until all single particle's don't overlap
% and the particle sizes make sense!

Dreal = Dtrue/scaleFactor; %mean(D_real); units: Dtrue = um, Dreal = pixels
dev = Dstd/scaleFactor;
DrealUpper = Dreal+dev;
DrealLower = Dreal-dev;

fLower = (DrealLower - nanmean(D(:)))/nanmean(w(:));
fUpper = (DrealUpper - nanmean(D))/nanmean(w);

%f should only be a positive
if(fLower <= 0)
    fLower = 0; 
end

Dcorr = (fLower + (fUpper-fLower).*rand(size(pksBest,1),1)).*w + D; %units: pixels
% DcorrOld = Dcorr;

% Adjust Dcorr of overlapping particles by exhausting all values of
% correction factor 
[distance, pksBestOrder] = orderPeaks(pksBest);
f = linspace(fLower,fUpper,100);
dimerBool = zeros([size(pksBest,1),size(pksBest,1)],'int8');
numParticle = randperm(size(pksBest,1)); %random listing of all particles 1 to size(pksBest,1)
for j = 1:size(pksBest,1)
    np = numParticle(j);
    
    neighbors1 = find(distance(:,np) <= Dcorr(np)); %neighboring particles whose centers are located <= Dcorr(np) away
    neighbors2 = pksBestOrder(1:12,np); %max of 12 neighbors/contacts in a crystalline packing of mono-disperse hard spheres

    neighbors = union(neighbors1,neighbors2);
    neighbors(neighbors == np) = []; %remove np from list of neighbors
        
    disp(sprintf('----------------<pksBest(%d,:), Dcorr = %.04f>----------------',np,Dcorr(np)))
    for jj = 1:numel(neighbors)
        npp = neighbors(jj);
        if(distance(npp,np) >= ((Dcorr(npp)+Dcorr(np))/2))
            disp(sprintf('\t(%d) GOOD! No adjustment! (Dcorr(%d) + Dcorr(%d))/2 = %.04f <= %.04f',jj,np,npp,(Dcorr(np)+Dcorr(npp))/2,distance(npp,np)))
        end

        if(distance(npp,np) < ((Dcorr(npp)+Dcorr(np))/2)) 
            DcorrTemp = D(np) + (f*w(np));
            DcorrTempNeighbor = D(npp) + (f*w(npp));

            %recall, D < Dreal < D + 4*w
            %D is a gross underestimate of real particle size
            %Therefore, 
            %if
            %distance(npp,np) <= ((D(npp)+D(np))/2)
            %then pksBest(np,:) and pksBest(npp,:) are on the same particle 
            if(distance(npp,np) <= ((D(npp)+D(np))/2))
                dimerBool(npp,np) = 1;
                disp(sprintf('\t(%d) Multimer! (D(%d) + D(%d))/2 = %.04f >= %.04f',jj,np,npp,(D(np)+D(npp))/2,distance(npp,np)))
            end
                        
            %if pksBest(np,:) and pksBest(npp,:) haven't been identified as
            %belonging to the same particle...
            if(dimerBool ~=1) 
                %...then continue adjusting Dcorr(np) and
                %Dcorr(npp) until they don't overlap or until
                %all possible values of DcorrTemp() and DcorrTempNeighbor()
                %have been tried  
                iterNotDimer = 0;

                idxTemp = find(DcorrTemp < Dcorr(np));
                idxTempNeighbor = find(DcorrTempNeighbor < Dcorr(npp));

                iterMax = numel(idxTemp)+numel(idxTempNeighbor); 
                
                while((distance(npp,np) < ((Dcorr(npp)+Dcorr(np))/2)) && (iterNotDimer < iterMax))
                    %change Dcorr for one particle at a time...and only change
                    %Dcorr if distance(npp,np) <
                    %((Dcorr(npp)/2)+(Dcorr(np)/2))   

                    %first reduce Dcorr(np)...
                    if(~isempty(idxTemp) && (distance(npp,np) < ((Dcorr(npp)+Dcorr(np))/2)))
                        Dcorr(np) = DcorrTemp(randi(max(idxTemp),1)); %randomly choose a particle size less than the original Dcorr(npp)
                    end

                    %...then reduce Dcorr(npp) if (dist(npp,np) is still < ((Dcorr(npp)+Dcorr(np))/2)) after reducing Dcorr(np)
                    if(~isempty(idxTempNeighbor) && (distance(npp,np) < ((Dcorr(npp)+Dcorr(np))/2)))
                        Dcorr(npp) = DcorrTempNeighbor(randi(max(idxTempNeighbor),1)); %randomly choose a particle size less than the original Dcorr(npp)
                    end
                    
                    %reduce the number of particle sizes that you choose
                    %from the next iteration
                    idxTemp = find(DcorrTemp < Dcorr(np));
                    idxTempNeighbor = find(DcorrTempNeighbor < Dcorr(npp));

                    iterNotDimer = iterNotDimer + 1;
                end
                
                %after the while-loop...
                if(distance(npp,np) < ((Dcorr(npp)+Dcorr(np))/2))
                    disp(sprintf('\t(%d) Adjustment FAILED! (Dcorr(%d) + Dcorr(%d))/2 = %.04f > %.04f',jj,np,npp,(Dcorr(np)+Dcorr(npp))/2,distance(npp,np)))
                end
                
                if(distance(npp,np) > ((Dcorr(npp)+Dcorr(np))/2))
                    disp(sprintf('\t(%d) Adjustment SUCCESS! (Dcorr(%d) + Dcorr(%d))/2 = %.04f < %.04f',jj,np,npp,(Dcorr(np)+Dcorr(npp))/2,distance(npp,np)))
                end
            end            
        end
    end
    disp(sprintf('------------------------------x-------------------------------\n',np,Dcorr(np)))
end

% Remember,
%   Dcorr(np) = D(np) + f*w(np)
%       where f = random value between fLower and fUpper
% 
%   fUpper.*w(:) + D(:) = Max. estimate of particle size taking into accont the manufacturer's avg. and dispersity. Most particles overlap at this size! 
%   fLower.*w(:) + D(:) = Min. estimate of particle size taking into accont the manufacturer's avg. and dispersity. Most particles are floaters at this size! 
DcorrFinal = [fLower.*w + D, Dcorr,fUpper.*w + D]; 
dimerBoolFinal = dimerBool;
end