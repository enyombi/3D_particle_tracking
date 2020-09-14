function  [DcorrFinal, dimerBoolFinal] = getDimerBoolOLD(pksBest,D,w)
% Objective: To adjust particle size so that particles don't overlap.  And
% identify particles that do overlap as dimers/dumbells.

% Modifed 9-Nov-2014,
%   eliminated scaleFactor as in input variable
%   not to worry dimerBoolFinal was the same for 30may2014_2_29mayRLP_4x

% 5.06um = (mean(D(:)) + f*mean(w(:)))*scaleFactor %units: um
% 
% 20percent size dispersity. So,
%   4.05um = 0.8*mean(D(:))*scaleFactor
%   6.07um = 1.2*mean(D(:))*scaleFactor
% 
% assuming variability in particle size is entirely caused by 'w'... 
%   4.05um = 0.8*mean(D(:))*scaleFactor
%          = (mean(D(:)) + fLower*mean(w(:)))*scaleFactor
% 
% 0.8*mean(D(:))*scaleFactor = (mean(D(:)) + fLower*mean(w(:)))*scaleFactor
% 0.8*mean(D(:)) = mean(D(:)) + fLower*mean(w(:))
% -0.2*mean(D(:)) = fLower*mean(w(:))
% -0.2*mean(D(:))/mean(w(:)) = fLower
% 
% likewise,
% 1.2*mean(D(:))*scaleFactor = (mean(D(:)) + fUpper*mean(w(:)))*scaleFactor
% 0.2*mean(D(:)) = fUpper*mean(w(:))
% 0.2*mean(D(:))/mean(w(:)) = fUpper

% true particle size is composite of D and w that can never be smaller than
% D. So let fLower = 0
fLower = 0; %f should only be a positive adjust to D
fUpper = 0.2*mean(D(:))/mean(w(:));

Dcorr = (fLower + (fUpper-fLower).*rand(size(pksBest,1),1)).*w + D; %units: pixels

% % replace Dcorr with D for particles with unreasonably large values of w
% badP = find(w(:) > mean(w(:))*1.8);
% Dcorr(badP) = D(badP); 

% --(start: locate dimers/dumbells/trimers and correct Dcorr if you can)--
[distance, pksBestOrder] = orderPeaks(pksBest);
DcorrOld = Dcorr;

f = linspace(fLower,fUpper,100);

dimerBool = zeros([size(pksBest,1),size(pksBest,1)],'int8');
for np = 1:size(pksBest,1)
    neighbors1 = find(distance(:,np) <= Dcorr(np)); %neighboring particles whose centers are located <= Dcorr(np) away
    neighbors2 = pksBestOrder(1:12,np); %max of 12 neighbors/contacts in a crystalline packing of mono-disperse hard spheres

    neighbors = union(neighbors1,neighbors2);
    neighbors(neighbors == np) = []; %remove np from list of neighbors
        
    for npp = 1:numel(neighbors)        
        if(distance(neighbors(npp),np) < ((Dcorr(neighbors(npp))+Dcorr(np))/2)) 
            DcorrTemp = D(np) + (f*w(np));
            DcorrTempNeighbor = D(neighbors(npp)) + (f*w(neighbors(npp)));

            %particleSize = f(D,w) 
            %D is underestimate of particle size. So particles should not
            %touch using D/2 as the particle radius
            %
            %   i.e., 
            %       distance(neighbors(npp),np) > ((D(neighbors(npp))/2)+(D(np)/2))
            %
            %...but if particles are touching using only D, then this pair
            %of particle is most likely a dimer or dumbell shaped particle
            iterNotDimer = 0;
            while(distance(neighbors(npp),np) < ((Dcorr(neighbors(npp))+Dcorr(np))/2))
                %particleSize = f(D,w), espeically D.
                %...But D is still very small, so if 
                %(D(np)+D(neighbors(npp))/2) > distance(neighbors(npp),np) 
                %then these particles np and neighbors(npp) are more than
                %touching, they are different halves of the same particle!
                %(i.e., a dumbell or dimer)
                if(distance(neighbors(npp),np) <= ((D(neighbors(npp))+D(np))/2))
                    iterNotDimer = 0;
                    dimerBool(neighbors(npp),np) = 1;
                    break
                end
                
                %reduce the number of particle sizes that you choose
                %from with each iteration
                idxTemp = find(DcorrTemp < Dcorr(np));
                idxTempNeighbor = find(DcorrTempNeighbor < Dcorr(neighbors(npp)));
                
                %change Dcorr for one particle at a time...and only change
                %Dcorr if distance(neighbors(npp),np) <
                %((Dcorr(neighbors(npp))/2)+(Dcorr(np)/2))   
                %
                %'one at time' is the reason there are 2 conditions in the
                %if-statement
                if(~isempty(idxTemp) && (distance(neighbors(npp),np) < ((Dcorr(neighbors(npp))+Dcorr(np))/2)))
                    Dcorr(np) = DcorrTemp(randi(max(idxTemp),1)); %randomly choose a particle size less than the original Dcorr(neighbors(npp))
                end
                
                if(~isempty(idxTempNeighbor) && (distance(neighbors(npp),np) < ((Dcorr(neighbors(npp))+Dcorr(np))/2)))
                    Dcorr(neighbors(npp)) = DcorrTempNeighbor(randi(max(idxTempNeighbor),1)); %randomly choose a particle size less than the original Dcorr(neighbors(npp))
                end
                iterNotDimer = iterNotDimer + 1;
            end
        end
    end
end

DcorrFinal = Dcorr;
dimerBoolFinal = dimerBool;
end