% Author: Eru K.
% Objective: to measure particle size using Green images

addpath('~/poincareProgs/particleTrackMatlab/');

clearvars -except I*

% prepared a newly corrected list of particles' positions, 2mers, 3mers,
% etc. and wrote in .txt files 
% (see correctMultimers29may2014_2_28mayRefrig_4x2.m)

cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('ALLcorrected30may2014_1_30mayRCP_4x_id_row_col_slice_K_D_w_offset.txt','r');
temp = fscanf(fid,'%d %f %f %f %f %f %f %f\n',[8, Inf]);
fclose(fid);

pksBest(:,1) = temp(2,:);
pksBest(:,2) = temp(3,:);
pksBest(:,3) = temp(4,:);
K = temp(5,:)';
D = temp(6,:)';
w = temp(7,:)';
offset = temp(8,:)';


% 2mers:
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('2merALLcorrected30may2014_1_30mayRCP_4x_id_row_col_slice_K_D_w_offset.txt','r');
dimers = fscanf(fid,'%d %f %f %f %f %f %f %f\n',[8, Inf]);
fclose(fid);

npDimer = dimers(1,:)';                                                                                                                                                            

% 3mers:
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('3merALLcorrected30may2014_1_30mayRCP_4x_id_row_col_slice_K_D_w_offset.txt','r');
trimers = fscanf(fid,'%d %f %f %f %f %f %f %f\n',[8, Inf]);
fclose(fid);

npTrimer = trimers(1,:)';  

npMultimer = vertcat(npDimer,npTrimer);

npBad = find(isnan(pksBest(:,1)) == 1);

npGood = 1:size(pksBest,1);
npGood(npBad) = [];

npSingle = 1:size(pksBest);
npSingle(union(npBad,npMultimer)) = [];

% numel(npSingle)/numel(npGood)
% 
% ans =
% 
%     0.9159

% (numel(npSingle)/numel(npGood)) + (numel(npMultimer)/numel(npGood))
% 
% ans =
% 
%      1

% -------------------------(compute and save Dcorr)

numMeasurements = 10;
Dcorr = nan(size(pksBest,1),3,numMeasurements);
scaleFactor = 59.52381/1024; %um/pixels

tic;
for num = 1:numMeasurements
    [Dcorr(npGood,:,num), naught] = getDimerBool(pksBest(npGood,:),D(npGood),w(npGood),scaleFactor,5.06,0.44); %getDimerBool returns corrected particle size and the particles whose sizes could not be corrected and therefore qualify as dimers/dumbells/multi-mers
end
elapsedTime = toc; %595.9635seconds = 9.9327min for numMeasurements = 10

cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
for jj = 1:numMeasurements
    fid = fopen(sprintf('30may2014_1_30mayRCP_4x_id_DcorrLower_Dcorr_DcorrUpper_measurement%02d.txt',jj),'w');
    for np = 1:size(pksBest,1)
        fprintf(fid,'%d %.04f %.04f %.04f\n',np,Dcorr(np,1,jj),Dcorr(np,2,jj),Dcorr(np,3,jj));
    end
    fclose(fid);
end

% -------------------------(histogram, P(Dcorr))------------------------

clearvars -except I* pksBest npGood npBad npMultimer npSingle scaleFactor

% -------------------------(read Dcorr)

numMeasurements = 10; %number of files recording each measurement of Dcorr
DcorrM = nan([size(pksBest,1),numMeasurements]);
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
for num = 1:numMeasurements
    fid = fopen(sprintf('30may2014_1_30mayRCP_4x_id_DcorrLower_Dcorr_DcorrUpper_measurement%02d.txt',num));
    temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
    fclose(fid);
    order = temp(1,:)';
    DcorrM(order,num) = temp(3,:)';
end

Dcorr = nan([size(pksBest,1),1]);
for np = 1:size(pksBest,1)
    Dcorr(np) = nanmean(DcorrM(np,:));
end

cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('30may2014_1_30mayRCP_4x_id_Dcorr_CorrectionFactor.txt','r');
temp = fscanf(fid,'%d %f\n',[2,Inf]);
fclose(fid);

order = temp(1,:)';
f(order) = temp(2,:)';

for k = 1:numel(npGood)
    Dcorr(k) = Dcorr(k).*f(k);
end

% % f = ones([size(pksBest,1),1]);
% % f(npBad) = nan;
% % 
% % cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
% % fid = fopen('30may2014_1_30mayRCP_4x_id_Dcorr_CorrectionFactor.txt','w');
% % for np = 1:size(pksBest,1)
% %     fprintf(fid,'%d\t%.02f\n',np,f(np));
% % end
% % fclose(fid);

% % % view Dcorr to see if it makes sense:
% % for k = 1415:numel(npGood)
% %     np = npGood(k);
% %     [v,center] = getCutoutSectVertices(pksBest(np,:),round(Dcorr(np)*0.8),size(IeditedRed));
% %     cutoutSectGreen = IeditedGreen(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
% %     cutoutSectRed = IeditedRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
% %     
% %     figure(1) 
% %     hold off
% %     simage(cutoutSectGreen(:,:,round(center(3))));
% %     hold on
% %     diameter = [1.5:-0.1:0.5]'*Dcorr(np);
% %     plot(center(2),center(1),'k*','markersize',20)
% %     plot(center(2)-(diameter./2),center(1),'y*','markersize',20)
% %     plot(center(2)+(diameter./2),center(1),'y*','markersize',20)
% %     plot(center(2),center(1)-(diameter./2),'y*','markersize',20)
% %     plot(center(2),center(1) + (diameter./2),'y*','markersize',20)
% %         
% %     plot(center(2)-(Dcorr(np)/2),center(1),'mo','markersize',20)
% %     plot(center(2)+(Dcorr(np)/2),center(1),'mo','markersize',20)
% %     plot(center(2),center(1)-(Dcorr(np)/2),'mo','markersize',20)
% %     plot(center(2),center(1) + (Dcorr(np)/2),'mo','markersize',20)
% %     
% % 
% %     title(sprintf('pksBest(%d,:)',np),'fontsize',20)
% %     xlabel(sprintf('(y*) Dcorr*%.02f to Dcorr*%.02f \n(mo)1.00*Dcorr = %.02fum',diameter(1)/Dcorr(np),diameter(numel(diameter))/Dcorr(np),Dcorr(np)*scaleFactor),'fontsize',20)
% %     
% %     pause
% % end



DScaled = Dcorr.*scaleFactor;

% -----------------------x
% % P(Dscaled)
DrangeNum = 100; %number of particle size from nanmin(DScaled(:)) to nanmax(DScaled(:))
diffDrange = 1;
numTry = 1;
while((numTry < 100) && (diffDrange ~= 0))
    [Dfrequency, Drange] = hist(DScaled(npSingle),DrangeNum);
    DrangeNum = numel(Dfrequency(Dfrequency > 0));
    diffDrange = numel(Drange) - DrangeNum; %check for convergence 
    numTry = numTry + 1;
end

Dfrequency(Dfrequency == 0) = nan; %just in case, the while loop didn't converge

figure(1)
hold off
plot(Drange,Dfrequency,'*','markersize',20)
title(sprintf('P(Diamter) \n 30may2014-1-30mayRCP-4x \n %d single particles(%d total)',numel(npSingle),numel(npGood)),'fontsize',20)
xlabel(sprintf('Diameter \n mean = %.03f +/- %.03fum \n manaufacturer: 5.06 +/- 0.44um',nanmean(DScaled(npSingle)),nanstd(DScaled(npSingle),1)),'fontsize',20)
ylabel('Frequency','fontsize',20)
axis([0,nanmax(DScaled(npSingle)),0,nanmax(Dfrequency(:))])
set(gca,'fontsize',20)

% save datapts need to plot P(Dcorr)
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen(sprintf('30may2014_1_30mayRCP_4x_Drange_Dfreq.txt'),'w');
for jj = 1:numel(Drange)
    fprintf(fid,'%.04f %.04f\n',Drange(jj),Dfrequency(jj));
end
fclose(fid);

% -------------------------(packingFraction)------------------------------

clearvars -except I* Dcorr pksBest npGood npBad npMultimer npSingle scaleFactor

% -------------------STEP1
volSphere = (4/3)*pi*((Dcorr./2).^3);

volParticle = volSphere; %assumes all particles are spheres

% CORRECT volParticle(np) for all particles that are only partly
% inside the image (aka 'field of view')

% boundaries/dimensions of the 'field of view' containing all particles:
dimC(1,1) = nanmax(pksBest(npGood,1)); %x-coord, column
dimC(1,2) = nanmax(pksBest(npGood,2)); %y-coord, row
dimC(1,3) = nanmax(pksBest(npGood,3)); %z-coord, slice

dimC(2,1) = nanmin(pksBest(npGood,1));
dimC(2,2) = nanmin(pksBest(npGood,2));
dimC(2,3) = nanmin(pksBest(npGood,3));

% particles may have:
%   1 spherical cap if they only stick out of ONE side/face of the
%   container 
%   2 spherical caps if they are at the intersection of 2 edges
%   3 spherical caps if they are in a corner 

%      (+)slice
%     /
%    /----(+)col
%    |
%    |
%   (+)row
%               edge1L
%             ____|_____
%            /|   |edge3U 
%           / |   v   / |
% dimC(2,:)/__|______/  |
%          |  |edge3L| <----edge2U 
%  edge2L--|->|______|__|dimC(1,:)
%          |  /      |  /
%          | /    A  | /
%          |/_____|__|/
%                 |
%               edge1U

edgeBool = zeros([size(dimC),size(pksBest,1)],'uint8'); % edgeBool(:,:,np) = [1U, 2U, 3U; 1L, 2L, 3L]

% edgeBool(:,:,np) indicates which faces (aka edges or boundaries) of the
% container pksBest(np,:) is closest to 
%   e.g1, edgeBoo1(:,:,np) = [1,0,0;0,0,0] means pksBest(np,:) sticks out
%   of face 1U (i.e., pksBest(np,:) is close to just one face)
% 
%   e.g2, edgeBoo1(:,:,np) = [1,1,0;0,0,0] means pksBest(np,:) sticks out
%   of face 1U and 2U (i.e., pksBest(np,:) is at the meeting of 2 faces)
% 
%   e.g3, edgeBoo1(:,:,np) = [1,1,0;0,0,1] means pksBest(np,:) sticks out
%   of face 1U, 2U, and 3L. (i.e., pksBest(np,:) is at a corner of the
%   container)

edgeParticles = nan(1);
for dim = 1:3
    npU = find((pksBest(:,dim) + Dcorr(:)/2) > dimC(1,dim)); %particles whose edges stick out of the edge of the container/box at dimC(2,dim)
    npL = find((pksBest(:,dim) - Dcorr(:)/2) < dimC(2,dim)); %particles whose edges stick out of the edge of the container/box at dimC(2,dim)
    
    edgeBool(1,dim,npU) = 1;
    edgeBool(2,dim,npL) = 1;
    
    edgeParticles = union(edgeParticles,[npU;npL]);
end
% numel(edgeParticles)/numel(npGood) = 0.3006
% ~30% of all particles protude outside the boundaries (i.e., dimC) the
% container.  This makes sense because the faces have greater surface area

% Determine how many sides of the container a particle protrudes from
threeSides = nan([size(pksBest,1),1]);
twoSides = threeSides;
oneSide = threeSides;

% for np = 1:size(pksBest,1)
%     temp(np) = numel(find(edgeBool(:,:,np) == 1));
% end
% 
% % Checked and indeed none of the particles are intersected by more than 3
% % edges 
% % nanmax(temp(:))
% % 
% % ans =
% % 
% %      3

for np = 1:size(pksBest,1)
    if(numel(find(edgeBool(:,:,np) == 1)) == 3)
        threeSides(np) = np;
    end
    
    if(numel(find(edgeBool(:,:,np) == 1)) == 2)
        twoSides(np) = np;
    end    

    if(numel(find(edgeBool(:,:,np) == 1)) == 1)
        oneSide(np) = np;
    end
end

threeSides(find(isnan(threeSides) == 1)) = [];
twoSides(find(isnan(twoSides) == 1)) = [];
oneSide(find(isnan(oneSide) == 1)) = [];
% whos *Side*
%   Name              Size            Bytes  Class     Attributes
% 
%   oneSide         390x1              3120  double              
%   threeSides        2x1                16  double              
%   twoSides         33x1               264  double     

volCap = nan([size(pksBest,1),3]); %potential caps in all three spatial directions (i.e., row, column, slice)

% =======oneSide 
% (Correct volParticle(np) regardless if pksBest(np,:) is INSIDE or OUTSIDE
% bounds of container defined by dimC)

% particles with one spherical cap (i.e., particles protruding from only
% one side of the container)
for jj = 1:numel(oneSide)
    np = oneSide(jj);
    
    [rowTemp,colTemp] = ind2sub(size(dimC),find(edgeBool(:,:,np) == 1));
    
    if(rowTemp == 2) %dimC(2,:) is the LOWER boundary. Therefore, boolLowerBound = 1
        boolLowerBound = 1;
    end
    
    if(rowTemp == 1) %dimC(1,:) is the UPPER boundary. Therefore, boolLowerBound = 0
        boolLowerBound = 0;
    end
    
    [radiusCap,heightCap,x] = getOverlapWallSphere(pksBest(np,:),Dcorr(np)/2,dimC(rowTemp,:),colTemp,boolLowerBound); %boolLowerBound == 1

    if(x >= 0)%pksBest(np,:) is INSIDE container/box. Therefore, volParticle(np) = volSphere(np) - volCap(np)
        volCap(np,colTemp) = ((pi*heightCap)/6)*(3*(radiusCap^2) + (heightCap^2));
        volParticle(np) = volSphere(np) - volCap(np,colTemp);
    end

    if(x < 0)%pksBest(np,:) is OUTSIDE container/box. Therefore, volParticle(np) = volCap(np)
        volCap(np,colTemp) = ((pi*heightCap)/6)*(3*(radiusCap^2) + (heightCap^2));
        volParticle(np) = volCap(np,colTemp);
    end
end

% ========twoSides
% (ONLY properly corrects volParticle(np) if pksBest(np,:) is INSIDE
% container whose bouardies are defined by dimC) 

volCap(twoSides,:) = nan([numel(twoSides),3]);

% particles INSIDE the image with only two spherical caps outside the image
for j = 1:numel(twoSides)
    np = twoSides(j);
    
    [rowTemp,colTemp] = ind2sub(size(dimC),find(edgeBool(:,:,np) == 1));
    
    for jj = 1:numel(rowTemp)
        if(rowTemp(jj) == 2) %dimC(2,:) is the LOWER boundary. Therefore, boolLowerBound = 1
            boolLowerBound = 1;
        end

        if(rowTemp(jj) == 1) %dimC(1,:) is the UPPER boundary. Therefore, boolLowerBound = 0
            boolLowerBound = 0;
        end

        [rB(jj), h(jj), x(jj)] = getOverlapWallSphere(pksBest(np,:),Dcorr(np)/2,dimC(rowTemp(jj),:),colTemp(jj),boolLowerBound);
        volCap(np,colTemp(jj)) = ((pi*h(jj))/6)*(3*(rB(jj)^2) + (h(jj)^2));
    end

    %adjust size of smallest cap by exluding portion of cap outside the
    %image that intersects other caps:
    dist1 = dimC(rowTemp(1),colTemp(1)) - pksBest(np,colTemp(1));
    dist2 = dimC(rowTemp(2),colTemp(2)) - pksBest(np,colTemp(2));
    distCornerToParticle = sqrt((dist1^2) + (dist2^2));
    if((Dcorr(np)/2) > distCornerToParticle) %Checks to see that pksBest(np,:) wraps around a corner so that its caps intersect
        numSmallestCap = min(find(h == nanmin(h)));

        rowOrtho = rowTemp;
        colOrtho = colTemp;

        rowOrtho(numSmallestCap) = [];
        colOrtho(numSmallestCap) = [];%use the dimension(s) orthogonal to numSmallestcap, which is the dimension of the other cap(s)

        xDist(1) = pksBest(np,colOrtho) - rB(numSmallestCap); %distance to the left of center
        xDist(2) = pksBest(np,colOrtho) + rB(numSmallestCap); %distance to the right of center

        amtInside = abs(dimC(rowOrtho,colOrtho) - xDist(xDist >= dimC(2,colOrtho) & xDist <= dimC(1,colOrtho)));
        %when 2 caps intersect, a single segment of the circular base of
        %one cap is cut off so that the base has one straight edge like a
        %semi-circle: 
        %           
        %       . .
        %    .      .
        %   .        |
        %   .    x   |
        %   .        | 
        %    .      .
        %       . .
        %   |<------>| |
        %   | xDist    |
        %   |          |
        %   |<-------->|
        %       2*r = Diameter
        %
        %therefore, used a 1-dimensional/linear ratio between the
        %particle's diameter and length to the edge to approximate the
        %fraction of the particle inside the image (i.e., amtInside)
        %
        %   fracInside = xDist/(2*r)
        %   volCapNew = fracInside*volCap
        volCap(np,colTemp(numSmallestCap)) = (amtInside/(2*rB(numSmallestCap)))*((pi*h(numSmallestCap))/6)*(3*(rB(numSmallestCap)^2) + (h(numSmallestCap)^2));
    end
end

% apply correction to particle size
for jj = 1:numel(twoSides)
    np = twoSides(jj);
    volParticle(np) = volSphere(np) - nansum(volCap(np,:));
end


%============threeSides
% (ONLY properly corrects volParticle(np) if pksBest(np,:) is INSIDE
% container whose bouardies are defined by dimC)

volCap(threeSides,:) = nan([numel(threeSides),3]);

% particles INSIDE the image with 3 spherical caps outside the image
for j = 1:numel(threeSides)
    np = threeSides(j);

    [rowTemp,colTemp] = ind2sub(size(dimC),find(edgeBool(:,:,np) == 1));

    for jj = 1:numel(rowTemp)
        if(rowTemp(jj) == 2) %dimC(2,:) is the LOWER boundary. Therefore, boolLowerBound = 1
            boolLowerBound = 1;
        end

        if(rowTemp(jj) == 1) %dimC(1,:) is the UPPER boundary. Therefore, boolLowerBound = 0
            boolLowerBound = 0;
        end

        [rB(jj), h(jj), x(jj)] = getOverlapWallSphere(pksBest(np,:),Dcorr(np)/2,dimC(rowTemp(jj),:),colTemp(jj),boolLowerBound);
        volCap(np,colTemp(jj)) = ((pi*h(jj))/6)*(3*(rB(jj)^2) + (h(jj)^2));
    end

    %adjust size of smallest cap by exluding portion of cap outside the
    %image that intersects other caps:
    dist1 = dimC(rowTemp(1),colTemp(1)) - pksBest(np,colTemp(1));
    dist2 = dimC(rowTemp(2),colTemp(2)) - pksBest(np,colTemp(2));
    dist3 = dimC(rowTemp(3),colTemp(3)) - pksBest(np,colTemp(3));

    distCornerToParticle = sqrt((dist1^2) + (dist2^2) + (dist3^2));
    if((Dcorr(np)/2) > distCornerToParticle) %Checks to see that pksBest(np,:) wraps around a corner so that all three caps intersect
        numLargestCap = max(find(h == nanmax(h))); %exclude overlap with the largest cap from the 2 other caps

        numSmallCap = 1:3;
        numSmallCap(numSmallCap == numLargestCap) = [];%two smaller caps

        for kk = 1:numel(numSmallCap)
            rowOrtho = rowTemp;
            colOrtho = colTemp;

            rowOrtho(numSmallCap(kk)) = []; %closest boundary (i.e., if LOWER boundary is closest then dimC(rowOrtho,:) = [0,0,0] and if UPPER boundary is closest then dimC(rowOrtho,:) = [4,4,4])
            colOrtho(numSmallCap(kk)) = [];%use the dimensions orthogonal to numSmallcap, which is the dimension of the other caps

            %overlap with one cap in dimension/direction (row, column, or
            %slice) given by colOrtho(1)
            xDist(1) = pksBest(np,colOrtho(1)) - rB(numSmallCap(kk)); %distance to the left of center
            xDist(2) = pksBest(np,colOrtho(1)) + rB(numSmallCap(kk)); %distance to the right of center

            amtInsideX = abs(dimC(rowOrtho(1),colOrtho(1)) - xDist(xDist >= dimC(2,colOrtho(1)) & xDist <= dimC(1,colOrtho(1))));    

            %overlap with 2nd cap in the dimension/direction (row, column, or
            %slice) given by colOrtho(2)
            yDist(1) = pksBest(np,colOrtho(2)) - rB(numSmallCap(kk)); %distance to the left of center
            yDist(2) = pksBest(np,colOrtho(2)) + rB(numSmallCap(kk)); %distance to the right of center

            amtInsideY = abs(dimC(rowOrtho(2),colOrtho(2)) - yDist(yDist >= dimC(2,colOrtho(2)) & yDist <= dimC(1,colOrtho(2))));
            
            %Now when 3 caps intersect, 2 segment of the
            %circular base are cutoff so that it has 2 straight edges 
            %       
            %    ________   ___
            %   .        |   A
            %   .    x   |   |
            %   .        |   |yDist
            %    .      .    |
            %       . .     _V_
            %   |<------>| |
            %   | xDist    |
            %   |          |
            %   |<-------->|
            %       2*r = Diameter
            %
            %Therefore, use a 2dimensional ratio between the total base area
            %and the area of the square whose length/width is xDist/yDist
            %   note: the area of the square xDist*yDist is larger than the
            %   area of the circle pi*r^2.  You're not sure what a ratio of
            %   these 2 areas says because the shapes are so different!
            %   Instead, replace the circle with a square whose sides are
            %   length (2*r).  This square will always be larger than
            %   xDist*yDist.  The ratio of the 2 squares is the fraction of
            %   the particle inside the image
            %
            %   fracInside = (xDist*yDist)/((2*r)^2)
            %   volCapNew = fracInside*volCap
            %fracInside = (amtInsideX*amtInsideY)/(Dcorr(np)^2);
            volCap(np,colTemp(numSmallCap(kk))) = ((amtInsideX*amtInsideY)/(Dcorr(np)^2))*((pi*h(numSmallCap(kk)))/6)*(3*(rB(numSmallCap(kk))^2) + (h(numSmallCap(kk))^2));
        end
    end
end

% apply correction to particle size
for jj = 1:numel(threeSides)
    np = threeSides(jj);
    volParticle(np) = volSphere(np) - nansum(volCap(np,:));
end


% ---------------STEP3-get voronoiVol (computed in terminal using voro++)

clearvars -except I* volSphere volParticle Dcorr pksBest npGood npBad npMultimer npSingle scaleFactor

% ---------------STEP3A-print x,y,z,radius

% text file listing ALL particles except those whose centers are outside
% the image (i.e., particle listed in 'npReallyOutside')
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('pks30may2014_1_30mayRCP_4x','w');
for j = 1:numel(npGood)
    np = npGood(j);
    fprintf(fid,'%d %.04f %.04f %.04f %.04f\n',np,pksBest(np,2),pksBest(np,1),pksBest(np,3),Dcorr(np)/2);%recall, row == y-coord and column == x-coord
end
fclose(fid);

% ---------------STEP3B-print size of rectangular space with ALL particles

dimC(1,1) = nanmax(pksBest(npGood,2)); %x-coord, column
dimC(1,2) = nanmax(pksBest(npGood,1)); %y-coord, row
dimC(1,3) = nanmax(pksBest(npGood,3)); %z-coord, slice

dimC(2,1) = nanmin(pksBest(npGood,2));
dimC(2,2) = nanmin(pksBest(npGood,1));
dimC(2,3) = nanmin(pksBest(npGood,3));

dimContainer = abs(dimC(1,:) - dimC(2,:)); 

volBox = dimContainer(1)*dimContainer(2)*dimContainer(3);

cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('lbox_30may2014_1_30mayRCP_4x','w');
for j = 1:size(dimC,1)
    fprintf(fid,'%.04f\t%.04f\t%.04f\n',dimC(j,1),dimC(j,2),dimC(j,3));
end
fclose(fid);

% ---------------STEP3C-compute voronoiVol in terminal using voro++

% $ ./generateResults.exe pks29may2014_2_28mayRefrig_4x2 lbox_29may2014_2_28mayRefrig_4x2 

% ---------------STEP3D-read voronoiVol

voronoiVol = nan([size(pksBest,1),1]);

cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('pks30may2014_1_30mayRCP_4x_voro.txt','r');
B = fscanf(fid,'%d %f %f\n',[3, Inf]);
fclose(fid);

order = B(1,:)';
% radii(order) = B(2,:)';
voronoiVol(order) = B(3,:)';

% ---------------STEP3E-compute Local packingFrac
packingFrac = nan([size(pksBest,1),1]);

% packingFracMean = volParticle./voronoiVol;
for jj = 1:numel(npGood)
    np = npGood(jj);
    packingFrac(np) = volParticle(np)/voronoiVol(np);
end

% -----------------------x
% % P(packingFrac)
packingFracTemp = packingFrac;
packingFracTemp(packingFracTemp >= 1) = nan;

rangeNum = 100; %number of particle size from nanmin(DScaled(:)) to nanmax(DScaled(:))
diffRange = 1;
numTry = 1;

% % there's one particle with packingFrac = 3.2 that causes packingFracRange 
% % to converge to only 2 values. Therefore, abandon while loop:
% while((numTry < 100) && (diffRange ~= 0))
%     [packingFracFreq, packingFracRange] = hist(packingFrac(npSingle),rangeNum);
%     rangeNum = numel(packingFracFreq(packingFracFreq > 0));
%     diffRange = numel(packingFracRange) - rangeNum; %check for convergence 
%     numTry = numTry + 1;
% end

[packingFracFreq, packingFracRange] = hist(packingFrac(npSingle),rangeNum);
packingFracFreq(packingFracFreq == 0) = nan; %just in case, the while loop didn't converge

figure(3)
hold off
plot(packingFracRange,packingFracFreq,'*','markersize',20)
title(sprintf('P(packingFracLOCAL) \n 30may2014-1-30mayRCP-4x \n %d single particles (%d total)',numel(npSingle),numel(npGood)),'fontsize',20)
xlabel(sprintf('packingFracLocal \n mean = %.03f +/- %.03f \n mean(packingFrac < 1) = %.03f +/- %.03f',nanmean(packingFrac(npSingle)),nanstd(packingFrac(npSingle),1),nanmean(packingFracTemp(npSingle)),nanstd(packingFracTemp(npSingle),1)),'fontsize',20)
ylabel('Frequency','fontsize',20)
axis([0,nanmax(packingFrac(npSingle)),0,nanmax(packingFracFreq(:))])
set(gca,'fontsize',20)

% save datapts need to plot P(Dcorr)
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen(sprintf('npSingle30may2014_1_30mayRCP_4x_packingFracRange_packingFracFreq.txt'),'w');
for jj = 1:numel(packingFracRange)
    fprintf(fid,'%.04f %.04f\n',packingFracRange(jj),packingFracFreq(jj));
end
fclose(fid);


% nanmean(packingFrac(:))
% 
% ans =
% 
%     0.5566
% 
% nanmean(packingFrac(npSingle))
% 
% ans =
% 
%     0.5538
% 
% nanmean(packingFrac(packingFrac<1))
% 
% ans =
% 
%     0.5566
% 
% nanmean(packingFracTemp(npSingle))
% 
% ans =
% 
%     0.5538
% 
% nansum(volParticle(:))/volBox
% 
% ans =
% 
%     0.5530


% ------------------------(print all results)-----------------------
% (after run analysisCorrection...m and contacts_...m

clearvars -except I* packingFrac Dcorr pksBest npGood npBad npMultimer npSingle scaleFactor

% read 'boolContactMech':
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
boolContactMech = dlmread('30may2014_1_30mayRCP_4x_boolContactMech.txt',',');

neighborsTouching = cell([size(pksBest,1),1]);
for j = 1:numel(npGood)
    np = npGood(j);
    neighborsTouching{np} = find(boolContactMech(np,:)==1);
end

% correct neighborsTouching...
% the following adjust neighborsTouching so that if particle 47 touches
% 1355, then particle 1355 also lists particle 47 as a contact 
for j = 1:numel(npGood)
    np = npGood(j);
    neighborsTouchingTemp = neighborsTouching{np}; %np = 47, neighborsTemp = [36,49,51,91,185,399,693,982,1355]
    for jj = 1:numel(neighborsTouchingTemp)
        npp = neighborsTouchingTemp(jj);
        neighborsTouching{npp} = union(neighborsTouching{npp},np); %neighbors{1355} = 690,982 becomes neighbors{1355} = union([690,982],47)
    end
end

Zmechanical = nan([size(pksBest,1),1]);
for j = 1:numel(npGood)
    np = npGood(j);
    Zmechanical(np) = numel(neighborsTouching{np});
end

cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen(sprintf('30may2014_1_30mayRCP_4x_id_x_y_z_Dcorr_packingFrac_Zfluorophore.txt'),'w');
for np = 1:size(pksBest,1)
    fprintf(fid,'%d %.04f %.04f %.04f %.04f %.04f %.04f\n',np,pksBest(np,2),pksBest(np,1),pksBest(np,3),Dcorr(np),packingFrac(np),Zmechanical(np));
end
fclose(fid);


% ==================================x

clearvars -except I* 

scaleFactor = 59.52381/1024; %um/pixels

cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen(sprintf('30may2014_1_30mayRCP_4x_id_x_y_z_Dcorr_packingFrac_Zfluorophore.txt'),'r');
temp = fscanf(fid,'%d %f %f %f %f %f %f\n',[7, Inf]);
fclose(fid);

order = temp(1,:)';

pksBest(order,1) = temp(3,:); %row = y, edited 10-Nov-2015, g(r) was computed using (x,y,z) instead of (y,x,z)
pksBest(order,2) = temp(2,:); %col = x
pksBest(order,3) = temp(4,:);

Dcorr(order) = temp(5,:);

packingFrac(order) = temp(6,:);

Zmechanical(order) = temp(7,:);

% 2mers:
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('2merALLcorrected30may2014_1_30mayRCP_4x_id_row_col_slice_K_D_w_offset.txt','r');
dimers = fscanf(fid,'%d %f %f %f %f %f %f %f\n',[8, Inf]);
fclose(fid);

npDimer = dimers(1,:)';                                                                                                                                                            

% 3mers:
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen('3merALLcorrected30may2014_1_30mayRCP_4x_id_row_col_slice_K_D_w_offset.txt','r');
trimers = fscanf(fid,'%d %f %f %f %f %f %f %f\n',[8, Inf]);
fclose(fid);

npTrimer = trimers(1,:)';  

npMultimer = vertcat(npDimer,npTrimer);

npBad = find(isnan(pksBest(:,1)) == 1);

npGood = 1:size(pksBest,1);
npGood(npBad) = [];

npSingle = 1:size(pksBest);
npSingle(union(npBad,npMultimer)) = [];

%-------- step1: identify particles near the center of the packing
DcorrAvg = nanmean(Dcorr(npSingle));

rowRange = round(1:DcorrAvg:size(IeditedRed,1));
colRange = round(1:DcorrAvg:size(IeditedRed,2));
sliceRange = round(1:DcorrAvg:size(IeditedRed,3));

centerImg = size(IeditedRed)/2;

distRow = abs(rowRange-centerImg(1));
rowCenterIdx = find(distRow == nanmin(distRow));

distCol = abs(colRange-centerImg(2));
colCenterIdx = find(distCol == nanmin(distCol));

distSlice = abs(sliceRange-centerImg(3));
sliceCenterIdx = find(distSlice == nanmin(distSlice));


rowMed = find((pksBest(:,1) > rowRange(rowCenterIdx-1)) & (pksBest(:,1) < rowRange(rowCenterIdx+1)));
colMed = find((pksBest(:,2) > colRange(colCenterIdx-1)) & (pksBest(:,2) < colRange(colCenterIdx+1)));
sliceMed = find(pksBest(:,3) > sliceRange(sliceCenterIdx-1) & pksBest(:,3) < sliceRange(sliceCenterIdx+1));

npCenterImg = intersect(sliceMed,union(rowMed,colMed));

npDelete = intersect(npCenterImg,npMultimer);
if(~isempty(npDelete))
    for k = 1:numel(npDelete)
        npCenterImg(npCenterImg == npDelete(k)) = [];
    end
end

% -----------step2: find number of particles in coord shells

stepSize = DcorrAvg/10000; %Aste et al. Phys Rev E 71 (2005) used anything between D/10,000 and D/100

delta = 0:stepSize:5;

% %The following if-conditional adjusts stepSize so that if the largest
% %partice neighbors the smallest particle, the stepSize will not be so large
% %that the smallest particle is overlooked...   
% if(stepSize >= (min(Dcorr(npGood))/max(Dcorr(npGood)))) 
%     stepSize = 0.5*(min(Dcorr(npGood))/max(Dcorr(npGood))); %reduce stepSize so that no particle is overlooked! Even in extreme cases where the smallest particle neighbors the largest particle  
% end
% delta = 0:stepSize:5;

numberDensity = nan([numel(delta),size(pksBest,1)]);% numberDensity(step,np) = # of particles/avgShellVol(step)
Ravg = nan([numel(delta),1]);
avgShellVol = nan([numel(delta),1]);

for step = 2:numel(delta)
    %sphericalShell thickness:
    %   from DcorrAvg*delta(step-1) to DcorrAvg*delta(step)
    %   therefore, median is: 
    %   DcorrAvg*delta(step) - DcorrAvg*((delta(step)-delta(step-1))/2)
    %
    %       stepSize = delta(step)-delta(step-1) 
    %   so,
    %   median simplified is:
    %   DcorrAvg*delta(step) - DcorrAvg*(stepSize/2)
    Ravg(step) = DcorrAvg*delta(step)-(DcorrAvg*(stepSize/2)); %average radial distance from particle np
    avgShellVol(step) = 4*pi*(Ravg(step)^2)*(DcorrAvg*stepSize); %volume of spherical shell, units: length^3
end

[dist, pksOrder] = orderPeaks(pksBest);

for j = 1:numel(npCenterImg)
    np = npCenterImg(j);
%     disp(sprintf('---------------pksBest(%d,:)-----------------',np))
    for step = 2:numel(delta)
        npTemp = find((dist(:,np) >= (delta(step-1)*DcorrAvg)) & (dist(:,np) < (delta(step)*DcorrAvg))); %particle centers located between a radial distance of delta(step-1)*DcorrAvg and delta(step)*DcorrAvg away from particle np
        
        npTemp(npTemp == np) = []; %delete instances of self-identification
        numberDensity(step,np) = numel(npTemp)/avgShellVol(step);
                
%         disp(sprintf('%d paricles in coordination shell from %.03f*<Dcorr> to %.03f*<Dcorr>, density = %.03f',numel(npTemp),delta(step-1),delta(step),numberDensity(step,np)))
    end
%     disp(sprintf('-----------------------x---------------------\n'))
end

% for jj = 1:numel(npCenterImg)
%     np = npCenterImg(jj);
%     
%     figure(2)
%     hold off
%     plot(delta,numberDensity(:,np))
%     title(sprintf('for particle %d \n numberDensity = #particles/avgShellVol \n where avgShellVol = (4*pi*r^2*deltaR)',np),'fontsize',20)
%     xlabel(sprintf('r/<Dcorr> \n <Dcorr> =%.04f',DcorrAvg),'fontsize',20)
% %     axis([delta(1),delta(numel(delta)),min([0,nanmin(pcf(npSingle))]),nanmax(pcf(npSingle))])
%     set(gca,'fontsize',20)
%     
%     pause
% end

% radial distribution function (aka pair correlation function) is just
% normalized numberDensity
% 
%   note: you may divide avgNumberDensity by anything but it's best to
%   divide by avgNumberDensity(numel(delta)) so that the normalized
%   function converges to 1 at delta(numel(delta))
% 
%   Saadtafar et al's JofMechanicsOfSolids_Vol62 used:
%       nanmean(packingFrac(npSingle))/volParticle


avgNumberDensity = nanmean(numberDensity(:,npCenterImg),2); 
avgPCF = avgNumberDensity./avgNumberDensity(numel(delta));

% % g(r)
figure(1)
hold off
plot(delta,avgPCF,'b-')
title(sprintf('average g_2(r) for %d single particles at center of packing (%d total particles) \n30may2014-1-30mayRCP-4x',numel(npCenterImg),numel(npGood)),'fontsize',20)
xlabel(sprintf('r/<Dcorr> \n <Dcorr> = %.04fum \n <Z-fluorophore> = %.04f \n <packingFrac> = %.04f',DcorrAvg*scaleFactor,nanmean(Zmechanical(npCenterImg)),nanmean(packingFrac(npCenterImg))),'fontsize',20)
ylabel('g_2(r)','fontsize',20)
axis([delta(1),delta(numel(delta)),min([0,nanmin(avgPCF(:))]),nanmax(avgPCF(:))+1])
set(gca,'fontsize',20)

% save datapts need to plot g2(r)
cd('~/poincareProgs/analysis/analysis30may2014_1_30mayRCP_4x/')
fid = fopen(sprintf('30may2014_1_30mayRCP_4x_fracOfDcorrAvg_avgPCF.txt'),'w');
for jj = 1:numel(delta)
    fprintf(fid,'%.04f %.04f\n',delta(jj),avgPCF(jj));
end
fclose(fid);


% -----------------------x
% % P(Dscaled)
DrangeNum = 100; %number of particle size from nanmin(DScaled(:)) to nanmax(DScaled(:))
diffDrange = 1;
numTry = 1;
while((numTry < 100) && (diffDrange ~= 0))
    [Dfrequency, Drange] = hist(DScaled(npCenterImg),DrangeNum);
    DrangeNum = numel(Dfrequency(Dfrequency > 0));
    diffDrange = numel(Drange) - DrangeNum; %check for convergence 
    numTry = numTry + 1;
end

Dfrequency(Dfrequency == 0) = nan; %just in case, the while loop didn't converge

figure(2)
hold off
plot(Drange,Dfrequency,'*','markersize',20)
title(sprintf('P(Diamter) \n 29may2014-2-28mayRefrig-4x2 \n %d single particles at center of packing (%d total)',numel(npCenterImg),numel(npGood)),'fontsize',20)
xlabel(sprintf('Diameter \n mean = %.03f +/- %.03fum \n manaufacturer: 5.06 +/- 0.44um',nanmean(DScaled(npCenterImg)),nanstd(DScaled(npCenterImg),1)),'fontsize',20)
ylabel('Frequency','fontsize',20)
axis([0,nanmax(DScaled(npSingle)),0,nanmax(Dfrequency(:))])
set(gca,'fontsize',20)

% -----------------------x
% % P(packingFrac)
packingFracTemp = packingFrac;
packingFracTemp(packingFracTemp >= 1) = nan;

rangeNum = 100; %number of particle size from nanmin(DScaled(:)) to nanmax(DScaled(:))
diffRange = 1;
numTry = 1;
while((numTry < 100) && (diffRange ~= 0))
    [packingFracFreq, packingFracRange] = hist(packingFrac(npCenterImg),rangeNum);
    rangeNum = numel(packingFracFreq(packingFracFreq > 0));
    diffRange = numel(packingFracRange) - rangeNum; %check for convergence 
    numTry = numTry + 1;
end

packingFracFreq(packingFracFreq == 0) = nan; %just in case, the while loop didn't converge

figure(3)
hold off
plot(packingFracRange,packingFracFreq,'*','markersize',20)
title(sprintf('P(packingFracLOCAL) \n 29may2014-2-28mayRefrig-4x2 \n %d single particles at center of packing (%d total)',numel(npCenterImg),numel(npGood)),'fontsize',20)
xlabel(sprintf('packingFracLocal \n mean = %.03f +/- %.03f \n mean(packingFrac < 1) = %.03f +/- %.03f',nanmean(packingFrac(npCenterImg)),nanstd(packingFrac(npCenterImg),1),nanmean(packingFracTemp(npCenterImg)),nanstd(packingFracTemp(npCenterImg),1)),'fontsize',20)
ylabel('Frequency','fontsize',20)
axis([0,nanmax(packingFrac(npSingle)),0,nanmax(packingFracFreq(:))])
set(gca,'fontsize',20)

% -----------------------x
% % P(Zmechanical)
ZmechTemp = Zmechanical;
ZmechTemp(ZmechTemp<2 | ZmechTemp >11) = nan;

% % % This is NOT necessary because there are not any fractions
% % rangeNum = 100; %number of particle size from nanmin(DScaled(:)) to nanmax(DScaled(:))
% % diffRange = 1;
% % numTry = 1;
% % while((numTry < 100) && (diffRange ~= 0))
% %     [ZmechFreq, ZmechRange] = hist(Zmech(npCenterImg),rangeNum);
% %     rangeNum = numel(ZmechFreq(ZmechFreq > 0));
% %     diffRange = numel(ZmechRange) - rangeNum; %check for convergence 
% %     numTry = numTry + 1;
% % end

ZmechRange = nanmin(Zmechanical(npCenterImg)):nanmax(Zmechanical(npCenterImg))+1;
[ZmechFreq,naught] = hist(Zmechanical(npCenterImg),ZmechRange);

ZmechFreq(ZmechFreq == 0) = nan; %just in case, the while loop didn't converge

figure(4)
hold off
plot(ZmechRange,ZmechFreq,'*','markersize',20)
title(sprintf('P(Z-fluorophore) \n 29may2014-2-28mayRefrig-4x2 \n %d single Particles near center of packing (% total)',numel(npCenterImg),numel(npGood)),'fontsize',20)
xlabel(sprintf('Z-fluorophore \n max = %.04f \n mean = %.04f, mean(2<=Z<=11) = %.04f \n min = %.04f \n %d have Z == 0',nanmax(Zmechanical(npCenterImg)),nanmean(Zmechanical(npCenterImg)),nanmean(ZmechTemp(npCenterImg)),nanmin(Zmechanical(npCenterImg)),numel(find(Zmechanical(npCenterImg)==0))),'fontsize',20)
ylabel(sprintf('# of particles'),'fontsize',20)
axis([0,nanmax(Zmechanical(npSingle)),0,nanmax(ZmechFreq(:))])
set(gca,'fontsize',20)
