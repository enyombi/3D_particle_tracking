% Author: Eru K.
% Date: 24-Aug-2017
% objective: List all the contacts found using on InormGreen

addpath('~/poincareProgs/particleTrackMatlab/');

% ----(start: load Red fluorescent bkgrd image from hard disk into RAM)----

cd('/Volumes/FreeAgent GoFlex Drive/Backups.backupdb/MacBook Pro/2017-03-24-140731/Macintosh HD/Users/ekyey001/poincareProgs/30may2014_2_29mayRLP_4x_Iedited');
dr=dir('*30may2014_2_29mayRLP_4x_Iedited*.txt');
sampleName = dr(1).name(1:23);
filesIedited = strcat(sampleName,'_Iedited');

IeditedRed = zeros([1024,1024,size(dr,1)],'single'); %size(I) == size(Iedited)

tic;
for slice = 1:size(dr,1)
    IeditedRed(:,:,slice) = dlmread(strcat(filesIedited,'_z',sprintf('%03d',slice),'.txt'),';');
    disp(sprintf('IeditedRed(:,:,%d), finished',slice))
end
elapsedTime = toc; %1.8747e+03seconds, 31.24min

% normalization for contact analysis...
max_value = max(IeditedRed(:));
min_value = min(IeditedRed(:));
InormRed = (IeditedRed - min_value)/(max_value-min_value);

InormRed(:,:,741:size(InormRed,3)) = []; %note: if min_value and/or max_value are only located in the deleted images, then min(InormRed(:)) > 0 and/or max(InormRed(:)) < 1

% ----(end: load Red fluorescent bkgrd image from hard disk into RAM)----

% ---(start: load Green fluorescent bkgrd image from hard disk into RAM)---

clearvars -except Inorm*

cd('/Volumes/FreeAgent GoFlex Drive/Backups.backupdb/MacBook Pro/2017-03-24-140731/Macintosh HD/Users/ekyey001/poincareProgs/30may2014_2_29mayRLP_4x_IeditedGreen');
dr=dir('*30may2014_2_29mayRLP_4x_IeditedGreen*.txt');
sampleName = dr(1).name(1:23);
filesIedited = strcat(sampleName,'_IeditedGreen');

IeditedGreen = zeros([1024,1024,size(dr,1)],'single'); %size(I) == size(Iedited)

tic;
for slice = 1:size(dr,1)
    IeditedGreen(:,:,slice) = dlmread(strcat(filesIedited,'_z',sprintf('%03d',slice),'.txt'),';');
    disp(sprintf('IeditedGreen(:,:,%d), finished',slice))
end
elapsedTime = toc; %1.8747e+03seconds, 31.24min

% normalization for contact analysis...
max_value = max(IeditedGreen(:));
min_value = min(IeditedGreen(:));
InormGreen = (IeditedGreen - min_value)/(max_value-min_value);


% ------------------------(particle data)-----------------------

clearvars -except Inorm* 

scaleFactor = 59.52381/1024; %um/pixels

cd('~/poincareProgs/analysis/analysis30may2014_2_29mayRLP_4x/')
fid = fopen(sprintf('30may2014_2_29mayRLP_4x_id_x_y_z_Dcorr_packingFrac_Zfluorophore.txt'),'r');
temp = fscanf(fid,'%d %f %f %f %f %f %f\n',[7, Inf]);
fclose(fid);

order = temp(1,:)';

pksBestOLD(order,1) = temp(3,:); %row = y, edited 10-Nov-2015, g(r) was computed using (x,y,z) instead of (y,x,z)
pksBestOLD(order,2) = temp(2,:); %col = x
pksBestOLD(order,3) = temp(4,:);

DcorrOLD(order) = temp(5,:);
packingFracOLD(order) = temp(6,:);

cd('~/poincareProgs/analysis/analysis30may2014_2_29mayRLP_4x/')
fid = fopen('pks30may2014_2_29mayRLP_4x_measurement00_voro.txt','r');
temp = fscanf(fid,'%d %f %f\n',[3, Inf]);
fclose(fid);

order = temp(1,:)';
% radii(order) = temp(2,:)';
voronoiVolOLD(order) = temp(3,:)';

% 2mers:
cd('~/poincareProgs/analysis/analysis30may2014_2_29mayRLP_4x/')
fid = fopen('2merALLcorrected30may2014_2_29mayRLP_4x_id_row_col_slice_K_D_w_offset.txt','r');
dimers = fscanf(fid,'%d %f %f %f %f %f %f %f\n',[8, Inf]);
fclose(fid);

npDimerOLD = [dimers(1,1:2:end)',dimers(1,2:2:end)']; % each row is one dimer

% 3mers:
cd('~/poincareProgs/analysis/analysis30may2014_2_29mayRLP_4x/')
fid = fopen('3merALLcorrected30may2014_2_29mayRLP_4x_id_row_col_slice_K_D_w_offset.txt','r');
trimers = fscanf(fid,'%d %f %f %f %f %f %f %f\n',[8, Inf]);
fclose(fid);

npTrimerOLD = [trimers(1,1:3:end)',trimers(1,2:3:end)',trimers(1,3:3:end)']; % each row is one trimer

particleID = 1:size(pksBestOLD,1);
particleID(isnan(pksBestOLD(:,1))) = []; %old npGood

% ---------------------------x

npReassigned = [particleID',[1:numel(particleID)]']; %[old listing, new listing]

% correct npDimer...
npDimer = npDimerOLD;
for ii = 1:numel(npDimer)
    npDimer(ii) = npReassigned(npReassigned(:,1) == npDimer(ii),2); %reassign each old particleID in npDimer to its new value  
end

% correct npTrimer...
npTrimer = npTrimerOLD;
for ii = 1:numel(npTrimer)
    npTrimer(ii) = npReassigned(npReassigned(:,1) == npTrimer(ii),2); %reassign each old particleID in npDimer to its new value  
end

pksBest = pksBestOLD(npReassigned(:,1),:);
Dcorr = DcorrOLD(npReassigned(:,1));
voronoiVol = voronoiVolOLD(npReassigned(:,1));

npSingle = 1:size(pksBest,1);
npSingle(union(npDimer(:),npTrimer(:))) = [];

% ------- Identify particles near the center of the packing ----------

clearvars -except Inorm* Dcorr pksBest voronoiVol npReassigned npSingle npDimer npTrimer scaleFactor map

DcorrAvg = nanmean(Dcorr(npSingle));

rowRange = round(1:DcorrAvg:1024);
colRange = round(1:DcorrAvg:1024);
sliceRange = round(1:DcorrAvg:740);

rowMed = find((pksBest(:,1) > rowRange(2)) & (pksBest(:,1) < rowRange(numel(rowRange)-1)));
colMed = find((pksBest(:,2) > colRange(2)) & (pksBest(:,2) < colRange(numel(colRange)-1)));
sliceMed = find((pksBest(:,3) > sliceRange(2)) & (pksBest(:,3) < sliceRange(numel(sliceRange)-1)));

npCenterImg = intersect(sliceMed,intersect(rowMed,colMed));
% numel(npCenterImg)
% 
% ans =
% 
%    583 <----------- N_{total}

npDelete = intersect(npCenterImg,union(npDimer(:),npTrimer(:))); % numel(npDelete) = 59 <----- N_{total}-N_{spherical}
if(~isempty(npDelete))
    for k = 1:numel(npDelete)
        npCenterImg(npCenterImg == npDelete(k)) = [];
    end
end

npGood = 1:npReassigned(end,2);
npGood(find(voronoiVol == 0)) = [];

tmpID = intersect(npGood,npCenterImg);

% [numel(tmpID),numel(npCenterImg)]
% 
% ans =
% 
%    524   524 <--------N_{spherical}

% ------------(D, packingFracLocal, and packingFracGlobal)-------------

clearvars -except Inorm* Dcorr pksBest voronoiVol npReassigned npSingle npDimer npTrimer scaleFactor map tmpID

% Computing packingFracGlobal----------
cd('~/poincareProgs/analysis/analysis30may2014_2_29mayRLP_4x/')
fid = fopen('30may2014_2_29mayRLP_4x_ID_IDneighbor_radiusID.txt','r');
B = fscanf(fid,'%d %d %f\n',[3, Inf]);
fclose(fid);

radiusMeasured(:,1) = B(1,:)';
radiusMeasured(:,2) = B(2,:)';
radiusMeasured(:,3) = B(3,:)';

numSample = nan([size(pksBest,1),1]);
radiusMeanSTD = nan([size(pksBest,1),2]); % mean, standard dev.
corrections = zeros([size(pksBest,1),1]); % number of outliers deleted
for ii = 1:size(npReassigned,1)
    np = npReassigned(ii,1); % OLD particle ID
    npNew = npReassigned(ii,2); % NEW particle ID
    
    numSample(npNew) = numel(radiusMeasured(radiusMeasured(:,1) == np,3));
    radiusSample = radiusMeasured(radiusMeasured(:,1) == np,3);

    radiusMeanSTD(npNew,1) = nanmean(radiusMeasured(radiusMeasured(:,1) == np,3));%nanmean(radiusSample)
    radiusMeanSTD(npNew,2) = nanstd(radiusMeasured(radiusMeasured(:,1) == np,3),0); %nanstd(x) = sqrt((1/(numel(x)-1))*nansum((x-mean(x)).^2)

    if(numSample(npNew)>3)%remove 1st outlier if there's at least >3 measurements
        errorRadius = abs(radiusSample-radiusMeanSTD(npNew,1))./radiusMeanSTD(npNew,1);
        
        radiusSample(errorRadius == nanmax(errorRadius)) = [];%shortens list of measurements by deleting 1st outlier        
        radiusMeanSTD(npNew,1) = nanmean(radiusSample);
        radiusMeanSTD(npNew,2) = nanstd(radiusSample);
        corrections(npNew) = corrections(npNew)+1;
        
        if(numSample(npNew)>4)%remove 2nd outlier if there's at least >4 measurements
            errorRadius = abs(radiusSample-radiusMeanSTD(npNew,1))./radiusMeanSTD(npNew,1);

            radiusSample(errorRadius == nanmax(errorRadius)) = [];%shortens list of measurements by deleting 2nd outlier       
            radiusMeanSTD(npNew,1) = nanmean(radiusSample);
            radiusMeanSTD(npNew,2) = nanstd(radiusSample);
            corrections(npNew) = corrections(npNew)+1;
        end 
    end
end

% change STD to standard error: 
% see g_r.m, normally standard_error =(standard_deviation)/sqrt(sampleSize)
% with Flaviano used standard_error = standard_deviation/sqrt(sampleSize-1) 
radiusMeanSTD(:,2) = radiusMeanSTD(:,2)./sqrt(numSample(:)-1);

% [nanmean(radiusMeanSTD(tmpID,1)),nanmean(radiusMeanSTD(tmpID,2))]*scaleFactor*2
% 
% ans =
% 
%     5.0512    0.0818
    
volParticle(:,1) = (4/3)*pi*((radiusMeanSTD(:,1)-radiusMeanSTD(:,2)).^3);%lower
volParticle(:,2) = (4/3)*pi*(radiusMeanSTD(:,1).^3);%mean
volParticle(:,3) = (4/3)*pi*((radiusMeanSTD(:,1)+radiusMeanSTD(:,2)).^3);%upper

packingFracGlobal(1) = nansum(volParticle(tmpID,1))/nansum(voronoiVol(tmpID));
packingFracGlobal(2) = nansum(volParticle(tmpID,2))/nansum(voronoiVol(tmpID));
packingFracGlobal(3) = nansum(volParticle(tmpID,3))/nansum(voronoiVol(tmpID));

% packingFracGlobal =
% 
%     0.5322    0.5591    0.5870
%     
% tmp = packingFracGlobal([1,3])-packingFracGlobal(2)
% 
% tmp =
% 
%    -0.0269    0.0279
% 
% nanmean(abs(tmp))
% 
% ans =
% 
%     0.0274
    
% Computing packingFracLocal----------
%       IMPORTANT: Dcorr./2 ~= radiusMeanSTD because 
%       Dcorr = radiusMeanSTD(:,1).*correctionFactor 
% 
%       You determined the correctionFactor by actually looking at 
%       each particle individually when you attend MRS in Boston. 
volSphere = (4/3)*pi*((Dcorr./2).^3);

packingFracLocal = volSphere./voronoiVol;
packingFracLocal(find(voronoiVol == 0)) = nan;

% [nanmean(packingFracLocal(tmpID)),nanstd(packingFracLocal(tmpID),0)/sqrt(numel(tmpID)-1)]
% 
% ans =
% 
%     0.6524    0.0062
    
% ---------------------(Read-in contact information)----------------------

clearvars -except Inorm* Dcorr pksBest voronoiVol npReassigned npSingle npDimer npTrimer scaleFactor map tmpID packingFrac*

cd('~/poincareProgs/analysis/analysis30may2014_2_29mayRLP_4x/')
fileID = fopen('30may2014_2_29mayRLP_4x_id_idTouch_contactIdx_contactRegSize_profileGreen_22-Jun-2017.txt');

clear contactListGreen
ii = 1; % row number
while ~feof(fileID)
    tmpLine = strsplit(fgetl(fileID),','); %note: fgetl() reads all data as string-type
    
    for jj = 1:numel(tmpLine) % jj == column number
        contactListGreen(ii,jj) = str2num(tmpLine{jj});
    end
    ii = ii + 1;
end
fclose(fileID);

% remove cross-listed values:
for row = 1:size(contactListGreen,1)
    np = contactListGreen(row,1);
    if(~isnan(np))
        npp = contactListGreen(row,2);
        idxCol = find(contactListGreen(:,1) == npp);
        rowDuplicate = idxCol(find(contactListGreen(idxCol,2) == np));
        contactListGreen(rowDuplicate,:) = nan;
    end
end
contactListGreen(isnan(contactListGreen(:,1)),:) = [];

% Remove contacts between particles that different parts of the SAME
% particle/mulitmer
for ii = 1:size(contactListGreen,1)
    np = contactListGreen(ii,1); % old particleID
    npTouch = contactListGreen(ii,2); % old particleID of neighboring particle

    npNew = npReassigned(npReassigned(:,1) == np,2); % NEW particleID
    npTouchNew = npReassigned(npReassigned(:,1) == npTouch,2); % NEW particleID of neighboring particle

    boolMultimer = 0; % default setting

    [rowD,colD] = find(npDimer == npNew);
    if(intersect(npDimer(rowD,:),npTouchNew)) % i.e., if np and npTouch occupy the same row of npDimer, then they are 2 parts of the SAME particle
        boolMultimer = 1;
    end

    [rowT,colT] = find(npTrimer == npNew);
    if(intersect(npTrimer(rowT,:),npTouchNew)) % i.e., if np and npTouch occupy the same row, then they are 2 parts of the SAME particle
        boolMultimer = 1;
    end
    
    if(boolMultimer == 1)
        contactListGreen(ii,:) = nan;
    end
end
contactListGreen(isnan(contactListGreen(:,1))==1,:) = [];

contactShortList = contactListGreen(contactListGreen(:,4) == 1,[1,2]); % the particle and its neighboring particle for contactRegion = 1


Zgreen = nan([size(npReassigned,1),1]);

% [nanmean(Zgreen(tmpID)),nanstd(Zgreen(tmpID),0)./sqrt(numel(tmpID)-1)]
% 
% ans =
% 
%     7.4122    0.0568

% % [ZfreqCenterImg,ZrangeCenterImg] = hist(Zgreen(npCenterImg),nanmin(Zgreen(npCenterImg)):nanmax(Zgreen(npCenterImg))); % Zrange = nanmin(Zgreen(:)):nanmax(Zgreen(:));
% % 
% % [Zfreq,Zrange] = hist(Zgreen(npSingle),nanmin(Zgreen(npSingle)):nanmax(Zgreen(npSingle))); % Zrange = nanmin(Zgreen(:)):nanmax(Zgreen(:));
% % 
% % figure(3)
% % hold off
% % h(1) = plot(Zrange,Zfreq,'b-*');
% % hold on
% % h(2) = plot(ZrangeCenterImg,ZfreqCenterImg,'r-*');
% % legendH = legend(h,'all spherical particles','all spherical particel AT center');
% % 
% % xlabel('Z','Interpreter','latex','fontsize',20)
% % ylabel('P(Z)','Interpreter','latex','fontsize',20)
% % 
% % legendH.FontName = 'Times New Roman';
% % legendH.FontSize = 12;
% % 
% % set(gca,'fontname','Times New Roman','fontsize',15) %change font of tickmarks and all other labels
% % 
% % cd('~/poincareProgs/analysis/packing-code/')
% % saveas(figure(3),'P_Z_30may2014_2_29mayRLP_4x','epsc')