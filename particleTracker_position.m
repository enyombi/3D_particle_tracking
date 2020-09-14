% Author: Eru K.
% Date: 22-May-2015
% objective: locate particles using Red Images (i.e., 12may2015_6may2015_3_1920A1000r_13xS3_Iedited)

addpath('~/poincareProgs/particleTrackMatlab/');

% -----------------(load Iedited from hard disk into RAM)-----------------
% load Iedited quickly via parallel computing (ref: particleTracker3.m)

cd('~/poincareProgs/12may2015_6may2015_3_1920A1000r_13xS3_Iedited/');
dr=dir('*12may2015_6may2015_3_1920A1000r_13xS3_Iedited*.txt');
sampleName = dr(1).name(1:37);
filesIedited = strcat(sampleName,'_Iedited');

clear job
clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

slicePerJob = 1:100:size(dr,1); %each task of a job is to load on image and each job consists of 100 tasks

for j = 1:length(slicePerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think good memory allocation to create only one job with many tasks delete the job after its complete to free of memory on the MJS queue

    startSlice = slicePerJob(j);
    stopSlice = startSlice+100-1; %subtract to account for MATLAB indexing displacement
   
    if(stopSlice > size(dr,1)) 
        stopSlice = size(dr,1);
    end
    
    for slice = startSlice:stopSlice
        createTask(job(j),@dlmread,1,{strcat(filesIedited,'_z',sprintf('%03d',slice),'.txt'),';'});
        disp(sprintf('job%02d,slice%03d',j,slice));
    end
end %2min24sec to finish 1 job and copy 100 images onto RAM

for j = 1:length(job)
    submit(job(j)); %~3min for 61 slices
    disp(sprintf('job%02d, submitted',j));
end %takes ~20min to load 'Iedited'

IeditedRed = zeros([1024,1024,size(dr,1)],'single'); %size(I) == size(IeditedRed)
% populate 'c' with the convolution matrix data and delete jobs from memory
for j = 1:length(job) %length(job) == number of jobs
    state = job(j).State;
    
    startSlice = slicePerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        results = fetchOutputs(job(j));
        for slice = 1:size(results,1)
            IeditedRed(:,:,startSlice+slice-1) = results{slice}; %minus 1 becuase MATLAB indexing starts at 1
            disp(sprintf('IeditedRed(:,:,%d)',startSlice+slice-1))
        end
    end
end 

% Shattuck's normalization...
max_value = max(IeditedRed(:));
min_value = min(IeditedRed(:));
InormRed = (max_value - IeditedRed)/(max_value-min_value);

% ----(end: load Red fluorescent bkgrd image from hard disk into RAM)----

% ---(start: load Green fluorescent bkgrd image from hard disk into RAM)---

clearvars -except Inorm*

cd('~/poincareProgs/12may2015_6may2015_3_1920A1000r_13xS3_IeditedGreen/');
dr=dir('*12may2015_6may2015_3_1920A1000r_13xS3_IeditedGreen*.txt');
sampleName = dr(1).name(1:37);
filesIedited = strcat(sampleName,'_IeditedGreen');

clear job
clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

slicePerJob = 1:100:size(dr,1); %each task of a job is to load on image and each job consists of 100 tasks

for j = 1:length(slicePerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think good memory allocation to create only one job with many tasks delete the job after its complete to free of memory on the MJS queue

    startSlice = slicePerJob(j);
    stopSlice = startSlice+100-1; %subtract to account for MATLAB indexing displacement
   
    if(stopSlice > size(dr,1)) 
        stopSlice = size(dr,1);
    end
    
    for slice = startSlice:stopSlice
        createTask(job(j),@dlmread,1,{strcat(filesIedited,'_z',sprintf('%03d',slice),'.txt'),';'});
        disp(sprintf('job%02d,slice%03d',j,slice));
    end
end %2min24sec to finish 1 job and copy 100 images onto RAM

for j = 1:length(job)
    submit(job(j)); %~3min for 61 slices
    disp(sprintf('job%02d, submitted',j));
end %takes ~20min to load 'Iedited'

IeditedGreen = zeros([1024,1024,size(dr,1)],'single'); %size(I) == size(IeditedGreen)
% populate 'c' with the convolution matrix data and delete jobs from memory
for j = 1:length(job) %length(job) == number of jobs
    state = job(j).State;
    
    startSlice = slicePerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        results = fetchOutputs(job(j));
        for slice = 1:size(results,1)
            IeditedGreen(:,:,startSlice+slice-1) = results{slice}; %minus 1 becuase MATLAB indexing starts at 1
            disp(sprintf('IeditedGreen(:,:,%d)',startSlice+slice-1))
        end
    end
end 

% Shattuck's normalization...
max_value = max(IeditedGreen(:));
min_value = min(IeditedGreen(:));
InormGreen = IeditedGreen/(max_value-min_value);

% ---(end: load Green fluorescent bkgrd image from hard disk into RAM)---

% ----------(Start: select one particle to use as maskImg)---------------

% use getPositionGreen.m to locate the center of the select particle

clearvars -except Iedited* Inorm*

for slice = 1:size(InormRed,3)
    figure(1)
    subplot(1,2,1)
    simage(InormRed(:,:,slice));
    title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)

%     subplot(1,2,2)
%     simage(InormGreen(:,:,slice));
%     title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
    pause
end

% center = [611,353,607];
center = [344,492,135];

[v,centerCutout] = getCutoutSectVertices(center,100,size(InormRed));
cutoutSectGreen = InormGreen(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));

[centerNew, radiusNew, validPtsEgdes] = getPositionGreen(v(1,:),centerCutout,cutoutSectGreen,85); %   350   499   140

clearvars -except Iedited* Inorm* center v centerCutout cutoutSectGreen centerNew radiusNew validPtsEgdes

centerCutoutNew = centerNew - v(1,:) + ones(1,3);

% [vNew,centerCutoutNew] = getCutoutSectVertices(centerNew,70,size(InormRed));
% cutoutSectGreenNew = InormGreen(vNew(1,1):vNew(7,1),vNew(1,2):vNew(7,2),vNew(1,3):vNew(7,3));

figure(1)        
subplot(1,2,1)
hold off
simage(cutoutSectGreen(:,:,centerCutoutNew(3)))
hold on
plot(centerCutoutNew(2),centerCutoutNew(1),'k*','markersize',20)
title(sprintf('cutoutSectGreen(:,:,%d)',centerCutoutNew(3)),'fontsize',20);
xlabel(sprintf('centerCutoutNew = [%d, %d, %d] (k*)',centerCutoutNew(1), centerCutoutNew(2), centerCutoutNew(3)),'fontsize',20);

subplot(1,2,2)
hold off
simage(cutoutSectGreen(:,:,centerCutout(3)))
hold on
plot(centerCutout(2),centerCutout(1),'k*','markersize',20)
title(sprintf('cutoutSectGreen(:,:,%d)',centerCutout(3)),'fontsize',20);
xlabel(sprintf('centerCutout = [%d, %d, %d] (k*)',centerCutout(1), centerCutout(2), centerCutout(3)),'fontsize',20);

% ====Left off here
for slice = 1:size(cutoutSectGreen,3)
    figure(2)        
    subplot(1,2,1)
    hold off
    simage(cutoutSectGreen(:,:,slice))
    hold on
    plot(centerCutoutNew(2),centerCutoutNew(1),'k*','markersize',20)
    title(sprintf('cutoutSectGreen(:,:,%d of %d)', slice, size(cutoutSectGreen,3)),'fontsize',20);
    xlabel(sprintf('centerCutoutNew = [%d, %d, %d] (k*)',centerCutoutNew(1), centerCutoutNew(2), centerCutoutNew(3)),'fontsize',20);
    
    if(slice == centerCutoutNew(3))
        plot(centerCutoutNew(2),centerCutoutNew(1),'m*','markersize',20)
    end
    
    subplot(1,2,2)
    hold off
    simage(cutoutSectGreen(:,:,slice))
    hold on
    plot(centerCutout(2),centerCutout(1),'k*','markersize',20)
    title(sprintf('cutoutSectGreen(:,:,%d of %d)', slice, size(cutoutSectGreen,3)),'fontsize',20);
    xlabel(sprintf('centerCutout = [%d, %d, %d] (k*)',centerCutout(1), centerCutout(2), centerCutout(3)),'fontsize',20);
    
    if(slice == centerCutout(3))
        plot(centerCutout(2),centerCutout(1),'m*','markersize',20)
    end
    pause
end

% ---------------------------------

clearvars -except Iedited* Inorm*

% now use the center of the selected particle to find its profile in
% InormRed 

center = [344,492,135];

[v,centerCutout] = getCutoutSectVertices(center,70,size(InormRed));
cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));

[v,naught] = getCutoutSectVertices(center,70,size(InormRed));
cutoutSectGreen = InormGreen(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));


for slice = 1:size(cutoutSectGreen,3)
    figure(3)        
    subplot(1,2,1)
    hold off
    simage(cutoutSectGreen(:,:,slice))
    hold on
    plot(centerCutout(2),centerCutout(1),'k*','markersize',20)
    title(sprintf('cutoutSectGreen(:,:,%d of %d)', slice, size(cutoutSectGreen,3)),'fontsize',20);
    xlabel(sprintf('centerCutout = [%d, %d, %d] (k*)',centerCutout(1), centerCutout(2), centerCutout(3)),'fontsize',20);
    
    if(slice == centerCutout(3))
        plot(centerCutout(2),centerCutout(1),'m*','markersize',20)
    end
    
    subplot(1,2,2)
    hold off
    simage(cutoutSect(:,:,slice))
    hold on
    plot(centerCutout(2),centerCutout(1),'k*','markersize',20)
    title(sprintf('cutoutSect(:,:,%d of %d)', slice, size(cutoutSect,3)),'fontsize',20);
    xlabel(sprintf('centerCutout = [%d, %d, %d] (k*)',centerCutout(1), centerCutout(2), centerCutout(3)),'fontsize',20);
    
    if(slice == centerCutout(3))
        plot(centerCutout(2),centerCutout(1),'m*','markersize',20)
    end
    pause
end


% step1: 
% look at the image and guess the diameter of the particle 
% There are markers on the image located at a distance 'diameter' from
% 'centerCutout' to help with this  
diameter = 110; 

figure(4) 
hold off
simage(cutoutSect(:,:,centerCutout(3)));
hold on
plot(centerCutout(2),centerCutout(1),'k*','markersize',20)
plot(centerCutout(2)-(diameter/2),centerCutout(1),'m*','markersize',20)
plot(centerCutout(2)+(diameter/2),centerCutout(1),'m*','markersize',20)
plot(centerCutout(2),centerCutout(1)-(diameter/2),'m*','markersize',20)
plot(centerCutout(2),centerCutout(1) + (diameter/2),'m*','markersize',20)

% step2: Guessing 'w' in ssf.m 
% a) get 'avgProfile' of the particle using 'diameter'.
% b) plot 'avgProfile'
% c) looking at plot of avgProfile, adjust w until it properly marks off 
%    the steepest pts of avgProfile as in figure 1 of
%    http://gibbs.engr.ccny.cuny.edu/technical/Tracking/ChiTrack.php   

% step2a
avgProfile = getSphericalProfile(centerCutout,diameter,0,cutoutSect);

% step2b,2c
w = 40/2; %adjust this value
stopW = 4; %distance from the edge of the particle where the region (2*w) of 'avgProfile' ends

figure(5) 
hold off
g(1) = plot(avgProfile,'r-*');
hold on
g(2) = plot(stopW:ceil(w),avgProfile(stopW:ceil(w)),'b*');
g(3) = plot((numel(avgProfile)-ceil(w)):(numel(avgProfile)-stopW),avgProfile((numel(avgProfile)-ceil(w)):(numel(avgProfile))-stopW),'b*');



D = diameter;
K = max(avgProfile(:))-min(avgProfile(:));
% Imax = (K/2)*(1-tanh(abs(r_max)-D/2)/w) + Imin
%   I is symmetric so r_max = 0 and tanh(abs(0)-D/2)/w) ~ -1 
%
% so,
% Imax = (K/2)*(1-(-1)) + Imin
% Imax - Imin = (K/2)*2
% K = Imax - Imin
offset = (min(avgProfile((length(avgProfile)-stopW):length(avgProfile)))+min(avgProfile(1:stopW)))/2;
b = ceil((length(avgProfile)-1)/2);
r = linspace(-b,b,length(avgProfile));

param(1,:) = nlinfit(r,avgProfile,'ssf',[K,D,w,offset]); 

iter = 2;
while((param(iter-1,4) < offset*0.8 || param(iter-1,4) > offset*1.8)  && iter < 150)
    param(iter,:) = nlinfit(r,avgProfile,'ssf',[K,param(iter-1,2),param(iter-1,3),offset]); 
    disp(sprintf('iteration%d; K = %.04f, D = %.04f, w = %.04f, offset = %.04f',iter,param(iter,1),param(iter,2),param(iter,3),param(iter,4)));
    iter = iter + 1; 
end

error = abs(param(:,4) - offset)./offset;
iterBest = min(find(error(:) == min(error)));


iter = 2;
while((param(iter-1,4) < offset*0.8 || param(iter-1,4) > offset*1.8)  && iter < 150)
    param(iter,:) = nlinfit(r,avgProfile,'ssf',[K,param(iter-1,2),param(iter-1,3),offset]); 
    disp(sprintf('iteration%d; K = %.04f, D = %.04f, w = %.04f, offset = %.04f',iter,param(iter,1),param(iter,2),param(iter,3),param(iter,4)));
    iter = iter + 1; 
end

error = abs(param(:,4) - offset)./offset;
iterBest = min(find(error(:) == min(error)));

% param =
% 
%    0.1248355  70.5648041  12.4865561   0.5861360


Kactual = param(iterBest,1);
Dactual = param(iterBest,2);
wActual = param(iterBest,3);
offsetActual = param(iterBest,4);

ss=2*fix((Dactual/2)+(2*wActual))-1;  % size of ideal particle image NOT the size of the particle in the image
os=(ss-1)/2; % (size-1)/2 of ideal particle image
center = [os+1 os+1 os+1]; %center of ideal particle image
[voronoiVol, Rmask] = peakPlacementB(center,[ss ss ss],Dactual,wActual);
maskExp = getCalcImg(Rmask,Kactual,Dactual,wActual,offsetActual,voronoiVol);%the 'actual' maskImg to be use in particleTracking algorithm

% Crucial Note: 
%   peakPlacementB().m uses the default: particleRadius = (D/2)+(2*w)
%   but this is largely an ARBITRARY empiricism (i.e., you use this because
%   it's what Shattuck used in computing the size of the mask image 'ss'
%   and because it usually provides an accurate intensity of a real image
%   of a particle) 
avgProfileActual = getSphericalProfile(center,Dactual,wActual,maskExp); %default 
avgProfileActual2 = getSphericalProfile(center,diameter,0,maskExp); %default 

avgProfileCorr = interp1(1:length(avgProfile),avgProfile,linspace(1,length(avgProfile),length(avgProfileActual)));
b = ceil((max(numel(avgProfile),numel(avgProfileActual))-1)/2);
r = linspace(-b,b,max(numel(avgProfile),numel(avgProfileActual)));
diff = (avgProfileActual - avgProfileCorr).^2;
Rsquared = sum(diff(:));
    
figure(6)
h(1) = plot(r,avgProfileActual,'b--')
hold on
h(2) = plot(r,avgProfileCorr,'k*')
legend(h,'maskExp, obtained using defualt (Dactual/2)+(2*wActual)','real data, obtained using ONLY D')
title(sprintf('Intensity Profile \n sum((maskExp - rawExp).^2) = %.04f',Rsquared),'fontsize',20)

b = ceil((numel(avgProfile)-1)/2); %numel(avgProfile) == numel(avgProfileActual2) becuase getSpericlaProfile(..,diameter,0,...) 
r = linspace(-b,b,numel(avgProfile));
diff = (avgProfileActual2 - avgProfile).^2;
Rsquared = sum(diff(:));

figure(7)
h(1) = plot(r,avgProfileActual2,'b--')
hold on
h(2) = plot(r,avgProfile,'k*')
legend(h,'maskExp, obtained using ONLY D','real data, obtained using ONLY D')
title(sprintf('Intensity Profile \n sum((maskExp - rawExp).^2) = %.04f',Rsquared),'fontsize',20)























% --------------------(start: thresholding)-------------------------
clearvars -except Inorm*

Kactual = 0.1248355; 
Dactual = 70.5648041;  
wActual = 12.4865561;   
offsetActual = 0.5861360;

ss=2*fix((Dactual/2)+(2*wActual))-1;  % size of ideal particle image NOT the size of the particle in the image
os=(ss-1)/2; % (size-1)/2 of ideal particle image
center = [os+1 os+1 os+1]; %center of ideal particle image
[voronoiVol, Rmask] = peakPlacementB(center,[ss ss ss],Dactual,wActual);
maskExp = getCalcImg(Rmask,Kactual,Dactual,wActual,offsetActual,voronoiVol);%the 'actual' maskImg to be use in particleTracking algorithm

N = size(InormRed) + size(maskExp) - 1;
a = single(fftn(InormRed,N));
b = single(fftn(maskExp,N));

prod = @(x,y) x.*y;

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

numTasks = 100; %one slice == one task
tasksPerJob = 1:numTasks:N(3);
for j = 1:length(tasksPerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue

    startSlice = tasksPerJob(j);
    stopSlice = startSlice+numTasks-1; %subtract to account for MATLAB indexing displacement
   
    if(stopSlice > N(3)) 
        stopSlice = N(3);
    end
    
    for slice = startSlice:stopSlice %slice == task #
        createTask(job(j),prod,1,{a(:,:,slice),b(:,:,slice)});
        disp(sprintf('job%02d,slice%03d',j,slice))
    end
end

% submit one job at a time to be computed by all available slave nodes
clear a b %first clear RAM of 'a' and 'b' to make room for computation of 'c'
for j = 1:length(job)
    submit(job(j)); %~3min for 61 slices
    disp(sprintf('job%02d, submitted',j))
end 

% populate 'c' with the convolution matrix data and delete jobs from memory
cFourier = zeros(N,'single');
for j = 1:length(job)
    state = job(j).State;
    
    startSlice = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        results = fetchOutputs(job(j));
        for slice = 1:size(results,1)
            cFourier(:,:,startSlice+slice-1) = results{slice}; %minus 1 becuase MATLAB indexing starts at 1
            disp(sprintf('c(:,:,%d)',startSlice+slice-1))
        end
    end
end

c = ifftn(cFourier,N);

clear cFourier

RowColStart = ((size(maskExp,1)-1)/2)+1;%plus 1 to correct for matlab indexing/counting which starts at 1 NOT 0 
RowColStop = N(1)-(RowColStart-1); %minus 1 to correct for matlab indexing/counting which starts at 1 NOT 0 

SliceStart = RowColStart;
SliceStop = N(3)-(RowColStart-1);

cScaled = c(RowColStart:RowColStop,RowColStart:RowColStop,SliceStart:SliceStop);

% ------------(stop: convolution of Inorm & maskExp)--------------------


% --------------------(start: thresholding)-------------------------

clearvars -except cScaled Inorm* maskExp Kactual Dactual wActual offsetActual 

sizeSect = round(size(maskExp,1)); 
rangeRow = 1:sizeSect:size(InormRed,1); %size(InormRed,1) = 1024pixels
rangeCol = 1:sizeSect:size(InormRed,2); %size(InormRed,1) = 1024pixels
rangeSlice = 1:sizeSect:size(InormRed,3); %size(InormRed,3) = 795pixels

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

counter = 0;
for slice = 1:length(rangeSlice)
    job(slice) = createJob(clusterD);
    disp(sprintf('\n------------job%02d, created with the following tasks------------',slice));
    for col = 1:length(rangeCol)  
        for row = 1:length(rangeRow)
            stopRow = rangeRow(row)+sizeSect-1;
            stopCol = rangeCol(col)+sizeSect-1;
            stopSlice = rangeSlice(slice)+sizeSect-1;
            
            if(stopRow > size(InormRed,1))
                stopRow = size(InormRed,1);
            end
            
            if(stopCol > size(InormRed,2))
                stopCol = size(InormRed,2);
            end
            
            if(stopSlice > size(InormRed,3))
                stopSlice = size(InormRed,3);
            end
            counter = counter + 1;
            disp(sprintf('thresholding...section %d, InormRed(%d:%d,%d:%d,%d:%d)',counter,rangeRow(row),stopRow,rangeCol(col),stopCol,rangeSlice(slice),stopSlice))
            cutoutSect = cScaled(rangeRow(row):stopRow,rangeCol(col):stopCol,rangeSlice(slice):stopSlice);
            createTask(job(slice),@threshold,1,{cutoutSect,0.999,size(InormRed),[rangeRow(row), rangeCol(col), rangeSlice(slice)]});
        end
    end
end
clear maskExp

% submit one job at a time to computed by all available slave nodes
for num = 1:length(job)
    submit(job(num));
    disp(sprintf('job%02d, submitted',num));
end

% populate tPks...
numJob = 1;
numSect = 1;
pks = fetchOutputs(job(numJob));
thresholdPks = pks{numSect};
disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
for numSect = 2:length(pks)
    thresholdPks = vertcat(thresholdPks,pks{numSect});
    disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
end

for numJob = 1:length(job)
    pks = fetchOutputs(job(numJob));
    for numSect = 1:length(pks)
        thresholdPks = vertcat(thresholdPks,pks{numSect});
        disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
    end
end
thresholdPks = single(thresholdPks);
[tPks(:,1), tPks(:,2), tPks(:,3)] = ind2sub(size(InormRed),thresholdPks);
% --------------------(end: thresholding)-------------------------



% ----(start: eliminate mulitple peaks identifying the same particle)----
clearvars -except tPks sizeSect Inorm* Kactual Dactual wActual offsetActual 

% sizeSect = size(maskExp,1);
numBoxes = 2; %OVERESTIMATE THE NUMBER OF PARTICLES BY MAKING THE number of boxes with particles smaller (i.e., sizeSect) smaller so that there are more boxes!  round(size(maskExp,1)/2) makes each box half the size of the mask image and thereby assumes that there are 2 particles (mask images) per box

rangeRow = 1:(numBoxes*sizeSect):size(InormRed,1); %size(InormRed,1) = 1024pixels
rangeCol = 1:(numBoxes*sizeSect):size(InormRed,2); %size(InormRed,1) = 1024pixels
rangeSlice = 1:(numBoxes*sizeSect):size(InormRed,3); %size(InormRed,3) = 795pixels

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

counter = 0;

for slice = 1:length(rangeSlice)
    job(slice) = createJob(clusterD); %1024x1024x795pixels == 9x9x7sections for size(maskExp) = [119,119,119]. 1section == 1 task. 9x9sections or 81 sections per job 
    disp(sprintf('\n------------job%02d, created with the following tasks------------',slice));
    for col = 1:length(rangeCol)  
        for row = 1:length(rangeRow)
            stopRow = rangeRow(row)+(sizeSect*numBoxes)-1;
            stopCol = rangeCol(col)+(sizeSect*numBoxes)-1;
            stopSlice = rangeSlice(slice)+(sizeSect*numBoxes)-1;
            
            if(stopRow > size(InormRed,1))
                stopRow = size(InormRed,1);
            end
            
            if(stopCol > size(InormRed,2))
                stopCol = size(InormRed,2);
            end
            
            if(stopSlice > size(InormRed,3))
                stopSlice = size(InormRed,3);
            end
            
            corrSlice = find(tPks(:,3) >= rangeSlice(slice) & tPks(:,3) <= stopSlice);
            corrCol = find(tPks(:,2) >= rangeCol(col) & tPks(:,2) <= stopCol);
            corrRow = find(tPks(:,1) >= rangeRow(row) & tPks(:,1) <= stopRow);
        
            idxPks = intersect(intersect(corrSlice,corrCol),corrRow);%thresholding peaks located inside cScaled(rangeRow(row):stopRow,rangeCol(col):stopCol,rangeSlice(slice):stopSlice)
            
            counter = counter + 1;
            disp(sprintf('locating peaks %0.2fpixels apart inside...LARGE section %d, InormRed(%d:%d,%d:%d,%d:%d)',Dactual,counter,rangeRow(row),stopRow,rangeCol(col),stopCol,rangeSlice(slice),stopSlice));
            createTask(job(slice),@getNeighborsValid,1,{tPks(idxPks,:),Dactual});
        end
    end
end

% submit one job at a time to computed by all available slave nodes
for num = 1:length(job)
    submit(job(num));
    disp(sprintf('job%02d, submitted',num));
end

% populate validPks...
numJob = 1;
numSect = 1;
pks = fetchOutputs(job(numJob));
validPks = pks{numSect};
disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
for numSect = 2:length(pks)
    validPks = vertcat(validPks,pks{numSect});
    disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
end

for numJob = 1:length(job)
    pks = fetchOutputs(job(numJob));
    for numSect = 1:length(pks)
        validPks = vertcat(validPks,pks{numSect});
        disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
    end
end
validPks = single(validPks);

% validPks =[
% 
%    117    88    60
%    225    60    60
%     65    60   236
%    129   119   214
%    161   221   238
%    241    60    60
%    409    98    60
%    473   238    60
%    353   119   237
%    321   236   238
%    597    60    60
%    477   237    60
%    713   224    60
%    481   119   238
%    593    60   236
%    829    60    60
%    749   216    60
%    837   237    60
%    769   116   238
%    833    60   238
%    721   203   238
%    965    60    60
%    965   198    60
%    961   119   236
%     61   355    60
%    237   336    60
%    237   472    60
%    113   239   238
%     97   475   238
%    241   337    60
%    465   246    60
%    241   464    60
%    361   453    60
%    321   242   238
%    241   475   216
%    385   476   237
%    477   244    60
%    593   319    60
%    665   280    60
%    481   299   237
%    593   470   238
%    641   390   237
%    749   242    60
%    841   265    60
%    829   476    60
%    881   358    60
%    721   241   238
%    833   279   204
%    721   389   238
%    945   476   238
%    965   241    60
%    965   470    60
%    961   357   237
%    961   476   235
%    121   544    60
%    237   487    60
%    177   606    60
%     97   480   236
%     65   712   238
%    177   633   238
%    241   482    60
%    357   530    60
%    357   630    60
%    241   486   215
%    321   507   209
%    385   499   238
%    241   599   238
%    353   596   194
%    593   517    60
%    573   714    60
%    657   713    60
%    577   529   238
%    593   693   236
%    833   482    60
%    945   550    60
%    813   707    60
%    949   643    60
%    785   595   230
%    737   714   238
%    953   554    60
%    953   641    60
%    961   565   238
%     61   763    60
%    213   720    60
%     61   834    60
%    221   919    60
%     65   755   235
%     65   834   237
%    241   717    60
%    469   833    60
%    357   952    60
%    241   795   183
%    241   952   234
%    477   833    60
%    577   726    60
%    657   740    60
%    713   835    60
%    481   833   238
%    593   716   222
%    481   927   237
%    817   719    60
%    949   721    60
%    829   952    60
%    737   726   238
%    721   883   237
%    833   932   238
%    953   736    60
%    953   949    60
%    961   715   216
%    961   934   234
%    121   965    60
%    233   953    60
%    113   965   237
%    241   955    60
%    353   965    60
%    241   965   238
%    481   954    60
%    601   965    60
%    481   953   238
%    593   965   238
%    833   965    60
%    833   953   237
%    953   965    60
%    961   953   229
%    117    88    60
%    225    60    60
%     65    60   236
%    129   119   214
%    161   221   238
%    241    60    60
%    409    98    60
%    473   238    60
%    353   119   237
%    321   236   238
%    597    60    60
%    477   237    60
%    713   224    60
%    481   119   238
%    593    60   236
%    829    60    60
%    749   216    60
%    837   237    60
%    769   116   238
%    833    60   238
%    721   203   238
%    965    60    60
%    965   198    60
%    961   119   236
%     61   355    60
%    237   336    60
%    237   472    60
%    113   239   238
%     97   475   238
%    241   337    60
%    465   246    60
%    241   464    60
%    361   453    60
%    321   242   238
%    241   475   216
%    385   476   237
%    477   244    60
%    593   319    60
%    665   280    60
%    481   299   237
%    593   470   238
%    641   390   237
%    749   242    60
%    841   265    60
%    829   476    60
%    881   358    60
%    721   241   238
%    833   279   204
%    721   389   238
%    945   476   238
%    965   241    60
%    965   470    60
%    961   357   237
%    961   476   235
%    121   544    60
%    237   487    60
%    177   606    60
%     97   480   236
%     65   712   238
%    177   633   238
%    241   482    60
%    357   530    60
%    357   630    60
%    241   486   215
%    321   507   209
%    385   499   238
%    241   599   238
%    353   596   194
%    593   517    60
%    573   714    60
%    657   713    60
%    577   529   238
%    593   693   236
%    833   482    60
%    945   550    60
%    813   707    60
%    949   643    60
%    785   595   230
%    737   714   238
%    953   554    60
%    953   641    60
%    961   565   238
%     61   763    60
%    213   720    60
%     61   834    60
%    221   919    60
%     65   755   235
%     65   834   237
%    241   717    60
%    469   833    60
%    357   952    60
%    241   795   183
%    241   952   234
%    477   833    60
%    577   726    60
%    657   740    60
%    713   835    60
%    481   833   238
%    593   716   222
%    481   927   237
%    817   719    60
%    949   721    60
%    829   952    60
%    737   726   238
%    721   883   237
%    833   932   238
%    953   736    60
%    953   949    60
%    961   715   216
%    961   934   234
%    121   965    60
%    233   953    60
%    113   965   237
%    241   955    60
%    353   965    60
%    241   965   238
%    481   954    60
%    601   965    60
%    481   953   238
%    593   965   238
%    833   965    60
%    833   953   237
%    953   965    60
%    961   953   229
%     65    60   350
%    129   119   240
%    129   224   330
%    129    60   380
%    353    60   355
%    321   238   254
%    385   146   284
%    353   156   359
%    481    60   356
%    705    60   357
%    481   122   281
%    705   220   314
%    481    81   446
%    545    60   474
%    705    60   441
%    673   120   464
%    769   111   306
%    833    79   328
%    737   184   271
%    865    80   475
%    929   237   476
%    961   118   261
%    961   115   475
%     97   276   357
%    161   239   309
%     97   476   264
%    129   476   357
%    129   366   442
%    129   476   428
%    321   239   254
%    353   476   242
%    449   474   357
%    353   239   358
%    385   334   363
%    289   405   431
%    481   291   327
%    705   239   297
%    481   475   357
%    705   409   354
%    481   465   430
%    609   423   476
%    737   239   298
%    833   357   340
%    769   398   356
%    929   262   476
%    961   357   357
%    961   476   356
%    961   256   476
%     97   513   293
%    129   516   356
%     97   686   353
%    193   632   312
%    129   482   419
%    225   551   476
%    129   714   472
%    353   505   241
%    449   482   355
%    257   596   357
%    257   553   473
%    417   596   466
%    481   478   357
%    609   479   258
%    481   596   357
%    609   678   357
%    673   708   340
%    673   494   474
%    609   713   476
%    769   477   357
%    737   511   472
%    929   688   476
%    961   511   357
%    961   596   262
%    961   525   476
%    961   679   474
%     65   778   356
%    225   952   357
%    129   754   473
%    129   834   464
%    193   951   473
%    353   801   357
%    257   952   357
%    353   746   475
%    385   830   455
%    481   832   357
%    673   715   326
%    481   918   329
%    705   951   357
%    609   729   476
%    705   952   476
%    929   833   357
%    833   922   310
%    833   833   473
%    929   830   473
%    897   914   360
%    961   833   357
%    961   829   453
%    961   913   475
%    129   965   357
%    225   965   357
%    129   953   475
%    225   965   476
%    257   965   357
%    481   953   306
%    609   965   296
%    705   965   476
%    833   953   297
%    737   965   476
%    961   953   242
%    961   965   476
%     65    60   477
%    193   119   594
%    129   174   595
%     65   117   714
%    193   182   626
%    385   119   595
%    257   176   595
%    449   215   593
%    257    60   714
%    449    78   713
%    321   142   662
%    449   238   712
%    481    85   477
%    545    60   498
%    705    60   478
%    577   151   594
%    577   119   712
%    705   122   705
%    833    60   551
%    769   120   525
%    897   238   595
%    769    60   597
%    897    60   711
%    897   237   687
%    961    98   567
%    961   238   594
%    961    60   711
%    961   238   706
%    129   357   477
%    129   291   594
%    225   296   511
%     65   322   714
%    193   357   714
%    257   288   525
%    449   356   595
%    321   476   594
%    257   356   714
%    449   290   713
%    321   421   712
%    385   476   714
%    449   442   714
%    577   357   594
%    705   357   478
%    481   457   477
%    577   419   546
%    641   357   713
%    705   270   714
%    705   469   714
%    897   265   594
%    737   476   479
%    897   476   594
%    897   256   672
%    961   270   595
%    961   405   594
%    961   258   670
%    961   407   713
%    129   478   477
%    225   545   511
%    129   714   495
%    193   713   595
%    129   477   714
%    129   698   713
%    257   548   513
%    321   565   594
%    385   596   595
%    321   577   679
%    385   572   714
%    449   679   713
%    513   516   594
%    641   593   595
%    577   696   575
%    577   528   714
%    705   542   706
%    513   684   713
%    737   513   490
%    833   714   595
%    897   669   714
%    961   531   506
%    961   597   540
%    961   477   707
%    961   688   714
%    129   755   494
%    129   834   477
%    193   952   483
%    129   768   713
%     65   887   713
%    193   952   596
%    353   747   505
%    385   827   479
%    385   720   596
%    449   716   714
%    257   952   596
%    449   834   714
%    577   717   519
%    577   938   595
%    705   951   594
%    577   716   712
%    577   877   714
%    833   825   523
%    929   828   477
%    769   893   539
%    833   724   614
%    833   920   713
%    961   826   480
%    961   906   546
%    961   715   714
%    129   953   519
%    193   965   486
%    129   953   596
%    385   965   585
%    257   965   596
%    577   953   595
%    705   965   528
%    705   965   695
%    737   965   491
%    897   965   589
%    833   953   714
%    961   965   544
%     65   116   774
%    257    60   774
%    449    70   743
%    449   238   774
%    577   117   774
%    705   192   774
%    833    60   774
%    897    95   774
%    769   196   774
%    897   237   717
%    961   118   774
%    961   235   716
%    961   194   774
%     65   356   774
%    193   357   774
%    257   357   774
%    449   287   774
%    385   474   774
%    449   444   722
%    577   356   774
%    705   244   774
%    513   453   728
%    897   242   717
%    769   475   774
%    961   244   715
%    961   357   774
%    129   477   774
%     65   713   774
%    193   714   716
%    257   594   774
%    321   572   718
%    385   495   774
%    449   712   774
%    577   535   738
%    513   713   774
%    577   714   739
%    833   594   774
%    897   668   774
%    961   579   774
%    961   487   774
%    961   664   774
%     65   717   774
%    193   773   718
%     65   882   774
%    257   729   774
%    449   765   774
%    257   834   774
%    449   836   774
%    513   716   774
%    577   766   727
%    577   877   737
%    833   833   753
%    897   715   774
%    833   936   774
%    961   715   774
%    961   952   774
%    129   953   774
%    257   965   774
%    449   953   774
%    577   953   715
%    705   965   774
%    769   965   774
%    961   965   774];

% ----(end: eliminate mulitple peaks identifying the same particle)----




% --(start: move peaks to minimize error and locate at pixel resolution)--

clearvars -except pksNewer K D w offset Inorm* Kactual Dactual wActual offsetActual

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

numTasks = 6; %one particle == one task
tasksPerJob = 1:numTasks:size(validPks,1);
for j = 1:length(tasksPerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue
    disp(sprintf('\n------------job%02d, created with the following tasks------------',j));

    startTask = tasksPerJob(j);
    stopTask = startTask+numTasks-1; %subtract to account for MATLAB indexing displacement
   
    if(stopTask > size(validPks,1)) 
        stopTask = size(validPks,1);
    end
    
    for np = startTask:stopTask %slice == task #
        [v,center] = getCutoutSectVertices(validPks(np,:),2*Dactual,size(InormRed)); %vertices of rectangular portion of 'InormRed' that is (2*Dactual)^3 and centered on validPks(np,:) 
        createTask(job(j),@getPosition,2,{v(1,:),center,Kactual,Dactual,wActual,offsetActual,InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3))});
        disp(sprintf('locating validPks(%d,:) at pixel resolution...',np));
    end
end

% submit one job at a time to be computed by all available slave nodes
for j = 1:length(job)
    submit(job(j));
    disp(sprintf('job%02d, submitted',j))
end 
% Actual time:
%   first job started at: Mon Jul 27 20:19:44 EDT 2015
%   last job finished at: Tue Jul 28 13:20:01 EDT 2015
%   total time: 17hrs01min17sec for locating 332 particles

pksNew = zeros(size(validPks),'single');
Rsquared = zeros([8,size(validPks,1)],'single'); 
for j = 1:length(job)
    state = job(j).State;
    
    startTask = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        results = fetchOutputs(job(j));
        disp(sprintf('\n------------job%02d, finished with the following tasks------------',j));
        for np = 1:size(results,1)
            pksNew(startTask+np-1,:) = results{np,1}; %minus 1 becuase MATLAB indexing starts at 1
            Rsquared(:,startTask+np-1) = results{np,2};
            disp(sprintf('validPks(%d,:) becomes pksNew(%d,:)',startTask+np-1,startTask+np-1));
        end
    end
end 

pksNew = zeros(size(validPks),'single');
Rsquared = zeros([8,size(validPks,1)],'single'); 
for j = 1:length(job)
    state = job(j).State;
    
    startTask = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        try 
            results = fetchOutputs(job(j));
            disp(sprintf('\n------------job%02d, COMPLETED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksNew(startTask+np-1,:) = results{np,1}; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = results{np,2};
                disp(sprintf('validPks(%d,:) becomes pksNew(%d,:)',startTask+np-1,startTask+np-1));
            end
        catch err 
            disp(sprintf('\n------------job%02d, FAILED with the following tasks------------',j));
            for np = 1:numel(job(j).Tasks)%fetchOutputs(job(j)) doesn't work so you have to read the number of tasks using numel(job(j).Tasks)
                pksNew(startTask+np-1,:) = NaN; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = NaN;
                disp(sprintf('pksNew(%d,:) = NaN',startTask+np-1));
            end
        end%end of try
    end
end

temp = isnan(pksNew(:,1));
npRedo = find(temp(:) == 1);

clear pks*

pksNew =[

          97         119          68
         216           1          66
          45         105         254
         153         113         221
         133         198         285
         270          62          68
         400          83          71
         493         250          70
         354         116         276
         368         208         250
         603           1          80
         493         249          70
         703         234          73
         434         149         209
         624           1         228
         848           1          82
         795         219          76
         795         219          76
         758         145         216
         862          47         240
         731         230         240
        1024           1          80
        1024         194         130
        1014         103         197
           1         341         157
         228         369          66
         240         479          72
         133         198         283
         130         501         267
         228         369          65
         493         249          70
         239         477          71
         344         500          73
         369         208         250
         243         502         244
         411         452         259
         493         249          70
         629         335          72
         704         235          73
         448         311         257
         589         513         253
         674         364         218
         796         219          76
         794         220          76
         818         455          67
         879         319          77
         731         231         242
         809         296         220
         674         364         218
         950         540         229
        1024         189         137
         970         495          73
         915         326         226
         950         541         230
         115         499          65
         240         478          71
         155         605          69
         130         501         268
           1         704         229
         182         638         234
         240         478          71
         344         501          74
         422         705           1
         243         502         244
         362         530         230
         362         530         230
         240         593         326
         362         530         230
         628         512          66
         626         627           1
         643         588           1
         589         513         253
         626         675         287
         851         513          70
         954         595          74
         820         730          65
         954         595          74
         793         584         215
         763         745         230
         954         595          74
         954         595          74
         950         541         227
          82         776          28
         210         743          64
           1         913           2
          83         908           7
          91         734         195
           1         808         222
         209         755          28
         497         829          66
         342        1024         112
         246         828         184
         228         945         253
         497         829          66
         624         607           6
         695         754          67
         749         866          66
         485         895         221
         624         693         209
         485         895         221
         819         731          65
         925         695          73
         840        1024          89
         762         746         234
         739         858         258
         847         976         248
         925         695          73
        1024         958          78
        1022         690         207
         944         957         255
          56        1024           4
         198        1024           1
         236         945         252
         176        1024           1
         338        1024         130
         227         945         250
         472        1024           1
         637        1024           1
         461         958         281
         545         972         232
         840        1024          90
         847         976         248
        1024         958          79
         943         957         255
          96         119          68
         215           1          66
          44         105         255
         153         114         220
         133         198         285
         270          62          68
         401          83          71
         493         249          70
         355         117         276
         369         207         250
         604           1          79
         495         248          71
         704         236          73
         452          89         291
         623           1         234
         848           1          82
         793         219          76
         795         219          77
         758         145         216
         862          47         240
         731         231         243
        1024           1          75
        1024         190         144
         968         156         254
           1         344          80
         228         369          65
         239         477          71
         133         198         286
          62         490         222
         228         369          65
         493         249          70
         239         477          71
         344         501          74
         369         208         250
         243         502         244
         411         452         260
         493         249          70
         630         334          72
         703         236          73
         448         312         257
         588         513         254
         674         364         218
         703         235          73
         795         219          76
         818         455          67
         879         318          77
         730         231         240
         809         295         218
         674         364         218
         950         541         230
        1024         191         103
        1024         518          68
         916         327         227
         949         541         230
         116         498          65
         240         478          71
         156         605          69
          61         490         223
           1         704         228
         182         638         234
         239         478          71
         344         501          74
         435         690           1
         243         502         244
         362         530         230
         362         530         230
         240         592         328
         362         530         230
         627         512          66
         605         606           1
         627         605           4
         588         513         254
         625         693         209
         851         513          70
         954         595          74
         819         731          65
         954         595          74
         793         584         215
         760         765         312
         954         595          74
         954         595          74
         950         541         228
           2         704           3
         221         850           7
           1         883           1
         119         858           1
          91         734         194
           1         808         224
         220         817           6
         497         829          66
         342        1023         152
         246         828         184
         227         944         251
         498         828          66
         600         632           6
         695         753          67
         750         866          66
         485         895         220
         625         694         208
         485         896         221
         818         730          66
         925         695          73
         840        1024          89
         759         764         315
         738         859         259
         847         976         248
         926         695          73
        1024         958          79
        1023         691         207
         943         957         255
          62        1024           1
         200        1024           1
         235         945         252
         196        1024           1
         339        1024         159
         227         944         251
         521        1024           1
         596        1024           3
         461         958         281
         545         972         232
         840        1024          93
         847         976         247
        1024         958          79
         944         957         254
          32          61         370
         152         114         220
         105         248         357
         105           1         401
         316          23         388
         369         208         250
         411         163         311
         339         137         350
         469           1         338
         688          99         381
         452          89         290
         706         261         333
         485          53         465
         571          14         471
         740          25         401
         645         100         467
         777          72         272
         789          95         356
         721         143         302
         850          87         459
         939         247         436
         967         156         255
         921         143         491
         105         249         356
         189         255         314
         130         501         269
         144         502         377
         139         370         493
         100         449         447
         369         208         250
         362         530         230
         467         506         387
         325         204         377
         390         366         321
         310         425         464
         435         281         329
         706         260         334
         467         507         387
         721         380         355
         467         505         388
         612         452         442
         747         221         326
         849         332         309
         805         408         338
         939         247         435
        1024         352         353
         964         431         390
         939         247         436
         131         500         269
         144         503         377
          51         720         350
         162         648         317
         144         502         375
         248         568         481
         128         698         475
         363         529         231
         467         506         388
         240         593         324
         248         568         481
         378         601         488
         467         506         388
         589         513         252
         471         603         427
         640         693         365
         704         684         318
         677         446         510
         624         695         474
         801         473         406
         729         544         453
         909         687         485
         937         510         328
        1013         616         268
         916         540         526
        1003         661         445
         102         822         367
         258         962         357
         128         782         533
         117         869         464
         184         920         527
         388         842         378
         258         961         357
         366         719         474
         339         855         444
         452         856         321
         704         684         317
         465         937         372
         689        1020         378
         625         695         473
         727         946         436
         954         821         356
         815         892         300
         835         789         502
         913         814         436
         866         875         374
         956         821         356
        1010         845         432
         942         862         514
         183        1009         380
         258         961         355
         145         969         464
         234         969         453
         258         961         357
         461         958         281
         598         974         330
         727         946         437
         872         956         337
         727         946         438
         943         957         255
         906         953         441
          39          18         467
         160         144         634
         158         143         638
          62         123         661
         159         143         637
         373         120         546
         255         159         602
         447         254         631
         254          23         664
         460         106         691
         338         129         660
         423         215         772
         484          52         465
         550         103         509
         741          85         484
         543         140         601
         565         116         693
         696          91         683
         832          10         560
         809         103         552
         873         222         580
         757          62         603
         884          74         731
         885         278         691
        1001          99         614
         949         213         644
        1002          98         599
         950         253         787
         139         371         491
         148         320         571
         222         329         487
          36         317         741
         200         334         741
         222         329         488
         433         367         600
         345         450         581
         282         393         733
         497         267         731
         359         407         663
         354         451         743
         454         429         658
         603         365         617
         730         369         469
         514         458         470
         572         456         547
         647         328         752
         728         239         689
         695         508         709
         873         222         579
         791         501         502
         867         482         576
         885         278         690
         927         296         586
         994         377         574
         950         213         642
        1024         437         624
         177         486         460
         248         568         482
         127         698         475
         224         703         589
         190         414         702
         100         720         714
         248         568         481
         276         570         620
         368         559         613
         315         539         707
         367         622         689
         442         690         684
         476         548         591
         618         643         615
         575         672         541
         608         501         767
         695         508         711
         528         674         649
         729         544         454
         848         682         566
         916         682         732
         915         541         526
         959         644         568
         981         542         706
         990         703         659
         128         782         533
         116         869         464
         184         920         527
         135         816         711
          28         885         718
         196         961         617
         366         719         476
         340         855         445
         390         675         607
         442         690         683
         263         972         543
         434         801         756
         624         695         474
         547         955         626
         724         941         563
         591         726         678
         549         906         699
         835         789         501
         914         814         436
         794         871         561
         816         764         605
         803         940         714
         942         862         513
         946         933         578
         990         703         660
         113         948         584
         145         969         464
         114         948         583
         382        1024         598
         263         972         543
         547         955         625
         724         941         562
         722         959         656
         727         945         436
         945         933         580
         803         940         714
         946         933         578
          88          82         833
         270          80         827
         480          28         752
         424         215         772
         628          98         762
         732         151         745
         801           1         715
         885          76         731
         762         238         766
         908         179         728
         993          76         813
         950         252         784
        1024         241         833
          37         316         740
         200         334         741
         283         392         734
         497         267         731
         354         451         745
         454         429         657
         569         386         833
         668         228         801
         523         465         746
         885         278         690
         769         498         766
         950         252         785
         934         346         754
         149         441         771
         100         721         712
         197         750         692
         264         618         708
         315         539         707
         408         476         816
         433         706         833
         608         502         768
         507         739         734
         564         681         772
         824         635         805
         917         683         733
        1024         641         733
         981         542         707
        1024         643         735
         100         721         712
         197         749         692
          27         884         717
         217         767         833
         434         801         756
         171         851         833
         434         801         756
         507         739         732
         618         754         770
         605         828         703
         890         804         771
         812         782         833
         803         940         714
         973         732         833
         893         940         833
         111         914         763
         152         941         833
         487         962         833
         549         906         699
         726         924         808
         802         989         819
         976        1024         711];

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('pksAdded1stPass_row_col_slice.txt','r');
temp = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksAdded(:,1) = temp(1,:)';
pksAdded(:,2) = temp(2,:)';
pksAdded(:,3) = temp(3,:)';

npAdded = (size(pksNew,1)+1):(size(pksAdded,1)+size(pksNew,1));
pksNew = vertcat(pksNew,pksAdded);

pksNew = double(pksNew);
for slice = 50:size(InormRed,3)
    npTemp1 = find(round(pksNew(:,3)) >= (slice-35));
    npTemp2 = find(round(pksNew(:,3)) <= (slice + 35));
    npp = intersect(npTemp1,npTemp2);
    
    clear label 
    disp(sprintf('---------------------------x-----------------------------'));
    for jj = 1:numel(npp)
        label{jj} = num2str(sprintf('%d',npp(jj)));
        disp(sprintf('\tpksNew(%d,:) = [%d, %d, %d]',npp(jj),pksNew(npp(jj),1),pksNew(npp(jj),2),pksNew(npp(jj),3)));
    end
    disp(sprintf('---------------------------x-----------------------------\n'));  
    
    figure(1)
    subplot(1,2,1)
    hold off
    simage(InormRed(:,:,slice));
    hold on
    plot(pksNew(npp,2),pksNew(npp,1),'k*','markersize',20)
    text(pksNew(npp,2),pksNew(npp,1),label,'fontsize',12)
    title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20);

    subplot(1,2,2)
    hold off
    simage(InormGreen(:,:,slice));
    hold on
    plot(pksNew(npp,2),pksNew(npp,1),'k*','markersize',20)
    text(pksNew(npp,2),pksNew(npp,1),label,'fontsize',12)
    title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)

    pause
end


[dist, pksOrder] = orderPeaks(pksNew);

% view multimers:
for j = 376:numel(npGood)
    np = npGood(j);
    npp = pksOrder(1:3,np);
    
%     %for 2mers:
%     npp = union(pksOrder(1:3,np),[np:(np+1)]);
    
%     npp = pksOrder(1:3,np);
    
%     %for 5mers:
%     npp = union(pksOrder(1:6,np),[np:(np+4)]);

%     %for 6mers
%     npp = union(pksOrder(1:6,np),[np:(np+5)]);
    
    clear label 
    disp(sprintf('---------------------------x-----------------------------'));
    for jj = 1:numel(npp)
        label{jj} = num2str(sprintf('%d',npp(jj)));
        disp(sprintf('\tpksNew(%d,:) = [%.03f, %.03f, %.03f], distance(%d,%d) = %.04f',npp(jj),pksNew(npp(jj),1),pksNew(npp(jj),2),pksNew(npp(jj),3),npp(jj),np,dist(npp(jj),np)));
    end
    disp(sprintf('---------------------------x-----------------------------\n'));
   
    
    startSlice = round(min(pksNew(npp,3))-40);
    stopSlice = round(max(pksNew(npp,3))+40);

    if(startSlice < 1)
        startSlice = 1;
    end

    if(stopSlice > size(InormRed,3))
        stopSlice = size(InormRed,3);
    end

    nppClose = npp;
    nppClose(nppClose == np) = [];
    
    for slice = startSlice:stopSlice
        figure(1)
        subplot(1,2,1)
        hold off
        simage(InormRed(:,:,slice));
        hold on
        plot(pksNew(npp,2),pksNew(npp,1),'k*','markersize',20)
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        text(pksNew(npp,2),pksNew(npp,1),label,'fontsize',12)
        title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNew(%d of %d,:) = [%.04f, %.04f, %.04f] \n %d particles plotted',np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3),numel(npp)),'fontsize',20)

        subplot(1,2,2)
        hold off
        simage(InormGreen(:,:,slice));
        hold on
        plot(pksNew(npp,2),pksNew(npp,1),'k*','markersize',20)
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        text(pksNew(npp,2),pksNew(npp,1),label,'fontsize',12)
        title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNew(%d of %d,:) = [%.04f, %.04f, %.04f] \npksNew(%d of %d) = [%.04f, %.04f, %.04f]',nppClose(1),size(pksNew,1),pksNew(nppClose(1),1),pksNew(nppClose(1),2),pksNew(nppClose(1),3),nppClose(2),size(pksNew,1),pksNew(nppClose(2),1),pksNew(nppClose(2),2),pksNew(nppClose(2),3)),'fontsize',20)

        pause
    end
end

% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('sinlge_pksNew1stPass_row_col_slice.txt','w');
% for np = 1:(npDimer(1)-1)
%     fprintf(fid,'%d %d %d %d\n',np,pksNew(np,1),pksNew(np,2),pksNew(np,3));
% end
% fclose(fid);

clear pks*

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('pksNew1stPass_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4, Inf]);
fclose(fid);

order = temp(1,:)';
pksNew(order,1) = temp(2,:)';
pksNew(order,2) = temp(3,:)';
pksNew(order,3) = temp(4,:)';

npDelete1 = find(isnan(pksNew(:,1)) == 1); 
npDelete2 = find(pksNew(:,3) <= 90); %get rid of all particles found automatically in particleTracking that have z-coord < 100

npDelete = union(npDelete1,npDelete2);
pksNew(npDelete,:) = [];

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers1stPass_row_col_slice.txt','r');
dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksDimers(:,1) = dimers(1,:)';
pksDimers(:,2) = dimers(2,:)';
pksDimers(:,3) = dimers(3,:)';

npDimer = (size(pksNew,1)+1):(size(pksDimers,1)+size(pksNew,1));
pksNew = vertcat(pksNew,pksDimers);

% ----------correct pksNew

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('sinlge_pksNew1stPass_row_col_slice.txt','r');
temp = fscanf(fid,'%d %d %d %d\n',[4,Inf]);
fclose(fid);

npGood = temp(1,:)';

npBad = 1:(npDimer(numel(1))-1);
npBad(npGood) = [];

pksNew(npGood,1) = temp(2,:)';
pksNew(npGood,2) = temp(3,:)';
pksNew(npGood,3) = temp(4,:)';

pksNew(npBad,:) = nan([numel(npBad),3]);

% % print all good particles:
% npDeleteTot = find(isnan(pksNew(:,1)) == 1); 
% pksNew(npDeleteTot,:) = [];
% 
% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('sinlge2_pksNew1stPass_row_col_slice.txt','w');
% for np = 1:(size(pksNew,1)-numel(npDimer))
%     fprintf(fid,'%d %d %d %d\n',np,pksNew(np,1),pksNew(np,2),pksNew(np,3));
% end
% fclose(fid);


% pksNew =[
% 
%           47         114         245
%          161         106         235
%          355         104         260
%          368         204         243
%          441         146         190
%          756         137         200
%          869          46         240
%          731         230         240
%         1002         188         132
%          994          89         181
%           15         342         129
%          407         440         248
%          459         318         244
%          601         501         248
%          686         355         210
%          806         289         215
%          947         540         229
%          939         316         217
%           21         692         215
%          195         634         227
%          356         526         240
%          240         593         326
%          621         665         283
%          810         574         200
%          763         745         230
%           30         811         210
%          346         996         130
%          255         825         170
%          219         936         248
%          487         888         216
%          623         688         180
%          738         845         256
%          847         976         243
%          956         957         246
%          461         958         281
%          556         965         220
%          455          75         275
%          634           1         235
%          957         155         254
%           50         488         215
%          773         762         310
%          842         986         133
%           40          63         375
%          105         248         357
%          119           5         395
%          311          34         380
%          422         165         315
%          469           1         329
%          686          92         370
%          478          50         463
%          576          12         455
%          748          23         401
%          643          94         459
%          772          58         262
%          796          93         356
%          720         137         292
%          852          83         457
%          939         244         420
%          927         138         493
%          202         260         314
%          148         498         372
%          145         375         488
%           83         444         452
%          467         506         387
%          388         355         310
%          310         425         464
%          435         257         328
%          717         368         355
%          612         452         442
%          852         324         309
%          818         417         322
%         1024         352         353
%          964         427         384
%           59         716         355
%          162         648         317
%          248         568         481
%          128         698         475
%          378         605         483
%          471         603         427
%          637         685         380
%          719         665         307
%          684         439         511
%          630         695         473
%          813         470         407
%          736         550         447
%          909         687         485
%         1005         612         268
%          916         540         526
%         1003         657         451
%          110         832         357
%          266         947         366
%          122         780         533
%          117         869         464
%          191         906         533
%          373         830         348
%          368         704         472
%          450         846         314
%          465         937         372
%          690         996         370
%          717         931         460
%          957         815         356
%          814         881         307
%          834         782         510
%          916         804         436
%          870         865         374
%         1008         844         438
%          942         862         514
%          175         994         375
%          134         966         458
%          234         969         453
%          613         948         340
%          872         956         337
%          906         953         453
%          162         135         634
%           63         108         653
%          370         108         546
%          264         154         593
%          447         254         631
%          254          23         664
%          462          98         693
%          338         129         660
%          424         204         768
%          550         103         512
%          741          85         484
%          543         140         601
%          565         116         693
%          696          91         683
%          832          10         560
%          815         108         560
%          864         211         572
%          884          74         731
%          884         272         690
%         1003          92         600
%          949         213         644
%          950         253         784
%          148         313         571
%          235         323         487
%           36         317         741
%          209         319         735
%          434         359         600
%          347         447         571
%          283         384         730
%          500         259         743
%          359         395         663
%          353         447         743
%          454         429         670
%          609         354         615
%          514         458         470
%          578         456         547
%          644         319         758
%          728         239         689
%          695         495         701
%          792         497         509
%          860         472         576
%          927         296         586
%          995         368         565
%          182         473         457
%          190         409         701
%          100         720         705
%          276         566         618
%          372         555         617
%          318         533         704
%          367         622         689
%          442         690         684
%          476         545         591
%          618         643         615
%          562         664         541
%          528         664         649
%          845         668         555
%          917         675         732
%          964         638         560
%          981         534         697
%          990         703         659
%          135         816         711
%           37         874         715
%          196         961         617
%          394         674         600
%          274         961         543
%          428         800         756
%          544         960         626
%          724         931         560
%          600         723         675
%          549         906         699
%          800         866         561
%          809         760         605
%          803         940         709
%          946         933         578
%          109         948         592
%           92          79         833
%          273          77         827
%          480          28         752
%          626          86         762
%          743         149         720
%          769         232         764
%          908         179         728
%          993          76         813
%          573         376         810
%          668         228         790
%          529         461         743
%          769         498         766
%          934         346         754
%          149         441         777
%          200         743         695
%          264         618         705
%          408         476         826
%          433         706         833
%          507         733         734
%          568         668         772
%          216         759         787
%          632         745         770
%          607         823         703
%          893         802         770
%          973         732         829
%          111         914         763
%          163         941         831
%          483         960         832
%          729         920         811
%          802         989         828
%          383         549         410
%          256         527         393
%          358         376         394
%          432         306         418
%          518         294         364
%          691         188         388
%          929          24         433
%          393           5         440
%           38         245         424
%          106         316         422
%           91         626         409
%          437         417         420
%          636         560         447
%          627         887         434
%          551         951         410
%          175          55         141
%           71          15         147
%          173         179         140
%          299         279         145
%          352         386         142
%          518         362         140
%          596         225         144
%          614          55         144
%          339         681         147
%          386         580         145
%          527         676         143
%          603         781         145
%         1006         788         138
%          982         954         138
%          204         906         141
%          429         902         140
%          571         890         145
%          128         341         329
%          243         695         345
%          379         969         341
%          185         901         395
%          266         819         376
%          334         749         374
%          417         693         365
%          541         434         379
%          534          41         390
%          922         338         399
%          958         911         374
%         1012          92         313
%          775         500         318
%          881          33         354
%          559          69         299
%          404         514         320
%          451         423         333
%          533         537         324
%           12         803         376
%          414          71         370
%          777         946         365
%          533        1015         322
%          700         897         340
%          962         652         355
%          616         788         443
%          424         771         436
%          436         868         431
%          527         135         414
%          294         628         405
%          106         542         450
%          139         226         442
%         1008         961         454
%         1007         563         477
%          776         996         495
%          413        1004         482
%          541         845         480
%          649         975         515
%         1022          99         495
%          964         407         470
%            7         663         520
%          205         148         507
%          454         143         485
%          614         369         495
%          858         387         483
%          874         291         504
%          723         839         495
%          269         715         488
%          401         241         550
%          502         242         542
%           44         890         543
%          240         819         516
%          508         766         543
%          556         564         530
%          705         534         540
%          675         626         533
%          318         532         532
%         1000         478         525
%          926         764         541
%          629         184         550
%          899        1017         536
%          441         471         538
%          825         575         559
%           66         216         590
%          427         159         603
%          607          19         590
%          907          61         600
%          992         549         595
%          683         816         595
%          322         740         599
%          569         246         615
%          853         921         618
%          165         237         638
%          157         480         643
%          624        1010         651
%          898         813         620
%          899        1002         635
%         1006         872         631
%            4         280         654
%          510         321         653
%          934           6         668
%          990         292         677
%          889         508         668
%          520         809         670
%          463         594         670
%          552         435         650
%          798         669         651
%          835         377         808
%          307          44         730
%          376         115         747
%          336         204         743
%          297         287         714
%          208          62         744
%          113          44         740
%           63         153         744
%           49         640         748
%           12          66         761
%          888         578         740
%          733         997         748
%          680           3         773
%          901         966         728
%         1021         339         805
%          505         300         827
%          311         331         807
%          213         376         820
%          566        1015         830
%          390          17         795
%          802          26         772
%          814         130         791
%          856         285         775
%          998         165         787
%          968         432         790
%          952         617         805
%          838         880         786
%          766         818         770
%          607         888         787
%          507         877         789
%          315         649         790
%           44         404         770
%          124         344         790
%           38         734         820
%          419         305         780
%          914         125         810
%         1001         903         817
%          164           6         826
%          680         833         820
%          489         531         820
%          215         557         832
%          102         603         825
%           82         766         133
%          108         714         180
%         1002         685         190
%          989         761         243
%          301         445         222
%          241         501         252
%          150         213         242
%          126         184         281
%          133         492         258
%          113         560         314
%          765         218         321
%          707         267         330
%          932         504         328
%          993         541         366
%          337         125         356
%          325         204         383
%          673         299         469
%          735         360         469
%          334         865         444
%          334         798         477
%          348         921         585
%          386         991         615
%          146         693         581
%          234         693         586
%          186         782         598
%          233         851         635
%          758          50         603
%          766           2         641
%          684         576         639
%          624         564         674
%          724         955         662
%          652         923         691
%          538         528         655
%          538         559         729
%          610         501         757
%          617         553         784
%          861         717         814
%          882         705         823
%          824         635         800
%          799         562         821
%         1023         536         824
%          952         520         824];


% -------(start: get all parameters K,D,w,offset for each particle)--------

clearvars -except Inorm* Kactual Dactual wActual offsetActual

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single2_pksNew1stPass_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
fclose(fid);

order = temp(1,:)';
pksNew(order,1) = temp(2,:)';
pksNew(order,2) = temp(3,:)';
pksNew(order,3) = temp(4,:)';

npDeleteTot = find(isnan(pksNew(:,1)) == 1); 
pksNew(npDeleteTot,:) = [];

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers1stPass_row_col_slice.txt','r');
dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksDimers(:,1) = dimers(1,:)';
pksDimers(:,2) = dimers(2,:)';
pksDimers(:,3) = dimers(3,:)';

npDimer = (size(pksNew,1)+1):(size(pksDimers,1)+size(pksNew,1));
pksNew = vertcat(pksNew,pksDimers);

numTry = 10; %# of tries to calculate param(np,:)

param = zeros([size(pksNew,1), 4]);
Rsquared = zeros([size(pksNew,1), 1],'single');

tic;
for np = 1:size(pksNew,1)
    [v,center] = getCutoutSectVertices(pksNew(np,:),0.75*Dactual,size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    %bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    %rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    [param(np,:), Rsquared(np)] = getParameters(avgProfile(startIdx:stopIdx),Dactual,wActual,numTry);
    
    if(~isnan(param(np,1)))
        disp(sprintf('SUCCESSFUL computation of parameters for pksNew(%d,:)',np));
    end
    
    if(isnan(param(np,1)))
        disp(sprintf('FAILED computation of parameters for pksNew(%d,:)',np));
    end
end
elapsedTime = toc; %16.5133669seconds for 285 particles

K = param(:,1);
D = param(:,2);
w = param(:,3);
offset = param(:,4);

% param =(1.0e+02)*[
% 
%    0.002399155795574   0.852902832031250   0.188171730041504   0.004636802077293
%    0.002699895203114   0.827876663208008   0.208574047088623   0.004477139413357
%    0.002513396739960   0.815522384643555   0.221820545196533   0.004821137189865
%    0.001834321767092   0.772924728393555   0.140175380706787   0.005159162878990
%    0.004339887201786   0.976313705444336   0.398254547119141   0.003454791307449
%    0.002063949704170   0.773710708618164   0.157455005645752   0.005129558444023
%    0.001747094541788   0.735041046142578   0.152868022918701   0.005624603033066
%    0.002339206337929   0.800012588500977   0.192875442504883   0.005311509966850
%    0.001644509136677   0.591372489929199   0.194154911041260   0.005335656404495
%    0.002832005918026   0.801290817260742   0.339936180114746   0.004448677003384
%    0.001711291670799   0.640976562500000   0.202699356079102   0.005150427222252
%    0.002110033333302   0.737530364990234   0.187679710388184   0.005089600086212
%    0.002558768689632   0.740640945434570   0.248750534057617   0.004887113273144
%    0.002119576632977   0.734271621704102   0.200057487487793   0.005126032829285
%    0.001911836862564   0.733770141601563   0.168457851409912   0.004924319982529
%    0.002684530615807   0.855536117553711   0.225781688690186   0.004478999972343
%    0.002020770609379   0.732295074462891   0.154375610351562   0.005499294400215
%    0.002664994001389   0.810934600830078   0.208932361602783   0.004205556511879
%    0.001902681589127   0.749789047241211   0.169145317077637   0.004904024899006
%    0.002912310361862   0.821709899902344   0.227163181304932   0.004593611657619
%    0.002716964185238   0.836206283569336   0.179299201965332   0.004821036458015
%    0.002050859332085   0.750887069702148   0.142740507125854   0.005230463147163
%    0.002053828537464   0.721777496337891   0.143964796066284   0.005369352102280
%    0.002708112597466   0.851037597656250   0.238950233459473   0.004271659851074
%    0.002306649833918   0.713111953735352   0.168099021911621   0.005127128362656
%    0.002054137587547   0.766894683837891   0.156467742919922   0.005026197433472
%    0.002117600589991   0.711743392944336   0.229026699066162   0.004807805120945
%    0.003237509131432   0.869193725585938   0.236697635650635   0.003555594384670
%    0.002486547827721   0.725330429077148   0.191789016723633   0.005009486675262
%    0.001912858784199   0.748124694824219   0.141952180862427   0.005543419122696
%    0.002326345890760   0.710005416870117   0.156199541091919   0.004882417619228
%    0.002072054594755   0.729846038818359   0.199922065734863   0.005125429630280
%    0.002726852595806   0.776755447387695   0.221079216003418   0.005118646025658
%    0.001876872479916   0.698100357055664   0.149284629821777   0.005717042088509
%    0.001856718808413   0.756748275756836   0.152344036102295   0.005339388251305
%    0.002265516370535   0.695430755615234   0.198785667419434   0.004919430315495
%    0.001405050009489   0.669663925170898   0.208562641143799   0.006319588422775
%    0.001962264776230   0.788521118164062   0.192760734558105   0.005109837055206
%    0.002204645574093   0.758237457275391   0.163722267150879   0.005165876746178
%    0.002192708402872   0.775413513183594   0.150029668807983   0.005219565033913
%    0.004848410785198   1.148626022338867   0.387155799865723   0.002323831021786
%    0.001885329186916   0.748891906738281   0.162431030273437   0.005690177679062
%    0.002277025431395   0.739242401123047   0.201624317169189   0.005350279211998
%    0.005301187634468   1.367215118408203   0.559642639160156   0.002039533108473
%    0.002581827044487   0.900925369262695   0.319296550750732   0.004767630398273
%    0.001903019100428   0.728562088012695   0.151852598190308   0.005262252092361
%    0.000975187644362   0.593207855224609   0.232037734985352   0.006661323904991
%    0.001846386492252   0.682949218750000   0.179519386291504   0.005550254583359
%    0.002182753533125   0.799698867797852   0.181519813537598   0.005357351303101
%    0.001884684115648   0.738108291625977   0.284587802886963   0.005536752939224
%    0.001344793289900   0.772734603881836   0.165622329711914   0.005781846642494
%    0.001914218366146   0.785532913208008   0.178040142059326   0.005475063920021
%    0.002406969070435   0.819778594970703   0.224235515594482   0.004851265549660
%    0.001990700364113   0.767936630249023   0.157584295272827   0.005507746338844
%    0.001859218031168   0.733214187622070   0.162952995300293   0.005548475384712
%    0.001809299886227   0.711666870117188   0.169266052246094   0.005588148236275
%    0.001854108124971   0.708142013549805   0.155819625854492   0.005426533222198
%    0.002159760743380   0.771605224609375   0.207688465118408   0.004993514418602
%    0.002704875171185   0.896424484252930   0.230220146179199   0.004367878735065
%    0.001915618032217   0.746324691772461   0.144796495437622   0.005658185482025
%    0.001621046066284   0.723254241943359   0.134689273834229   0.005684321522713
%    0.001836816966534   0.763336715698242   0.127180252075195   0.004881728291512
%    0.002538558244705   0.840088882446289   0.214350605010986   0.004841426014900
%    0.002583972215652   0.797050018310547   0.228114967346191   0.004727274775505
%    0.005140203237534   1.099832763671875   0.345637550354004   0.002459199428558
%    0.002672361135483   0.756466827392578   0.249688453674316   0.004588561654091
%    0.002210895270109   0.728952865600586   0.204640541076660   0.005299702286720
%    0.001518064737320   0.694083786010742   0.125068664550781   0.005740122795105
%    0.001947081685066   0.700905685424805   0.164821472167969   0.005218678712845
%    0.001614133566618   0.659304046630859   0.131490631103516   0.005367479920387
%    0.004812196791172   1.388297119140625   0.509119262695312   0.003045077919960
%    0.002352741211653   0.752821731567383   0.193698959350586   0.005198199748993
%    0.001897093057632   0.701995544433594   0.141393365859985   0.005235174298286
%    0.003193866908550   0.834001312255859   0.214131259918213   0.004470751285553
%    0.001829196810722   0.756867523193359   0.140361003875732   0.005725191831589
%    0.002250569611788   0.743613586425781   0.196859531402588   0.005199297666550
%    0.002064433097839   0.752961120605469   0.188178539276123   0.005328863859177
%    0.001644717454910   0.763934173583984   0.165890674591064   0.005437748432159
%    0.002272111773491   0.746573410034180   0.175078582763672   0.005454352498055
%    0.002599526643753   0.720845489501953   0.287084712982178   0.004617743194103
%    0.001587268263102   0.731637802124023   0.125180559158325   0.005484772920609
%    0.002764539718628   0.892342071533203   0.245151672363281   0.004831663072109
%    0.001784791797400   0.672115325927734   0.122768936157227   0.005299627184868
%    0.001760119050741   0.714375076293945   0.140515890121460   0.005357581973076
%    0.002088619470596   0.737830505371094   0.250051841735840   0.005510224699974
%    0.001311417818069   0.664466781616211   0.109498205184937   0.006034016013145
%    0.005269687771797   1.264479675292969   0.327985954284668   0.002274108380079
%    0.002128301262856   0.788474578857422   0.206680488586426   0.005495156645775
%    0.001971217095852   0.752721023559570   0.136845521926880   0.005300506353378
%    0.002162450253963   0.707376861572266   0.188166027069092   0.005428369045258
%    0.001903455108404   0.846447448730469   0.137851104736328   0.005597648620605
%    0.001674644201994   0.784506301879883   0.154307613372803   0.005517895221710
%    0.001797301620245   0.789207611083984   0.162578582763672   0.005461647510529
%    0.002400453090668   0.759789199829102   0.208927116394043   0.004399228394032
%    0.002242641448975   0.768721084594727   0.170476512908936   0.004988448321819
%    0.006999680995941   1.352394409179688   0.541523895263672   0.001020300462842
%    0.001906739622355   0.724737854003906   0.169026851654053   0.005614421963692
%    0.001959408074617   0.728236541748047   0.195403957366943   0.005173650383949
%    0.002355699986219   0.733626251220703   0.192397708892822   0.005258064270020
%    0.001876773536205   0.709136657714844   0.157455739974976   0.005673503279686
%    0.001962298750877   0.723394851684570   0.176996326446533   0.005589019060135
%    0.002061428874731   0.767487258911133   0.163947448730469   0.005528426170349
%    0.001925194710493   0.743201675415039   0.173249092102051   0.005701632499695
%    0.002286012470722   0.690790863037109   0.214681701660156   0.005439190268517
%    0.001916496306658   0.647575683593750   0.245525836944580   0.005714480280876
%    0.001928137838840   0.797372436523437   0.169994182586670   0.005679308772087
%    0.002198660373688   0.731158218383789   0.232807331085205   0.005354512333870
%    0.002079178243876   0.821771240234375   0.184691753387451   0.005288544893265
%    0.002922373116016   0.885136871337891   0.272032089233398   0.005278680324554
%    0.002050687968731   0.799799041748047   0.166426334381104   0.004460510611534
%    0.002922925949097   0.814297103881836   0.261519165039062   0.004977251589298
%    0.001986677795649   0.806185607910156   0.177275581359863   0.005484100580215
%    0.002290885448456   0.869483489990234   0.174292659759521   0.005189114809036
%    0.001521792709827   0.640970001220703   0.172802257537842   0.005601813197136
%    0.001441912055016   0.720282592773438   0.135819234848022   0.005590635538101
%    0.002182639688253   0.786178207397461   0.211058826446533   0.005109425187111
%    0.003345620930195   0.833331527709961   0.371011505126953   0.004504611790180
%    0.001133967414498   0.716469802856445   0.140933246612549   0.005969628691673
%    0.001589291244745   0.639809226989746   0.143297939300537   0.005602585673332
%    0.002542693912983   0.843278503417969   0.250056362152100   0.005009369254112
%    0.001425886601210   0.661991043090820   0.146538543701172   0.005736851096153
%    0.002279426157475   0.752091064453125   0.268955821990967   0.005390680432320
%    0.001772165447474   0.757567901611328   0.143982839584351   0.005716528892517
%    0.002177606374025   0.773919143676758   0.168795261383057   0.005563245415688
%    0.002506542205811   0.865020828247070   0.210013694763184   0.005204946398735
%    0.004897161424160   1.301836700439453   0.428854560852051   0.002792225480080
%    0.001121737360954   0.705978240966797   0.103687448501587   0.005980166792870
%    0.001748414784670   0.811486358642578   0.179731578826904   0.005478712916374
%    0.002274233847857   0.889659271240234   0.202166233062744   0.004751163721085
%    0.003487009704113   1.159634933471680   0.384177856445312   0.004225223362446
%    0.001192864254117   0.669517974853516   0.109100055694580   0.006420139074326
%    0.002078771293163   0.786946640014648   0.262237281799316   0.005562768578529
%    0.001685075014830   0.774405364990234   0.158293590545654   0.005758638381958
%    0.001539648175240   0.723678131103516   0.175440158843994   0.005969877243042
%    0.001987466514111   0.741544342041016   0.184753627777100   0.005511789321899
%    0.001871992051601   0.756637649536133   0.163608169555664   0.005088196396828
%    0.001974361389875   0.809013977050781   0.234358654022217   0.005453496575356
%    0.002597797214985   0.805160675048828   0.312958850860596   0.004719342291355
%    0.001907956302166   0.747864990234375   0.216576080322266   0.005427127480507
%    0.001952116936445   0.654330215454102   0.158718061447144   0.005274347662926
%    0.001385750621557   0.695482864379883   0.130480546951294   0.005937764644623
%    0.001489967852831   0.727862243652344   0.182378597259521   0.005996571779251
%    0.002238866984844   0.862671890258789   0.204958953857422   0.004886789917946
%    0.001215396150947   0.608744773864746   0.137219829559326   0.006123171448708
%    0.001817738115788   0.693330764770508   0.144463062286377   0.005896576642990
%    0.001718265265226   0.786983413696289   0.132352409362793   0.005411431789398
%    0.002092832624912   0.786269760131836   0.176482505798340   0.005521146059036
%    0.001957464069128   0.783438186645508   0.149920215606689   0.005614001154900
%    0.001700677871704   0.719710006713867   0.144664278030395   0.005778750777245
%    0.001342954039574   0.714648361206055   0.108157749176025   0.005705165266991
%    0.001760168075562   0.730490188598633   0.199453506469727   0.005379534959793
%    0.002975256741047   0.900896530151367   0.224727497100830   0.004254979491234
%    0.001723181307316   0.833144073486328   0.105890083312988   0.005201298594475
%    0.001659718602896   0.799931259155273   0.153713016510010   0.005964604616165
%    0.002854104638100   0.952137069702148   0.261803913116455   0.004613398313522
%    0.002094105929136   0.838207626342773   0.185561599731445   0.005104030966759
%    0.002078793793917   0.738973007202148   0.227171039581299   0.005408539772034
%    0.001728839427233   0.750025253295898   0.147013368606567   0.005580595135689
%    0.001686705648899   0.660696029663086   0.150828914642334   0.005855857729912
%    0.005705291032791   1.283733062744141   0.412009620666504   0.002188560515642
%    0.001759754270315   0.774657516479492   0.204484138488770   0.005697706937790
%    0.001585945338011   0.775795593261719   0.152391777038574   0.005799308419228
%    0.001445766836405   0.709722595214844   0.144542970657349   0.005789870619774
%    0.001321674287319   0.733368911743164   0.094848289489746   0.005858743786812
%    0.001996232718229   0.781024398803711   0.176732025146484   0.005210095047951
%    0.001686939299107   0.699524078369141   0.146019020080566   0.005690299868584
%    0.001862285137177   0.748355712890625   0.213214149475098   0.005492526888847
%    0.001656859368086   0.689699554443359   0.146922683715820   0.005711540579796
%    0.001596409678459   0.743970870971680   0.161814174652100   0.006083648800850
%    0.002855561077595   0.908047332763672   0.249184837341309   0.004704458117485
%    0.001859544962645   0.761653594970703   0.209201316833496   0.005765859484673
%    0.001488045752048   0.752424011230469   0.141064977645874   0.006088791489601
%    0.001923101842403   0.768267211914063   0.201367111206055   0.005523058772087
%    0.002283870279789   0.850175323486328   0.274153251647949   0.005122699141502
%    0.001947059184313   0.778550720214844   0.144434537887573   0.005466380119324
%    0.001383234560490   0.689563293457031   0.099564571380615   0.005807081460953
%    0.002386790663004   0.807865066528320   0.215204925537109   0.004988231360912
%    0.002068295329809   0.828012008666992   0.198069839477539   0.005462789535522
%    0.002629817724228   0.817054290771484   0.219646873474121   0.005042465925217
%    0.002281969934702   0.754136123657227   0.227376708984375   0.005444645285606
%    0.003517107963562   0.959750671386719   0.375519561767578   0.004263351261616
%    0.003065671920776   0.851871414184570   0.301124267578125   0.004926820397377
%    0.001878020167351   0.742850799560547   0.184634284973145   0.005625498294830
%    0.002529253661633   0.843559799194336   0.211300506591797   0.004769138693810
%    0.002524287104607   0.827814025878906   0.235195331573486   0.005318820476532
%    0.002187965661287   0.949174728393555   0.202419319152832   0.005226505994797
%    0.002168118655682   0.782188110351563   0.221133823394775   0.005107936859131
%    0.001372136324644   0.712898712158203   0.151272964477539   0.006343024373055
%    0.005933902859688   1.354922180175781   0.334792938232422   0.001740773469210
%    0.001195605173707   0.677215652465820   0.162268733978271   0.005899656414986
%    0.001648114621639   0.764463348388672   0.176473674774170   0.005647657513618
%    0.001199905127287   0.588506011962891   0.160334491729736   0.005734483599663
%    0.001929161250591   0.818841857910156   0.172611865997314   0.005583886504173
%    0.002253222614527   0.861528778076172   0.207680187225342   0.005337561964989
%    0.001174088567495   0.719706039428711   0.152408313751221   0.006229681968689
%    0.001790991425514   0.731771011352539   0.162601203918457   0.005781747698784
%    0.002374544739723   0.865602874755859   0.251217193603516   0.005048449039459
%    0.004509738087654   1.186035232543945   0.438493766784668   0.003386709094048
%    0.001724000573158   0.840342254638672   0.199953937530518   0.005445917844772
%    0.003361976444721   0.973890075683594   0.446709823608398   0.004609041512012
%    0.001519748419523   0.731218643188477   0.159354915618896   0.005908726453781
%    0.003209374845028   0.919884338378906   0.234185142517090   0.004571252465248
%    0.006823925375938   1.356649322509766   0.321403694152832   0.000657031089067
%    0.001492638885975   0.732007598876953   0.161913433074951   0.006077616214752
%    0.001807714402676   0.888740768432617   0.225183658599854   0.005859021544456
%    0.001633559316397   0.770017166137695   0.173909397125244   0.005782757401466
%    0.002432519644499   0.895694427490234   0.275720920562744   0.004764159619808
%    0.002211770266294   0.774246749877930   0.210384750366211   0.005671330690384
%    0.002262492626905   0.845641708374023   0.229219455718994   0.005029344558716
%    0.001881916522980   0.804489364624023   0.181052093505859   0.005368818044662
%    0.001847279965878   0.785155639648437   0.198571472167969   0.005410122275352
%    0.001637793928385   0.701531524658203   0.163441543579102   0.006042973995209
%    0.001431537717581   0.726471252441406   0.143964052200317   0.005865644216537
%    0.003035147786140   1.024355163574219   0.323041610717773   0.004375569820404
%    0.000958663374186   0.665099258422852   0.084783239364624   0.006084313392639
%    0.001280045509338   0.727485809326172   0.116500034332275   0.006152777075768
%    0.001359388977289   0.697558517456055   0.186647415161133   0.006456812620163
%    0.001929416805506   0.744260711669922   0.136641921997070   0.005394211411476
%    0.002556019425392   0.839663619995117   0.219248752593994   0.004774020016193
%    0.001904058754444   0.686859130859375   0.137443017959595   0.005342966914177
%    0.002388432025909   0.791474380493164   0.198754768371582   0.004880516231060
%    0.002295204699039   0.740650177001953   0.178284606933594   0.005002074241638
%    0.002760221660137   0.883936386108398   0.245310974121094   0.004430671334267
%    0.001482763886452   0.676083145141602   0.169156093597412   0.005700753927231
%    0.005537505149841   1.453038635253906   0.449084053039551   0.001940622925758
%    0.002089055329561   0.786831283569336   0.198242282867432   0.005164272785187
%    0.002049954831600   0.773386230468750   0.192629528045654   0.005337167978287
%    0.002023281604052   0.684960327148437   0.159837045669556   0.005193485021591
%    0.005006452798843   1.090516357421875   0.397425460815430   0.002779881954193
%    0.003408951461315   1.059174804687500   0.304338111877441   0.003863138556480
%    0.003763024806976   0.624922714233398   0.169728317260742   0.001308542937040
%    0.002496263831854   0.764753875732422   0.171162643432617   0.004993716478348
%    0.002409060895443   0.788797683715820   0.225107231140137   0.004286838471889
%    0.001907032430172   0.773743591308594   0.180354766845703   0.004680526554585
%    0.003763481676579   0.876801300048828   0.312540283203125   0.003424895107746
%    0.006472218632698   1.245659103393555   0.316524257659912   0.000748796686530
%    0.002484526932240   0.647615661621094   0.197497310638428   0.004441175162792
%    0.002318677902222   0.680341873168945   0.178785781860352   0.004528168439865
%    0.001965452581644   0.733414688110352   0.153433294296265   0.004526478350163
%    0.001373298615217   0.651773452758789   0.148413534164429   0.005159766077995
%    0.003974438905716   0.933150177001953   0.270165157318115   0.002979172766209
%    0.003144288659096   0.865784683227539   0.241701679229736   0.003753090500832
%    0.002987869381905   0.833986206054688   0.206262931823730   0.003866133689880
%    0.002790000140667   0.771882705688477   0.194768753051758   0.003936903476715
%    0.006625053286552   1.312572784423828   0.460691337585449   0.000949818864465
%    0.002425018996000   0.741629180908203   0.226651325225830   0.004805443286896
%    0.003045275211334   0.814701614379883   0.225617637634277   0.003946982324123
%    0.002667422294617   0.741512527465820   0.209121227264404   0.004423930048943
%    0.002257589101791   0.765710830688477   0.143796949386597   0.004081938564777
%    0.002341687679291   0.747829132080078   0.202084407806396   0.005069301724434
%    0.002187212109566   0.715945739746094   0.163935317993164   0.004835323095322
%    0.001585257351398   0.716608047485352   0.095920047760010   0.005513235926628
%    0.002089520692825   0.759936370849609   0.200063800811768   0.005108578205109
%    0.002540862262249   0.762660064697266   0.200694122314453   0.004460925459862
%    0.002403904348612   0.729814147949219   0.176127014160156   0.004629701375961
%    0.002529984414577   0.767191848754883   0.215460205078125   0.004401385188103
%    0.001750732064247   0.752673950195312   0.124574728012085   0.005513265728951
%    0.001885639876127   0.744115295410156   0.174471225738525   0.005879573822021
%    0.001699447184801   0.716886672973633   0.163328857421875   0.005146454572678
%    0.001763053536415   0.725673904418945   0.154617109298706   0.005676690936089
%    0.001497658938169   0.711598129272461   0.177964763641357   0.005869375467300
%    0.002033206075430   0.676338119506836   0.139472541809082   0.005008773207664
%    0.001923832595348   0.775416412353516   0.220050888061523   0.005360119938850
%    0.005785967111588   1.221675262451172   0.498089370727539   0.002166681587696
%    0.002117224633694   0.730390701293945   0.176712226867676   0.005484060645103
%    0.002050870060921   0.717255096435547   0.157643547058105   0.005279009938240
%    0.002377999275923   0.721220932006836   0.237550849914551   0.005087012648582
%    0.002036576867104   0.737556991577148   0.208373832702637   0.005492090582848
%    0.001414431482553   0.668167343139648   0.133546972274780   0.005712816119194
%    0.004915309250355   1.144874343872070   0.384694595336914   0.002703714668751
%    0.002838189005852   0.847455444335937   0.295263977050781   0.004821552634239
%    0.002206949740648   0.715609054565430   0.172367210388184   0.004784561991692
%    0.002226236611605   0.722875289916992   0.187367591857910   0.005335053801537
%    0.002790320217609   0.813449249267578   0.176257820129395   0.004585470855236
%    0.002109970301390   0.722828674316406   0.212856426239014   0.005298891067505
%    0.001888891458511   0.844702911376953   0.161231765747070   0.005416618585587
%    0.003054640889168   0.821763229370117   0.316394424438477   0.004539478719234
%    0.002155258506536   0.719949722290039   0.165649070739746   0.004922694563866
%    0.002775331735611   0.816731872558594   0.217353057861328   0.004536480009556
%    0.001908503025770   0.734573593139648   0.142734746932983   0.005324305295944
%    0.001861660331488   0.707695693969727   0.163556900024414   0.006014412045479
%    0.001542651057243   0.769871826171875   0.179113883972168   0.005840570926666
%    0.003443238437176   0.975120773315430   0.325404472351074   0.004184302985668
%    0.001556336283684   0.711073989868164   0.143864536285400   0.005355548262596
%    0.003021245598793   0.785602569580078   0.176368350982666   0.004513958096504
%    0.002342651486397   0.824602279663086   0.204196128845215   0.004922202229500
%    0.001572414040565   0.766477127075195   0.201204071044922   0.006139097213745
%    0.002865141928196   0.858105697631836   0.215746955871582   0.004914591014385
%    0.001655262857676   0.754964370727539   0.137858619689941   0.005272870063782
%    0.003317574858665   0.832974319458008   0.282690391540527   0.004104527831078
%    0.001655221730471   0.730835723876953   0.129226312637329   0.005489861369133
%    0.001501888334751   0.705823669433594   0.118443670272827   0.005496184825897
%    0.001745253503323   0.749375610351563   0.159931459426880   0.005186179876328
%    0.002833413183689   0.917103729248047   0.244452762603760   0.004538041055202
%    0.001814562976360   0.731273803710938   0.159843244552612   0.005404101014137
%    0.002169333100319   0.748608398437500   0.170282306671143   0.005170695781708
%    0.002378371655941   0.851425399780273   0.222569255828857   0.004656531810760
%    0.001924252510071   0.773948822021484   0.138558006286621   0.005209020376205
%    0.001729563623667   0.669891357421875   0.190460166931152   0.005375084877014
%    0.002670497298241   0.831127853393555   0.206176185607910   0.004545546472073
%    0.002599850594997   0.804112167358398   0.259652996063232   0.004586242139339
%    0.001977022439241   0.804528732299805   0.154396915435791   0.005126245617867
%    0.002198285162449   0.793072128295898   0.188170433044434   0.005242074131966
%    0.001774245053530   0.719280929565430   0.142155418395996   0.005410445332527
%    0.002150478214025   0.813845825195312   0.196425933837891   0.005203256011009
%    0.003138584792614   0.929234313964844   0.402746887207031   0.004339833557606
%    0.003435553610325   0.961252517700195   0.295714759826660   0.004150130152702
%    0.001850502938032   0.748863143920898   0.194873542785645   0.005135808587074
%    0.003251405060291   1.021809997558594   0.315628242492676   0.004576050639153
%    0.002401059269905   0.828313598632813   0.198826847076416   0.004956439733505
%    0.002299615591764   0.817791213989258   0.227214069366455   0.004661288857460
%    0.002396699339151   0.853402328491211   0.222364845275879   0.004554262757301
%    0.001572982966900   0.711357803344727   0.158125333786011   0.005641076564789
%    0.001111203208566   0.658942184448242   0.169923629760742   0.006264492273331
%    0.001591452360153   0.770391693115234   0.125347146987915   0.005699830055237
%    0.002129902243614   0.831130752563477   0.232228946685791   0.005696055293083
%    0.006801610589027   1.370974121093750   0.509335212707520   0.001084061637521
%    0.002682299911976   0.859430084228516   0.266433258056641   0.004534697830677
%    0.001476302444935   0.690260467529297   0.160835227966309   0.005384251475334
%    0.001633160263300   0.711151046752930   0.141865310668945   0.005915149450302
%    0.002783655524254   0.878065872192383   0.277172451019287   0.004544402956963
%    0.003536231517792   0.915965881347656   0.318607196807861   0.004053073227406
%    0.001542048901320   0.789835815429688   0.147705793380737   0.005423930883408
%    0.002958399355412   0.980300521850586   0.277327098846436   0.004790696799755
%    0.001671572625637   0.759680099487305   0.187838439941406   0.006136087179184
%    0.001932469904423   0.781555023193359   0.220888233184814   0.005722121596336
%    0.003150395452976   1.073191909790039   0.237876243591309   0.003951642513275
%    0.001964316517115   0.783867111206055   0.172754287719727   0.005444178581238
%    0.001688011884689   0.862825241088867   0.251033554077148   0.005880877971649
%    0.004348390102386   1.193459548950195   0.434842720031738   0.003456652760506
%    0.001482004821301   0.736992645263672   0.188162059783936   0.005881634354591
%    0.002684596478939   0.868831253051758   0.285498046875000   0.005042241215706
%    0.001674952208996   0.681926727294922   0.182930850982666   0.005705903768539
%    0.002085249274969   0.770729827880859   0.183519821166992   0.005712606906891
%    0.003130688667297   0.886683502197266   0.292692508697510   0.004368287324905
%    0.001351304650307   0.759714584350586   0.126982946395874   0.005890008211136
%    0.001689728647470   0.781366577148437   0.197519760131836   0.005742967128754
%    0.001981154233217   0.774453659057617   0.229484386444092   0.005258646607399
%    0.001530809551477   0.711581192016602   0.173983268737793   0.005328125357628
%    0.002667502462864   0.805622100830078   0.307303390502930   0.004485272467136
%    0.004595657885075   1.156373062133789   0.306610431671143   0.003054147958755
%    0.001751660257578   0.757230758666992   0.194198493957520   0.005534640550613
%    0.002666705846786   0.860817184448242   0.271520442962646   0.005112269520760
%    0.001451287120581   0.712890167236328   0.111643962860107   0.005439339876175
%    0.001581501960754   0.687333908081055   0.227140369415283   0.005787984728813
%    0.004297499358654   1.259583740234375   0.377326774597168   0.003619452714920
%    0.001131404116750   0.699374160766602   0.099446773529053   0.006406468153000
%    0.001436215937138   0.730950851440430   0.210014610290527   0.005906388759613
%    0.001546626985073   0.757558441162109   0.175382862091064   0.006026486754417
%    0.001070413142443   0.660922927856445   0.120144596099854   0.006558711528778
%    0.001358033865690   0.735707015991211   0.122331943511963   0.006109679937363
%    0.001516791731119   0.725734405517578   0.102719106674194   0.005852268934250
%    0.001757181286812   0.743041687011719   0.161681900024414   0.005865421295166
%    0.001577341556549   0.826918563842773   0.155433263778687   0.005330696702003
%    0.004090065658092   1.324483032226563   0.679926147460937   0.003458126783371
%    0.001851513087749   0.802531585693359   0.260002956390381   0.005667832493782
%    0.001484689116478   0.738187789916992   0.166085834503174   0.005758491754532
%    0.001803453117609   0.738980407714844   0.176348724365234   0.005875782370567
%    0.001366395801306   0.683758926391602   0.167425193786621   0.006228859424591
%    0.001285294592381   0.686601104736328   0.131614627838135   0.006514380574226
%    0.001474398970604   0.717035217285156   0.165671176910400   0.006264981031418
%    0.002603221535683   0.856117630004883   0.275933189392090   0.005273547768593
%    0.006409434080124   1.307674255371094   0.400382041931152   0.001262541860342
%    0.002394996881485   0.897696456909180   0.228619194030762   0.005322678089142
%    0.001722202301025   0.756017150878906   0.134340782165527   0.005969157218933
%    0.001689669936895   0.806633453369141   0.142598991394043   0.005459526181221
%    0.004571135044098   1.248732833862305   0.541756095886230   0.003115262985229
%    0.001337109655142   0.700169677734375   0.196944961547852   0.006015022397041
%    0.003956516683102   1.356957244873047   0.344356002807617   0.003923469185829
%    0.002911441624165   0.952355651855469   0.314488887786865   0.005140407085419
%    0.001601929366589   0.807462387084961   0.166953353881836   0.005988007187843
%    0.001883727610111   0.778503189086914   0.182665786743164   0.005935158133507
%    0.001048780903220   0.749705352783203   0.155408477783203   0.006326518058777
%    0.001673568189144   0.837528076171875   0.176486892700195   0.005555289983749
%    0.001861568093300   0.714581222534180   0.216781368255615   0.005653299093246
%    0.002846879661083   1.001155624389648   0.295301971435547   0.004465943276882
%    0.001816833466291   0.785151672363281   0.150845432281494   0.005398336648941
%    0.002752310037613   0.737258377075195   0.213875808715820   0.004130223393440
%    0.002052751332521   0.666280441284180   0.139309654235840   0.004700674116611
%    0.002613074779510   0.831892852783203   0.247707061767578   0.004474525749683
%    0.002171046584845   0.698366012573242   0.232350921630859   0.005392389297485
%    0.005053592324257   1.108731613159180   0.311181449890137   0.002320296615362
%    0.001587955504656   0.749605331420898   0.132744789123535   0.005806143283844
%    0.002385502308607   0.796761627197266   0.212878246307373   0.004674561619759
%    0.001577220857143   0.700052642822266   0.135017786026001   0.005554210543633
%    0.002254103422165   0.825198211669922   0.223241233825684   0.005402856469154
%    0.002384315580130   0.847536315917969   0.228363285064697   0.005013938546181
%    0.005482913255692   1.340830993652344   0.361918983459473   0.001860444396734
%    0.001760503798723   0.667162399291992   0.158011980056763   0.005223993062973
%    0.001623573005199   0.662318725585938   0.135490617752075   0.005489204525948
%    0.002238596230745   0.805527191162109   0.223960018157959   0.005134783983231
%    0.004863540232182   1.138839645385742   0.357170104980469   0.002897183895111
%    0.001576955616474   0.730380020141602   0.180856113433838   0.006158912181854
%    0.001782636642456   0.735988159179688   0.160094070434570   0.005467735528946
%    0.002236324697733   0.822190322875977   0.220435543060303   0.004950608611107
%    0.001699565500021   0.712417221069336   0.151588163375854   0.005185596346855
%    0.002099319547415   0.764851150512695   0.214069900512695   0.004968864321709
%    0.002165467590094   0.735257339477539   0.222224769592285   0.005156137347221
%    0.002390152364969   0.763685150146484   0.280283203125000   0.004745567739010
%    0.001511449664831   0.752103576660156   0.138678522109985   0.005183366537094
%    0.002862467467785   1.000421981811523   0.289650669097900   0.004120579659939
%    0.001979831457138   0.728271102905273   0.138630733489990   0.005479627251625
%    0.001614513397217   0.709712753295898   0.145726280212402   0.005542614459991
%    0.002069312334061   0.770458374023438   0.223784275054932   0.005547253489494
%    0.001785129904747   0.773252716064453   0.164297237396240   0.005266921520233
%    0.001217740625143   0.620577621459961   0.223407993316650   0.006000456213951
%    0.001102772802114   0.693284835815430   0.162398128509521   0.006367119550705
%    0.001595518141985   0.775068893432617   0.142764987945557   0.005944435000420
%    0.001294711232185   0.741084289550781   0.128252687454224   0.006135792732239
%    0.001382003873587   0.665839004516602   0.125552120208740   0.006138646006584
%    0.004080435633659   1.147267532348633   0.462971687316895   0.003795357048512
%    0.001151531338692   0.675037994384766   0.099849014282227   0.006176242232323
%    0.001724790334702   0.779597167968750   0.181343326568604   0.005944886803627
%    0.005223070383072   1.208223190307617   0.259674186706543   0.002052983343601
%    0.001063861027360   0.740437850952148   0.115399150848389   0.006436208486557
%    0.002135568708181   0.907623748779297   0.246829185485840   0.005545855760574
%    0.001355901211500   0.671952514648438   0.148086748123169   0.006635276675224
%    0.001790474951267   0.856012039184570   0.218916282653809   0.005663519501686
%    0.001165803447366   0.797534484863281   0.194498138427734   0.006000462770462
%    0.001325949132442   0.874669647216797   0.219196987152100   0.006480336785316
%    0.001883194744587   0.834361801147461   0.244827175140381   0.005860878825188];

% ----------------------(eliminate multiple peaks)--------------------

% criterion1: K must be POSITIVE (K>0)! Eliminate all pksNew with
% corresponding K values that are <=0 OR K == NaN (note: K == NaN means
% that nlinfit was not able to solve for parameters K, D, w, offset)
% 
% The brightest regions of a particle are towards it center...
% 
% K < 0 means that the sphere spread function(ssf) for the particle is 
% inverted so that the center of a particle is the darkest regions of the
% particle.  This is NOT possible in Inorm!
% 
% K = 0 means that the ssf of a particle is flat (i.e., the spherical
% profile is a horizontal line). This also is impossible! 
 
pksToEliminate = find(K(:)<=0 | isnan(K)); %there aren't any



% criterion2: Particle size is known (i.e., 5um +/-0.44um silica,
% manufactured by Bang Laboratories). All particles identified by 'pksNew'
% should have sizes ~Dactual +/- 10% of Dactual 
% 
% IMPORTANT NOTE: Dactual or D are really just parameters of ssf() NOT
% particle size.  You found empirically that particle size is usually some
% combination of D and w
% 
% But some particles are slightly less than 10% larger or smaller than
% Dactual. SO YOU MUST PLOT THE PARTICLE SIZES AS DONE BELOW, MANUALLY 
% LOOK AT THESE PARTICLES, AND ACCORDINGLY ADJUST lBound and uBound TO
% EITHER KEEP OR DELETE THEM. THIS MANUAL ADJUSTMENT MUST BE DONE WITH EACH
% PASS OF YOUR PARTICLE TRACKING ROUTINE!

lBound = 0.7*Dactual;
uBound = 1.3*Dactual;
lowerBound = lBound*ones([size(D,1),1],'single');
upperBound = uBound*ones([size(D,1),1],'single');

figure(1)
hold off
plot(D,'b-*')
hold on
plot(lowerBound,'r--')
plot(upperBound,'r--')
title(sprintf('D = 5.06 +/- 0.44 um according to Bang Labs',Dactual),'fontsize',20)
xlabel(sprintf('actually plotted:\n (%.03f*Dactual), r-- upper \n(%.03f*Dactual), r-- lower \n Dactual = %.02f pixels (size of ideal maskExp)',(uBound/Dactual),(lBound/Dactual),Dactual),'fontsize',20)

pksToEliminate = find(D>uBound | D<lBound); %49 particles

for np = 1:length(pksToEliminate)
    [v,center] = getCutoutSectVertices(pksNew(pksToEliminate(np),:),D(pksToEliminate(np)),size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    figure(2)
    hold off
    h(1) = plot(rNew,avgProfile(startIdx:stopIdx),'b*');
    hold on
    h(2) = plot(rNew,ssf([K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np))],rNew),'k-');
    legend(h,'InormRed','mask')
    title(sprintf('pksNew(%d of %d,:) = [%.04f, %.04f, %.04f] \n K = %.04f, D = %.04f, w = %.04f, offset = %.04f \n bad D because D < %.04f OR D > %.04f',pksToEliminate(np),size(pksNew,1),pksNew(pksToEliminate(np),1),pksNew(pksToEliminate(np),2),pksNew(pksToEliminate(np),3),K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),lowerBound(1),upperBound(1)),'fontsize',20)
    xlabel(sprintf('Rsquared = %.04f',Rsquared(pksToEliminate(np))),'fontsize',20)

    [voronoiVol, Rmask] = peakPlacementB(center,size(cutoutSect),D(pksToEliminate(np)),w(pksToEliminate(np)));
    mask = getCalcImg(Rmask,K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),voronoiVol);
    
    sphericalRegion = find(voronoiVol == 1);
    diff = cutoutSect;
    diff(sphericalRegion) = abs(cutoutSect(sphericalRegion)-mask(sphericalRegion));       
    
    b = 1:numel(mask); 
    b(sphericalRegion) = []; %indices to identify region of the image where there isn't a particle
    mask(b) = cutoutSect(b);
    
    figure(3)
    hold off
    simage([cutoutSect(:,:,round(center(3))), max(cutoutSect(:))*ones([size(cutoutSect,1) 1],'single'), mask(:,:,round(center(3))), max(cutoutSect(:))*ones([size(cutoutSect,1) 1],'single'), diff(:,:,round(center(3)))]);
    colorbar
    hold on
    plot(center(2),center(1),'k*','markersize',20)
    title(sprintf('InormRed \t mask \t abs(mask - InormRed)'),'fontsize',20)
    xlabel(sprintf('pksNew(%d of %d,:) = [%.04f, %.04f, %.04f] \n K = %.04f, D = %.04f, w = %.04f, offset = %.04f \n bad D because D < %.04f OR D > %.04f',pksToEliminate(np),size(pksNew,1),pksNew(pksToEliminate(np),1),pksNew(pksToEliminate(np),2),pksNew(pksToEliminate(np),3),K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),lowerBound(1),upperBound(1)),'fontsize',20)
    pause
end

% --------------------


% criterion 2B: look at particles with peculiar values of 'w' or more
% specifically, peculiar ratios of w to D.  Select those which are
% incorrect and delete them

ratio = w./D;

uBound = (wActual/Dactual)*1.35;
lBound = (wActual/Dactual)*0.65;

upperBound = uBound*ones([length(ratio) 1],'single');
lowerBound = lBound*ones([length(ratio) 1],'single');

figure(4)
plot(ratio,'b-*')
hold on
plot(upperBound,'r--')
plot(lowerBound,'r--')
title(sprintf('w/D \n wActual/Dactual = %.04f/%.04f = %.04f',wActual,Dactual,(wActual/Dactual)),'fontsize',20);
xlabel(sprintf('upperBound == (wActual/Dactual)*%.03f, r-- \n lowerBound == (wActual/Dactual)*%.03f, r--',uBound*(Dactual/wActual),lBound*(Dactual/wActual)),'fontsize',20);

pksToEliminate = find(ratio > uBound | ratio < lBound); %230. Too many to consider! But keep in mind the backgrd is very bright/yellow so the halo will be very large on most particles


% criterion3: Now identfity multiple peaks locating the same particle. 
% (i.e., peaks or centers of maskExp's separated by less than Dactaul) 
% This can happen because of Monte Carlo minimization of Rsquared. 
% The random movement of maskExp (i.e., peak location) may move a peak on 
% the edge of one particle in Inorm onto a neighboring particle. 
% Keep the particle with the lowest value of Rsquared. (note: you cannot
% just use getNeighborsValid.m for this because this doesn't compare
% Rsquared values of each peak it just keeps the first peak and eliminates
% all subsequent peaks/mask images identifying the same particle)
[dist, pksOrder] = orderPeaks(pksNew);

dummyVal = 0; %this should be orders of magnitude larger than any real Rsquared value
RsquaredCompare = dummyVal*ones([size(pksNew,1), size(pksNew,1)],'single');

% How to read RsquaredCompare...
% 
% provided RsquaredCompare ~= dummyVal
% 
% RsquaredCompare(1,np) = value of Rsquared for particle 1 NOT np
% RsquaredCompare(2,np) = value of Rsquared for particle 2 NOT np
% ...
% RsquaredCompare(np,np) = value of Rsquared for particle np
tic;
for np = 1:size(pksNew,1)
    idxClose = find((dist(:,np) < Dactual*0.4)); %positions in pksNew that are too close to pksNew(np,:) 
    disp(sprintf('-------------\nfinding peaks too close to pksNew(%d,:)...',np))
    for npClose = 1:length(idxClose)%goes through positions too close to pksNew(np,:) and computes Rsquared
        disp(sprintf('pksNew(%d,:) is too close',idxClose(npClose)))
        [v, peakNew] = getCutoutSectVertices(pksNew(idxClose(npClose),:),2*D(idxClose(npClose)),size(InormRed)); %vertices of rectangular portion of 'Inorm' that is (2*Dactual)^3 and centered on validPks(np,:) 
        [naught,RsquaredCompare(idxClose(npClose),np)] = removeParticle(peakNew,K(idxClose(npClose)),D(idxClose(npClose)),w(idxClose(npClose)),offset(idxClose(npClose)),InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3)),idxClose(npClose));
    end
    disp(sprintf('-------------\n'))
end 
elapsedTime = toc; %6.351172820000000e+02 = 10.709min for 1054 particles

RsquaredCompare1 = RsquaredCompare;
for np = 1:size(pksNew,1)
    idxClose = find((dist(:,np) < Dactual*0.4)); %positions in pksNew that are too close to pksNew(np,:) 
    if(numel(idxClose)>1)%>1 because dist(np,np) == 0 so idxClose will always be at least 1 
        disp(sprintf('particle %d has %d peaks: ',np,numel(idxClose)))
        disp(idxClose)
        minRsquared = min(RsquaredCompare(idxClose,np));
        numMinRsquared = find(RsquaredCompare(idxClose,np) == minRsquared);
        if (length(numMinRsquared) > 1)%in the event that 2 or more peaks are EXACTLY the same or have exactly the same minimum Rsquared value then just select one peak delete all others
            idxClose(numMinRsquared(1)) = [];%one of the multipe peaks having the same min. value of Rsquared is deleted from idxClose
            npDelete = 1:numel(idxClose);%ALL other peak in idxClose are tagged for deletion
        end
        if (length(numMinRsquared) == 1) %only one peak has the minimum Rsquared value
            npDelete = find(RsquaredCompare(idxClose,np) > min(RsquaredCompare(idxClose,np)));
        end
        RsquaredCompare(idxClose(npDelete),np) = NaN;
    end
end


pksToEliminate = find(isnan(diag(RsquaredCompare)));

pksToEliminate(pksToEliminate >= npDimer(1)) = []; %There aren't any!

% ===left off here


% --(start: move peaks to minimize error and locate at pixel resolution)--

clearvars -except K D w offset pksNew Inorm* Kactual Dactual wActual offsetActual 

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

numTasks = 6; %one particle == one task
tasksPerJob = 1:numTasks:size(pksNew,1);
for j = 1:length(tasksPerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue
    disp(sprintf('\n------------job%02d, created with the following tasks------------',j));

    startTask = tasksPerJob(j);
    stopTask = startTask+numTasks-1; %subtract to account for MATLAB indexing displacement
   
    if(stopTask > size(pksNew,1)) 
        stopTask = size(pksNew,1);
    end
    
    for np = startTask:stopTask %slice == task #
        [v,center] = getCutoutSectVertices(pksNew(np,:),0.4*D(np),size(InormRed)); %vertices of rectangular portion of 'InormRed' that is (2*Dactual)^3 and centered on pksNew(np,:) 
        createTask(job(j),@getPosition,2,{v(1,:),center,K(np),D(np),w(np),offset(np),InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3))});
        disp(sprintf('locating pksNew(%d,:) at pixel resolution...',np));
    end
end

% submit one job at a time to be computed by all available slave nodes
for j = 1:length(job)
    submit(job(j));
    disp(sprintf('job%02d, submitted',j))
end 
% Actual time for 0.6*D(np):
%   first job started at: Mon Aug 03 19:36:33 EDT 2015
%   last job finished at: Mon Aug 03 22:08:14 EDT 2015
%   total time: 2hr31min41sec for locating 420 particles


% Actual time for 0.4*D(np):
%   first job started at: Tue Aug 04 15:11:51 EDT 2015
%   last job finished at: Tue Aug 04 16:21:28 EDT 2015
%   total time: 1hr09min37sec for locating 421 particles

pksNewer = zeros(size(pksNew),'single');
Rsquared = zeros([8,size(pksNew,1)],'single'); 
for j = 1:length(job)
    state = job(j).State;
    
    startTask = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        try 
            results = fetchOutputs(job(j));
            disp(sprintf('\n------------job%02d, COMPLETED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksNewer(startTask+np-1,:) = results{np,1}; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = results{np,2};
                disp(sprintf('pksNew(%d,:) becomes pksNewer(%d,:)',startTask+np-1,startTask+np-1));
            end
        catch err 
            disp(sprintf('\n------------job%02d, FAILED with the following tasks------------',j));
            for np = 1:numel(job(j).Tasks)%fetchOutputs(job(j)) doesn't work so you have to read the number of tasks using numel(job(j).Tasks)
                pksNewer(startTask+np-1,:) = NaN; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = NaN;
                disp(sprintf('pksNewer(%d,:) = NaN',startTask+np-1));
            end
        end%end of try
    end
end

temp = isnan(pksNewer(:,1));
npRedo = find(temp(:) == 1);

% pksNewer =[
% 
%           53         116         244
%          157         110         229
%          353         112         264
%          370         205         246
%          437         150         195
%          754         141         208
%          865          52         239
%          728         233         241
%         1002         192         140
%          998          99         188
%           14         343         136
%          409         447         253
%          454         316         248
%          596         510         247
%          681         358         212
%          806         294         214
%          945         539         230
%          931         322         220
%           18         698         219
%          186         638         228
%          360         531         237
%          240         593         325
%          622         673         284
%          803         578         203
%          760         749         228
%           27         810         211
%          343        1002         140
%          250         829         172
%          220         940         249
%          486         892         219
%          623         691         184
%          736         852         256
%          848         977         243
%          458         958         285
%          550         971         226
%          454          82         279
%          621           1         229
%          963         157         255
%           54         492         220
%          764         763         309
%          837         991         145
%           38          59         370
%          105         247         357
%          104          20         398
%          320          30         382
%          417         167         313
%          473           1         348
%          683         100         376
%          485          51         463
%          568          22         462
%          742          31         405
%          646         100         456
%          774          69         268
%          791          94         352
%          721         142         291
%          849          87         460
%          940         247         426
%          923         142         490
%          194         260         313
%          147         502         374
%          141         375         485
%           88         446         451
%          467         507         386
%          386         362         315
%          309         425         463
%          434         270         328
%          719         376         353
%          610         456         443
%          851         329         309
%          811         416         327
%         1024         356         353
%          963         429         390
%           57         720         353
%          162         647         315
%          246         567         477
%          127         694         477
%          376         603         483
%          472         606         426
%          635         692         380
%          711         673         311
%          682         446         506
%          627         696         469
%          808         474         408
%          730         549         449
%          908         684         484
%         1005         615         269
%          918         539         525
%          999         659         451
%          104         828         360
%          261         957         368
%          125         783         528
%          115         871         462
%          187         918         528
%          377         838         359
%          367         711         473
%          454         851         317
%          462         940         369
%          689        1003         373
%          717         939         460
%          953         819         357
%          816         888         301
%          831         786         511
%          915         811         436
%          867         871         370
%         1003         846         435
%          942         861         511
%          177         997         382
%          144         972         460
%          235         967         453
%          605         963         336
%          872         958         333
%          905         953         447
%          159         143         632
%           60         116         656
%          369         118         547
%          257         159         597
%          448         254         628
%          254          29         667
%          459         102         691
%          338         128         660
%          424         212         769
%          549         108         513
%          740          85         482
%          544         138         601
%          562         117         694
%          695          94         684
%          830          14         562
%          810         108         550
%          869         220         573
%          887          73         735
%          883         278         689
%         1001         100         595
%          946         215         647
%          951         252         782
%          146         320         574
%          229         328         484
%           38         315         739
%          201         329         739
%          429         362         602
%          347         446         572
%          284         392         735
%          498         270         743
%          359         403         662
%          352         452         745
%          451         426         674
%          604         361         614
%          517         459         468
%          573         458         549
%          644         328         760
%          730         242         688
%          694         503         703
%          791         499         505
%          865         479         576
%          923         294         590
%          990         379         564
%          179         485         457
%          187         415         702
%           99         722         708
%          280         568         619
%          367         556         620
%          319         542         701
%          364         618         689
%          443         689         688
%          478         549         594
%          621         644         611
%          567         671         543
%          528         671         644
%          847         678         560
%          914         685         735
%          963         645         557
%          979         542         693
%          989         706         664
%          131         818         716
%           33         882         717
%          199         962         614
%          390         675         602
%          265         968         538
%          429         802         766
%          544         957         624
%          722         939         560
%          594         728         680
%          546         907         703
%          797         871         556
%          810         765         602
%          808         938         709
%          946         935         576
%          108         947         587
%           87          77         833
%          271          75         804
%          480          31         751
%          625          97         760
%          736         151         727
%          763         238         766
%          909         179         729
%          989          72         815
%          574         383         806
%          669         229         794
%          527         462         738
%          773         496         771
%          937         347         756
%          147         443         777
%          199         750         695
%          266         617         706
%          407         473         818
%          441         696         833
%          506         741         738
%          566         679         772
%          213         765         787
%          623         750         765
%          607         829         705
%          889         805         768
%          977         731         824
%          111         918         766
%          152         940         826
%          484         959         833
%          726         926         812
%          801         988         827
%          382         550         408
%          250         538         396
%          358         376         391
%          429         307         413
%          517         297         361
%          691         193         386
%          923          26         433
%          395          26         433
%           42         249         420
%          108         326         418
%           90         634         408
%          438         427         420
%          632         568         454
%          622         890         433
%          546         959         413
%          173          64         147
%           64          20         154
%          169         187         148
%          297         282         151
%          348         392         144
%          515         368         144
%          593         228         151
%          610          54         154
%          334         681         149
%          383         578         156
%          528         680         147
%          601         782         147
%         1005         787         145
%          974         961         147
%          201         913         145
%          427         908         146
%          567         892         152
%          127         345         335
%          241         695         345
%          377         969         338
%          185         905         400
%          264         824         375
%          334         749         376
%          414         697         369
%          542         442         383
%          532          49         389
%          918         343         403
%          954         915         374
%         1012          99         312
%          774         501         318
%          880          45         353
%          553          68         304
%          404         523         318
%          445         432         336
%          530         537         321
%           10         803         376
%          415          66         368
%          778         946         365
%          533        1007         316
%          703         897         336
%          962         655         350
%          614         788         444
%          421         772         434
%          427         874         436
%          525         136         419
%          287         628         407
%          110         541         446
%          137         232         442
%         1009         958         457
%         1003         563         480
%          769        1002         491
%          410        1009         484
%          533         845         479
%          649         980         514
%         1014         107         499
%          964         409         470
%           10         666         517
%          205         155         512
%          452         148         487
%          609         377         498
%          858         387         481
%          880         294         508
%          721         842         498
%          263         718         492
%          403         248         552
%          497         249         544
%           46         896         543
%          229         825         519
%          507         778         544
%          557         570         529
%          701         541         539
%          667         628         534
%          318         545         536
%          991         488         525
%          921         771         540
%          625         190         553
%          899        1017         539
%          442         474         536
%          824         575         559
%           58         221         588
%          423         161         602
%          598          22         590
%          907          64         603
%          983         552         596
%          685         824         596
%          319         745         593
%          563         244         615
%          851         922         624
%          166         247         634
%          163         482         640
%          626        1012         650
%          895         818         618
%          901        1004         634
%         1009         874         629
%           12         288         655
%          505         326         653
%          923          13         672
%          976         297         682
%          886         508         660
%          514         815         667
%          464         592         669
%          552         439         648
%          801         675         649
%          832         377         814
%          301          46         735
%          373         121         742
%          333         207         743
%          300         298         722
%          210          63         748
%          112          50         748
%           63         160         743
%           51         647         752
%           12          71         760
%          889         579         740
%          732         999         745
%          681           9         774
%          894         973         728
%         1017         350         803
%          499         301         823
%          309         335         806
%          210         377         815
%          564        1016         827
%          388          17         795
%          808          32         781
%          814         129         794
%          857         289         775
%          992         172         782
%          965         434         786
%          953         623         809
%          843         891         783
%          763         820         770
%          605         892         782
%          508         881         792
%          311         649         788
%           50         405         781
%          122         349         790
%           29         735         811
%          413         310         781
%          910         130         809
%         1004         907         813
%          163          17         822
%          674         840         817
%          490         535         814
%          213         562         832
%          102         610         822
%           82         769         139
%          104         721         182
%         1002         688         198
%          992         764         241
%          916         900         198
%          950         960         239
%          297         455         221
%          239         501         249
%          146         215         247
%          131         195         280
%          132         496         261
%          115         555         314
%          757         222         321
%          708         262         327
%          939         513         332
%          984         541         362
%          337         130         351
%          323         206         380
%          672         304         476
%          733         366         472
%          336         862         443
%          339         798         476
%          350         924         586
%          378         996         609
%          148         695         581
%          231         701         584
%          184         788         604
%          226         850         632
%          759          51         607
%          768           1         639
%          679         574         640
%          627         572         682
%          721         960         660
%          642         927         694
%          532         536         661
%          531         554         727
%          609         508         759
%          614         556         781
%          872         713         815
%          882         705         821
%          824         632         808
%          791         568         821
%         1023         540         829
%          953         528         823];

for np = 374:size(pksNewer,1)
    startSlice = round(nanmin([pksNew(np,3),pksNewer(np,3)])-4);
    stopSlice = round(nanmax([pksNew(np,3),pksNewer(np,3)])+4);

    if(startSlice < 1)
        startSlice = 1;
    end

    if(stopSlice > size(InormRed,3))
        stopSlice = size(InormRed,3);
    end
    
    for slice = startSlice:stopSlice
        figure(1)
        subplot(1,2,1)
        hold off
        simage(InormRed(:,:,slice));
        hold on
        plot(pksNew(np,2),pksNew(np,1),'k*','markersize',20)
        plot(pksNewer(np,2),pksNewer(np,1),'g*','markersize',20)
        title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNewer(%d of %d,:) = [%d, %d, %d](g*)',np,size(pksNewer,1),pksNewer(np,1),pksNewer(np,2),pksNewer(np,3)),'fontsize',20)

        subplot(1,2,2)
        hold off
        simage(InormGreen(:,:,slice));
        hold on
        plot(pksNew(np,2),pksNew(np,1),'k*','markersize',20)
        plot(pksNewer(np,2),pksNewer(np,1),'g*','markersize',20)
        title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNew(%d of %d,:) = [%d, %d, %d](k*)',np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)
        pause
    end
end

% correct pksNewer....
cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('npKeep_1stPass.txt','r');
npKeep = fscanf(fid,'%d\n',[1,Inf]);
fclose(fid);

pksNewer(npKeep,:) = pksNew(npKeep,:);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('pksCorrect1stPass_id_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
fclose(fid);

npCorrect = temp(1,:)';

pksNewer(npCorrect,1) = temp(2,:)';
pksNewer(npCorrect,2) = temp(3,:)';
pksNewer(npCorrect,3) = temp(4,:)';

% % print pksNewer....

% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('single_pksNewer1stPass_id_row_col_slice.txt','w');
% for np = 1:size(pksNewer,1)
%     fprintf(fid,'%d %d %d %d\n',np,pksNewer(np,1),pksNewer(np,2),pksNewer(np,3));
% end
% fclose(fid);


% -------(start: get all parameters K,D,w,offset for each particle)--------

clearvars -except pksNew Inorm* Kactual Dactual wActual offsetActual

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single_pksNewer1stPass_id_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
fclose(fid);

order = temp(1,:)';
pksNewer(order,1) = temp(2,:)';
pksNewer(order,2) = temp(3,:)';
pksNewer(order,3) = temp(4,:)';

npDeleteTot = find(isnan(pksNewer(:,1)) == 1); 
pksNewer(npDeleteTot,:) = [];

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers1stPass_row_col_slice.txt','r');
dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksDimers(:,1) = dimers(1,:)';
pksDimers(:,2) = dimers(2,:)';
pksDimers(:,3) = dimers(3,:)';

npDimer = (size(pksNewer,1)+1):(size(pksDimers,1)+size(pksNewer,1));
pksNewer = vertcat(pksNewer,pksDimers);

numTry = 10; %# of tries to calculate param(np,:)

param = zeros([size(pksNewer,1), 4]);
Rsquared = zeros([size(pksNewer,1), 1],'single');

tic;
for np = 1:size(pksNewer,1)
    [v,center] = getCutoutSectVertices(pksNewer(np,:),0.75*Dactual,size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    %bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    %rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    [param(np,:), Rsquared(np)] = getParameters(avgProfile(startIdx:stopIdx),Dactual,wActual,numTry);
    
    if(~isnan(param(np,1)))
        disp(sprintf('SUCCESSFUL computation of parameters for pksNewer(%d,:)',np));
    end
    
    if(isnan(param(np,1)))
        disp(sprintf('FAILED computation of parameters for pksNewer(%d,:)',np));
    end
end
elapsedTime = toc; %16.5133669seconds for 285 particles

K = param(:,1);
D = param(:,2);
w = param(:,3);
offset = param(:,4);

% param =(1.0e+02)*[
% 
%    0.003562898039818   0.918083572387695   0.226605930328369   0.003652048408985
%    0.001932534128428   0.720434875488281   0.146487293243408   0.005199065804482
%    0.002513396739960   0.815522384643555   0.221820545196533   0.004821137189865
%    0.002293787896633   0.841781997680664   0.185515747070312   0.004764998853207
%    0.001623463928699   0.721919326782227   0.129335927963257   0.005572109818459
%    0.002659511864185   0.847238922119141   0.201189632415771   0.004720102548599
%    0.001747094541788   0.735041046142578   0.152868022918701   0.005624603033066
%    0.002391870319843   0.830698776245117   0.157755823135376   0.005190610289574
%    0.005188450217247   1.153372879028320   0.457786293029785   0.002339069545269
%    0.002446760982275   0.865815048217773   0.251356868743897   0.004736838042736
%    0.001469092816114   0.641076431274414   0.136634769439697   0.005320386290550
%    0.003846722245216   0.902932662963867   0.370740776062012   0.003926971256733
%    0.002261549532413   0.744716949462891   0.193804721832275   0.005098503232002
%    0.003691112697124   0.836626052856445   0.310403671264648   0.004089068770409
%    0.001772504746914   0.731140594482422   0.143985881805420   0.005077682137489
%    0.001784213930368   0.712766571044922   0.135263166427612   0.005327902436256
%    0.002154643386602   0.735811004638672   0.163566417694092   0.005442760586739
%    0.002664994001389   0.810934600830078   0.208932361602783   0.004205556511879
%    0.001739692389965   0.703044586181641   0.161885871887207   0.005167341232300
%    0.002039238512516   0.741999816894531   0.121283092498779   0.005380040407181
%    0.002178759574890   0.761250686645508   0.142766313552856   0.005357314348221
%    0.002174675315619   0.771346130371094   0.150843734741211   0.005112498998642
%    0.002053828537464   0.721777496337891   0.143964796066284   0.005369352102280
%    0.003874946534634   0.948159713745117   0.321849441528320   0.003412044644356
%    0.004430649876595   0.873547821044922   0.312487277984619   0.003516063094139
%    0.002072436660528   0.745432281494141   0.160976543426514   0.005071168541908
%    0.003689958453178   0.903697662353516   0.341402168273926   0.003591495156288
%    0.006280149817467   1.170356063842773   0.307123107910156   0.000648962408304
%    0.002619991898537   0.718437652587891   0.227249126434326   0.004987476766109
%    0.002343272119761   0.780860137939453   0.172736473083496   0.005248964428902
%    0.002216872274876   0.691464309692383   0.149291391372681   0.005042942166328
%    0.002072054594755   0.729846038818359   0.199922065734863   0.005125429630280
%    0.002648636996746   0.778019561767578   0.202569084167480   0.005134487152100
%    0.002312978357077   0.731424026489258   0.169318981170654   0.005361485481262
%    0.001772082000971   0.684269256591797   0.137843952178955   0.005621888041496
%    0.002860338091850   0.741007537841797   0.267124996185303   0.004529403746128
%    0.001615057140589   0.684988708496094   0.224055309295654   0.006158733367920
%    0.001962264776230   0.788521118164062   0.192760734558105   0.005109837055206
%    0.002652666270733   0.794259033203125   0.200819473266602   0.004880713522434
%    0.002367719262838   0.737227401733398   0.171760787963867   0.005287778973579
%    0.002332021445036   0.720322265625000   0.200221271514893   0.004849617481232
%    0.002144317179918   0.743426208496094   0.197162837982178   0.005602274537086
%    0.002102290391922   0.724618377685547   0.179487915039063   0.005469843149185
%    0.005301187634468   1.367215118408203   0.559642639160156   0.002039533108473
%    0.002581827044487   0.900925369262695   0.319296550750732   0.004767630398273
%    0.002287541329861   0.727333374023437   0.205321750640869   0.005054565072060
%    0.000975187644362   0.593207855224609   0.232037734985352   0.006661323904991
%    0.001846386492252   0.682949218750000   0.179519386291504   0.005550254583359
%    0.002509103417397   0.806169204711914   0.239637508392334   0.005201125741005
%    0.001687073558569   0.736679077148438   0.240755805969238   0.005664960145950
%    0.001344793289900   0.772734603881836   0.165622329711914   0.005781846642494
%    0.002594554126263   0.850540161132812   0.214346141815186   0.004914277493954
%    0.002406969070435   0.819778594970703   0.224235515594482   0.004851265549660
%    0.001990700364113   0.767936630249023   0.157584295272827   0.005507746338844
%    0.001853330433369   0.730196533203125   0.163685989379883   0.005587133169174
%    0.001809299886227   0.711666870117188   0.169266052246094   0.005588148236275
%    0.001854108124971   0.708142013549805   0.155819625854492   0.005426533222198
%    0.001823631674051   0.706646804809570   0.159400272369385   0.005294198989868
%    0.003650313615799   0.933787612915039   0.386913032531738   0.003863624930382
%    0.002201397716999   0.760725097656250   0.170666160583496   0.005504624843597
%    0.001652873605490   0.712392654418945   0.125616331100464   0.005676262974739
%    0.002691122889519   0.894689407348633   0.223999195098877   0.004193076491356
%    0.003131627738476   0.918962249755859   0.250651912689209   0.004305494129658
%    0.002583972215652   0.797050018310547   0.228114967346191   0.004727274775505
%    0.003313705623150   0.876754226684570   0.269729843139648   0.004171521067619
%    0.002672361135483   0.756466827392578   0.249688453674316   0.004588561654091
%    0.002775359749794   0.812133636474609   0.230609836578369   0.004862931966782
%    0.002006234377623   0.774067840576172   0.197535705566406   0.005371289253235
%    0.005882137417793   1.194994049072266   0.355792961120605   0.001550783365965
%    0.001783006787300   0.625129241943359   0.154221210479736   0.005428149700165
%    0.002347054630518   0.889358215332031   0.290890655517578   0.005299174189568
%    0.002352741211653   0.752821731567383   0.193698959350586   0.005198199748993
%    0.001897093057632   0.701995544433594   0.141393365859985   0.005235174298286
%    0.004436828196049   0.972885742187500   0.268673744201660   0.003318995833397
%    0.001829196810722   0.756867523193359   0.140361003875732   0.005725191831589
%    0.002474085092545   0.814183578491211   0.221370811462402   0.004954202771187
%    0.002611064910889   0.828627243041992   0.283689117431641   0.004908200502396
%    0.001942449957132   0.821123046875000   0.181498031616211   0.005152390599251
%    0.002821293771267   0.852723999023438   0.183055286407471   0.004918677210808
%    0.003463510274887   0.957584762573242   0.264465599060059   0.003646545112133
%    0.001795925348997   0.763667297363281   0.156026639938354   0.005376699566841
%    0.002982469499111   0.868446502685547   0.270452842712402   0.004768515825272
%    0.002087339162827   0.669192047119141   0.175448474884033   0.005191667079926
%    0.002554253935814   0.829520797729492   0.210348701477051   0.004713914096355
%    0.001996951699257   0.734728393554688   0.255351772308350   0.005567609071732
%    0.001882081180811   0.759506072998047   0.178524494171143   0.005555369257927
%    0.002300853878260   0.866504287719727   0.230094280242920   0.005164309144020
%    0.001982410103083   0.771385498046875   0.201946258544922   0.005593833327293
%    0.001971217095852   0.752721023559570   0.136845521926880   0.005300506353378
%    0.002162450253963   0.707376861572266   0.188166027069092   0.005428369045258
%    0.001903455108404   0.846447448730469   0.137851104736328   0.005597648620605
%    0.001674644201994   0.784506301879883   0.154307613372803   0.005517895221710
%    0.001641494184732   0.734223937988281   0.123774700164795   0.005777201652527
%    0.002079573124647   0.751199264526367   0.178312873840332   0.004904397726059
%    0.002242641448975   0.768721084594727   0.170476512908936   0.004988448321819
%    0.006999680995941   1.352394409179688   0.541523895263672   0.001020300462842
%    0.002951220273972   0.772367935180664   0.270432777404785   0.004954878985882
%    0.001959408074617   0.728236541748047   0.195403957366943   0.005173650383949
%    0.002355699986219   0.733626251220703   0.192397708892822   0.005258064270020
%    0.001876773536205   0.709136657714844   0.157455739974976   0.005673503279686
%    0.002988758683205   0.778314056396484   0.254636173248291   0.004953200221062
%    0.002628991901875   0.826094818115234   0.209542407989502   0.005016428232193
%    0.001925194710493   0.743201675415039   0.173249092102051   0.005701632499695
%    0.002829387187958   0.739601974487305   0.246356620788574   0.005060132741928
%    0.001662267893553   0.651847991943359   0.175295524597168   0.005833370089531
%    0.002085985988379   0.822719573974609   0.178883380889893   0.005539190173149
%    0.002198660373688   0.731158218383789   0.232807331085205   0.005354512333870
%    0.002079178243876   0.821771240234375   0.184691753387451   0.005288544893265
%    0.002984271943569   0.884551010131836   0.298673324584961   0.005266079902649
%    0.002050687968731   0.799799041748047   0.166426334381104   0.004460510611534
%    0.002922925949097   0.814297103881836   0.261519165039062   0.004977251589298
%    0.002226405441761   0.853018188476563   0.185401229858398   0.005249646902084
%    0.002290885448456   0.869483489990234   0.174292659759521   0.005189114809036
%    0.001521792709827   0.640970001220703   0.172802257537842   0.005601813197136
%    0.001441912055016   0.720282592773438   0.135819234848022   0.005590635538101
%    0.002196076512337   0.732840576171875   0.220993900299072   0.005310125946999
%    0.002671604752541   0.763802337646484   0.287058296203613   0.005006905794144
%    0.001149259582162   0.678829650878906   0.131565656661987   0.005998502969742
%    0.001752616912127   0.720417938232422   0.131578397750855   0.005393113493919
%    0.002194170802832   0.794816894531250   0.215037994384766   0.005297750830650
%    0.001425886601210   0.661991043090820   0.146538543701172   0.005736851096153
%    0.002279426157475   0.752091064453125   0.268955821990967   0.005390680432320
%    0.001655032634735   0.718994293212891   0.136022129058838   0.005859014391899
%    0.001916005313396   0.749198532104492   0.157014760971069   0.005778966546059
%    0.002461523413658   0.834672241210938   0.208973846435547   0.005359222292900
%    0.004897161424160   1.301836700439453   0.428854560852051   0.002792225480080
%    0.001285890638828   0.635278091430664   0.175784111022949   0.006003211140633
%    0.001322377920151   0.707855224609375   0.134331836700439   0.005987162590027
%    0.002274233847857   0.889659271240234   0.202166233062744   0.004751163721085
%    0.001687625646591   0.867980728149414   0.204462184906006   0.005781116485596
%    0.001192864254117   0.669517974853516   0.109100055694580   0.006420139074326
%    0.002078771293163   0.786946640014648   0.262237281799316   0.005562768578529
%    0.001476773768663   0.739194412231445   0.154030475616455   0.005972854495049
%    0.001586660295725   0.749733657836914   0.170070018768311   0.005930283069611
%    0.001987466514111   0.741544342041016   0.184753627777100   0.005511789321899
%    0.001871992051601   0.756637649536133   0.163608169555664   0.005088196396828
%    0.002369395792484   0.828762893676758   0.269111080169678   0.005183232426643
%    0.006473037600517   1.399248809814453   0.537045707702637   0.001149839237332
%    0.001655175834894   0.733142700195312   0.172572002410889   0.005649855136871
%    0.001889985650778   0.639719390869141   0.155641498565674   0.005340615510941
%    0.001689569801092   0.755743942260742   0.132851381301880   0.005720228552818
%    0.001701992601156   0.748153762817383   0.181240444183350   0.005954551696777
%    0.001559901237488   0.763537139892578   0.128192253112793   0.005491324663162
%    0.001474397033453   0.680875854492187   0.154692239761353   0.005862284302711
%    0.001825481355190   0.725878753662109   0.152633934020996   0.005875266790390
%    0.001718265265226   0.786983413696289   0.132352409362793   0.005411431789398
%    0.001813886314631   0.762698364257812   0.137820625305176   0.005757212638855
%    0.002141964137554   0.740447235107422   0.186914806365967   0.005684841871262
%    0.001700677871704   0.719710006713867   0.144664278030395   0.005778750777245
%    0.001686055213213   0.745062332153320   0.137692613601685   0.005455988645554
%    0.001760168075562   0.730490188598633   0.199453506469727   0.005379534959793
%    0.001968698501587   0.776214294433594   0.129013099670410   0.005162475705147
%    0.001723181307316   0.833144073486328   0.105890083312988   0.005201298594475
%    0.003189887702465   0.993386077880859   0.283930988311768   0.004706431925297
%    0.002854104638100   0.952137069702148   0.261803913116455   0.004613398313522
%    0.001858792752028   0.768599548339844   0.167580890655518   0.005387793183327
%    0.002078793793917   0.738973007202148   0.227171039581299   0.005408539772034
%    0.001951714158058   0.712657394409180   0.193319854736328   0.005524438619614
%    0.001753225177526   0.649506912231445   0.178899002075195   0.005890708565712
%    0.001604214161634   0.750990295410156   0.148951673507690   0.005906263589859
%    0.001759754270315   0.774657516479492   0.204484138488770   0.005697706937790
%    0.001481887698174   0.696956710815430   0.161895179748535   0.005967107415199
%    0.001936347633600   0.770792541503906   0.191983451843262   0.005432505011559
%    0.001783862859011   0.809756622314453   0.162915115356445   0.005472788810730
%    0.001983243376017   0.791453323364258   0.166432018280029   0.005259096026421
%    0.002205379307270   0.734526977539062   0.198262329101562   0.005385054349899
%    0.001567820459604   0.716034698486328   0.150166521072388   0.005732630491257
%    0.001959304362535   0.700315399169922   0.190231761932373   0.005508047342300
%    0.001596409678459   0.743970870971680   0.161814174652100   0.006083648800850
%    0.002855561077595   0.908047332763672   0.249184837341309   0.004704458117485
%    0.001859544962645   0.761653594970703   0.209201316833496   0.005765859484673
%    0.001668637990952   0.747583923339844   0.173222808837891   0.006027838587761
%    0.002164547592402   0.794438858032227   0.191956787109375   0.005344584584236
%    0.002283870279789   0.850175323486328   0.274153251647949   0.005122699141502
%    0.002166275829077   0.806646652221680   0.144293365478516   0.005301138162613
%    0.001593843102455   0.768284988403320   0.130924892425537   0.005582553744316
%    0.002386790663004   0.807865066528320   0.215204925537109   0.004988231360912
%    0.002988650202751   0.958192214965820   0.239308032989502   0.004696541428566
%    0.002219777405262   0.749241409301758   0.181094169616699   0.005423290133476
%    0.002281969934702   0.754136123657227   0.227376708984375   0.005444645285606
%    0.002773035466671   0.860654830932617   0.286503829956055   0.004855567514896
%    0.003052037656307   0.886357269287109   0.261166744232178   0.004866128265858
%    0.002342762947083   0.828863220214844   0.205348968505859   0.005209243893623
%    0.005577626824379   1.136307296752930   0.483398742675781   0.002515273690224
%    0.002426760047674   0.854063796997070   0.193480224609375   0.005285267829895
%    0.001773692071438   0.846257324218750   0.183452167510986   0.005651514530182
%    0.002713506817818   0.900202102661133   0.257664718627930   0.004568333029747
%    0.001301806271076   0.602836303710938   0.170419578552246   0.006525883674622
%    0.005933902859688   1.354922180175781   0.334792938232422   0.001740773469210
%    0.001477205008268   0.726017608642578   0.210240211486816   0.005632732510567
%    0.001595606952906   0.753806457519531   0.156767368316650   0.005726218223572
%    0.001278574168682   0.698697052001953   0.172597846984863   0.005635942220688
%    0.001929161250591   0.818841857910156   0.172611865997314   0.005583886504173
%    0.002092048227787   0.841381835937500   0.189759654998779   0.005461319088936
%    0.001174088567495   0.719706039428711   0.152408313751221   0.006229681968689
%    0.001790991425514   0.731771011352539   0.162601203918457   0.005781747698784
%    0.001859407722950   0.748509063720703   0.201705093383789   0.005590065717697
%    0.006761362552643   1.447227783203125   0.455536079406738   0.001193955019116
%    0.001724000573158   0.840342254638672   0.199953937530518   0.005445917844772
%    0.002038236707449   0.785308532714844   0.291913433074951   0.005671285986900
%    0.002381518632174   0.866777267456055   0.275528182983398   0.005252603888512
%    0.002744978368282   0.838054275512695   0.230053234100342   0.005066308379173
%    0.002467882037163   0.895702056884766   0.257476406097412   0.004965988099575
%    0.001329510062933   0.700450439453125   0.129719924926758   0.006243470907211
%    0.001807714402676   0.888740768432617   0.225183658599854   0.005859021544456
%    0.001863719224930   0.741064605712891   0.178990840911865   0.005752902030945
%    0.002432519644499   0.895694427490234   0.275720920562744   0.004764159619808
%    0.002456250190735   0.850958633422852   0.192546310424805   0.005355806946754
%    0.002554353773594   0.886196365356445   0.270433177947998   0.004836653470993
%    0.001881916522980   0.804489364624023   0.181052093505859   0.005368818044662
%    0.001828738898039   0.751479873657227   0.160576877593994   0.005495369434357
%    0.002355676889420   0.832217864990234   0.259354877471924   0.005454779267311
%    0.001431537717581   0.726471252441406   0.143964052200317   0.005865644216537
%    0.003035147786140   1.024355163574219   0.323041610717773   0.004375569820404
%    0.000970790460706   0.667745132446289   0.102611627578735   0.006363573670387
%    0.001280045509338   0.727485809326172   0.116500034332275   0.006152777075768
%    0.001471983045340   0.729331359863281   0.203019886016846   0.006340550780296
%    0.002168832719326   0.777954254150391   0.152259263992310   0.005179882049561
%    0.002556019425392   0.839663619995117   0.219248752593994   0.004774020016193
%    0.001811704337597   0.688097839355469   0.135902414321899   0.005400695800781
%    0.002175256460905   0.747894821166992   0.183915939331055   0.005133824348450
%    0.002295204699039   0.740650177001953   0.178284606933594   0.005002074241638
%    0.003614018559456   1.007287979125977   0.289425964355469   0.003676078617573
%    0.001816951036453   0.716465454101562   0.210674476623535   0.005419515371323
%    0.001528788208961   0.734369888305664   0.178184986114502   0.005736764669418
%    0.002089055329561   0.786831283569336   0.198242282867432   0.005164272785187
%    0.002049954831600   0.773386230468750   0.192629528045654   0.005337167978287
%    0.001886971145868   0.714380264282227   0.137942409515381   0.005348479747772
%    0.002043371349573   0.732397613525391   0.167271080017090   0.005450758337975
%    0.001442035585642   0.745508880615234   0.106524810791016   0.005654758214951
%    0.004823257327080   0.678952865600586   0.223686275482178   0.000570209026337
%    0.002642083764076   0.792291107177734   0.174714069366455   0.004837693274021
%    0.002409060895443   0.788797683715820   0.225107231140137   0.004286838471889
%    0.001807956248522   0.652554321289063   0.179338665008545   0.004935242831707
%    0.003763481676579   0.876801300048828   0.312540283203125   0.003424895107746
%    0.006472218632698   1.245659103393555   0.316524257659912   0.000748796686530
%    0.002333583384752   0.628033828735352   0.180465240478516   0.004569447040558
%    0.002454920411110   0.702370605468750   0.171845989227295   0.004469799697399
%    0.002012547552586   0.778100662231445   0.150486431121826   0.004521095156670
%    0.002153686136007   0.797766113281250   0.263179264068604   0.004544389843941
%    0.003015795052052   0.788273773193359   0.216045379638672   0.003897154033184
%    0.002380941212177   0.736673355102539   0.191333713531494   0.004525539875031
%    0.002624871134758   0.815162811279297   0.175766239166260   0.004137767851353
%    0.002968344986439   0.819433746337891   0.194876518249512   0.003717902004719
%    0.003840587735176   0.947921981811523   0.296318016052246   0.003581053912640
%    0.002512595355511   0.784988403320312   0.218282546997070   0.004785437583923
%    0.003045275211334   0.814701614379883   0.225617637634277   0.003946982324123
%    0.002667422294617   0.741512527465820   0.209121227264404   0.004423930048943
%    0.002663233876228   0.790732498168945   0.180647659301758   0.003820178806782
%    0.002341687679291   0.747829132080078   0.202084407806396   0.005069301724434
%    0.002187212109566   0.715945739746094   0.163935317993164   0.004835323095322
%    0.002007242441177   0.772592163085937   0.144200086593628   0.005168022513390
%    0.002089520692825   0.759936370849609   0.200063800811768   0.005108578205109
%    0.002540862262249   0.762660064697266   0.200694122314453   0.004460925459862
%    0.002332108020782   0.735709915161133   0.163983345031738   0.004672377109528
%    0.002529984414577   0.767191848754883   0.215460205078125   0.004401385188103
%    0.001750732064247   0.752673950195312   0.124574728012085   0.005513265728951
%    0.001885639876127   0.744115295410156   0.174471225738525   0.005879573822021
%    0.002246018797159   0.783644027709961   0.239298076629639   0.004751084148884
%    0.001763053536415   0.725673904418945   0.154617109298706   0.005676690936089
%    0.001497658938169   0.711598129272461   0.177964763641357   0.005869375467300
%    0.002150685638189   0.687253570556641   0.149532880783081   0.004932304322720
%    0.001923832595348   0.775416412353516   0.220050888061523   0.005360119938850
%    0.005785967111588   1.221675262451172   0.498089370727539   0.002166681587696
%    0.002681432664394   0.830811004638672   0.202104949951172   0.004963640272617
%    0.002050870060921   0.717255096435547   0.157643547058105   0.005279009938240
%    0.001949970573187   0.678344345092773   0.185980777740479   0.005423024892807
%    0.002135311514139   0.775394744873047   0.214016895294189   0.005376824140549
%    0.001429180949926   0.537669486999512   0.156442937850952   0.005915979743004
%    0.003681410253048   0.992977752685547   0.310564994812012   0.003786623179913
%    0.002838189005852   0.847455444335937   0.295263977050781   0.004821552634239
%    0.002206949740648   0.715609054565430   0.172367210388184   0.004784561991692
%    0.002095760256052   0.703798141479492   0.168049926757812   0.005433177947998
%    0.003040748238564   0.848625793457031   0.184440822601318   0.004330361485481
%    0.003350003361702   0.824446945190430   0.376857528686523   0.004478645622730
%    0.001724840849638   0.796805267333984   0.145400848388672   0.005643475055695
%    0.002643606364727   0.806142120361328   0.241654567718506   0.004817925989628
%    0.002346525490284   0.747067260742188   0.190371246337891   0.004848866760731
%    0.002180344313383   0.767978668212891   0.174122066497803   0.004978697597980
%    0.001908503025770   0.734573593139648   0.142734746932983   0.005324305295944
%    0.001706966310740   0.684514312744141   0.148259601593018   0.006191613674164
%    0.001135899871588   0.602948493957519   0.133370780944824   0.006293917894363
%    0.003443238437176   0.975120773315430   0.325404472351074   0.004184302985668
%    0.001743954122066   0.780666732788086   0.175102882385254   0.005231199860573
%    0.003021245598793   0.785602569580078   0.176368350982666   0.004513958096504
%    0.002342651486397   0.824602279663086   0.204196128845215   0.004922202229500
%    0.001608573198318   0.744301071166992   0.172314624786377   0.006025561094284
%    0.002928221523762   0.865926437377930   0.225158710479736   0.004862830638885
%    0.001485214084387   0.746906585693359   0.126213703155518   0.005394389033318
%    0.002309937775135   0.749375228881836   0.189328937530518   0.004895468056202
%    0.001655221730471   0.730835723876953   0.129226312637329   0.005489861369133
%    0.003388043045998   0.957310028076172   0.281908779144287   0.004054368734360
%    0.001823593676090   0.773955993652344   0.165788803100586   0.005079376101494
%    0.001512543559074   0.704496841430664   0.120081920623779   0.005794232487679
%    0.002065167874098   0.780681991577148   0.175946578979492   0.005207548141479
%    0.002169333100319   0.748608398437500   0.170282306671143   0.005170695781708
%    0.002592226266861   0.908378295898438   0.227902297973633   0.004397542476654
%    0.002280317246914   0.816022186279297   0.182626991271973   0.004971090555191
%    0.001729563623667   0.669891357421875   0.190460166931152   0.005375084877014
%    0.002670497298241   0.831127853393555   0.206176185607910   0.004545546472073
%    0.002058706283569   0.770430221557617   0.188854961395264   0.005065277218819
%    0.001916195452213   0.743923797607422   0.160065898895264   0.005285060405731
%    0.002281710803509   0.787353744506836   0.189280395507813   0.005188954472542
%    0.001774245053530   0.719280929565430   0.142155418395996   0.005410445332527
%    0.001661729961634   0.765912780761719   0.110991392135620   0.005681448578835
%    0.001155662536621   0.597701911926270   0.137863388061523   0.005949981808662
%    0.003435553610325   0.961252517700195   0.295714759826660   0.004150130152702
%    0.001587293148041   0.695688018798828   0.152144203186035   0.005442665815353
%    0.002313528656960   0.843390502929688   0.241500701904297   0.005474021434784
%    0.002791703641415   0.910098037719727   0.200268325805664   0.004541390240192
%    0.002255451977253   0.808373336791992   0.207123985290527   0.004689464867115
%    0.002396699339151   0.853402328491211   0.222364845275879   0.004554262757301
%    0.001868878751993   0.823105545043945   0.202897529602051   0.005382898449898
%    0.001468473225832   0.748426589965820   0.198316249847412   0.005961800813675
%    0.001861634403467   0.805653533935547   0.179422187805176   0.005489081144333
%    0.001858018487692   0.765192565917969   0.195187263488770   0.005945045948029
%    0.002071690261364   0.762838897705078   0.174249725341797   0.005343416333199
%    0.002682299911976   0.859430084228516   0.266433258056641   0.004534697830677
%    0.001209836229682   0.692791976928711   0.121507539749146   0.005594403743744
%    0.001633160263300   0.711151046752930   0.141865310668945   0.005915149450302
%    0.002783655524254   0.878065872192383   0.277172451019287   0.004544402956963
%    0.002577842175961   0.816831207275391   0.236943740844727   0.004850858151913
%    0.004151795208454   1.188518829345703   0.322266960144043   0.003027047514915
%    0.002679159045219   0.947984161376953   0.253051300048828   0.005047816634178
%    0.001975617110729   0.812634582519531   0.215456314086914   0.005907616019249
%    0.001806811243296   0.808981933593750   0.161039009094238   0.005806286931038
%    0.003150395452976   1.073191909790039   0.237876243591309   0.003951642513275
%    0.002464054375887   0.840247344970703   0.231811714172363   0.005083907246590
%    0.001373518705368   0.759236145019531   0.205828952789307   0.006199684143066
%    0.006088508367538   1.354408569335938   0.374663658142090   0.001624601930380
%    0.001703419387341   0.775193252563477   0.202797889709473   0.005720340609550
%    0.002085740268230   0.819786453247070   0.201136360168457   0.005507355928421
%    0.002124632894993   0.724251327514648   0.277698955535889   0.005442286729813
%    0.002005890607834   0.744687118530273   0.184660282135010   0.005849013924599
%    0.001711651831865   0.731653060913086   0.143770837783813   0.005531733632088
%    0.002525928020477   0.862827606201172   0.306545047760010   0.005122088789940
%    0.001689728647470   0.781366577148437   0.197519760131836   0.005742967128754
%    0.001967422962189   0.753151168823242   0.231859741210937   0.005331518650055
%    0.001680722087622   0.718017730712891   0.188857860565186   0.005176895260811
%    0.002064693421125   0.767464904785156   0.178155288696289   0.004980225861073
%    0.002719770967960   0.940613937377930   0.282324104309082   0.004900139868259
%    0.001841060668230   0.752931365966797   0.199442539215088   0.005524509549141
%    0.002666705846786   0.860817184448242   0.271520442962646   0.005112269520760
%    0.001451287120581   0.712890167236328   0.111643962860107   0.005439339876175
%    0.001728678345680   0.764643020629883   0.239502010345459   0.005629344582558
%    0.005684335231781   1.441481018066406   0.427109909057617   0.002292271256447
%    0.001136224344373   0.691211090087891   0.118255605697632   0.006425678133965
%    0.004872546195984   1.296812133789063   0.284832363128662   0.002323878258467
%    0.002057639807463   0.867603759765625   0.185529041290283   0.005571421384811
%    0.001084020659328   0.643941879272461   0.128196544647217   0.006526120305061
%    0.001301031708717   0.697323989868164   0.112915735244751   0.006221711039543
%    0.001613764613867   0.734308853149414   0.139777717590332   0.005803414583206
%    0.001757181286812   0.743041687011719   0.161681900024414   0.005865421295166
%    0.000824753418565   0.713899230957031   0.100924978256226   0.005976176857948
%    0.002321444898844   0.919589157104492   0.365824050903320   0.004877621531487
%    0.002236461937428   0.969893112182617   0.246821479797363   0.005231838226318
%    0.001484689116478   0.738187789916992   0.166085834503174   0.005758491754532
%    0.001732496321201   0.713068695068359   0.164818058013916   0.005971375107765
%    0.001640295684338   0.704776763916016   0.199794635772705   0.006068143844604
%    0.001716905236244   0.794581451416016   0.177369480133057   0.006113269925117
%    0.001411285847425   0.701635894775391   0.166625213623047   0.006369696855545
%    0.002603221535683   0.856117630004883   0.275933189392090   0.005273547768593
%    0.004107212722301   1.008299789428711   0.394798965454102   0.003674463629723
%    0.002394996881485   0.897696456909180   0.228619194030762   0.005322678089142
%    0.003505154252052   0.959912033081055   0.293100299835205   0.004593545794487
%    0.001774685978889   0.747255935668945   0.162584018707275   0.005513183474541
%    0.001704669892788   0.734908523559570   0.186012592315674   0.005518343448639
%    0.005201699137688   1.435593109130859   0.529844894409180   0.002621654570103
%    0.003956516683102   1.356957244873047   0.344356002807617   0.003923469185829
%    0.002911441624165   0.952355651855469   0.314488887786865   0.005140407085419
%    0.001702211648226   0.732353134155273   0.222830867767334   0.006075572371483
%    0.001243066340685   0.702229232788086   0.094103717803955   0.006525958180428
%    0.001048780903220   0.749705352783203   0.155408477783203   0.006326518058777
%    0.001617356389761   0.830166854858398   0.166012477874756   0.005597482323647
%    0.001650708913803   0.703884353637695   0.148795413970947   0.005791711211205
%    0.001356547772884   0.733748855590820   0.141136522293091   0.005853787064552
%    0.001816833466291   0.785151672363281   0.150845432281494   0.005398336648941
%    0.002752310037613   0.737258377075195   0.213875808715820   0.004130223393440
%    0.002052751332521   0.666280441284180   0.139309654235840   0.004700674116611
%    0.002613074779510   0.831892852783203   0.247707061767578   0.004474525749683
%    0.002171046584845   0.698366012573242   0.232350921630859   0.005392389297485
%    0.005053592324257   1.108731613159180   0.311181449890137   0.002320296615362
%    0.001587955504656   0.749605331420898   0.132744789123535   0.005806143283844
%    0.002385502308607   0.796761627197266   0.212878246307373   0.004674561619759
%    0.001577220857143   0.700052642822266   0.135017786026001   0.005554210543633
%    0.002254103422165   0.825198211669922   0.223241233825684   0.005402856469154
%    0.002384315580130   0.847536315917969   0.228363285064697   0.005013938546181
%    0.005482913255692   1.340830993652344   0.361918983459473   0.001860444396734
%    0.001760503798723   0.667162399291992   0.158011980056763   0.005223993062973
%    0.001623573005199   0.662318725585938   0.135490617752075   0.005489204525948
%    0.002238596230745   0.805527191162109   0.223960018157959   0.005134783983231
%    0.007664164304733   1.397297668457031   0.410473175048828   0.000192923340946
%    0.001576955616474   0.730380020141602   0.180856113433838   0.006158912181854
%    0.001782636642456   0.735988159179688   0.160094070434570   0.005467735528946
%    0.002236324697733   0.822190322875977   0.220435543060303   0.004950608611107
%    0.002131938338280   0.790766448974609   0.186705112457275   0.004829219579697
%    0.004157951772213   1.036477050781250   0.416349411010742   0.003293715715408
%    0.002165467590094   0.735257339477539   0.222224769592285   0.005156137347221
%    0.002390152364969   0.763685150146484   0.280283203125000   0.004745567739010
%    0.001511449664831   0.752103576660156   0.138678522109985   0.005183366537094
%    0.004975223243237   1.313986206054687   0.516453170776367   0.002428002059460
%    0.001979831457138   0.728271102905273   0.138630733489990   0.005479627251625
%    0.001614513397217   0.709712753295898   0.145726280212402   0.005542614459991
%    0.001768473386765   0.720813293457031   0.173879737854004   0.005823337435722
%    0.001785129904747   0.773252716064453   0.164297237396240   0.005266921520233
%    0.001217740625143   0.620577621459961   0.223407993316650   0.006000456213951
%    0.001102772802114   0.693284835815430   0.162398128509521   0.006367119550705
%    0.001595518141985   0.775068893432617   0.142764987945557   0.005944435000420
%    0.001294711232185   0.741084289550781   0.128252687454224   0.006135792732239
%    0.001382003873587   0.665839004516602   0.125552120208740   0.006138646006584
%    0.004080435633659   1.147267532348633   0.462971687316895   0.003795357048512
%    0.001151531338692   0.675037994384766   0.099849014282227   0.006176242232323
%    0.001724790334702   0.779597167968750   0.181343326568604   0.005944886803627
%    0.005223070383072   1.208223190307617   0.259674186706543   0.002052983343601
%    0.001063861027360   0.740437850952148   0.115399150848389   0.006436208486557
%    0.002135568708181   0.907623748779297   0.246829185485840   0.005545855760574
%    0.001355901211500   0.671952514648438   0.148086748123169   0.006635276675224
%    0.001637690067291   0.769920501708984   0.208036441802979   0.005915223360062
%    0.001165803447366   0.797534484863281   0.194498138427734   0.006000462770462
%    0.001325949132442   0.874669647216797   0.219196987152100   0.006480336785316
%    0.001883194744587   0.834361801147461   0.244827175140381   0.005860878825188];


% ----------------------(eliminate multiple peaks)--------------------

% criterion1: K must be POSITIVE (K>0)! Eliminate all pksNew with
% corresponding K values that are <=0 OR K == NaN (note: K == NaN means
% that nlinfit was not able to solve for parameters K, D, w, offset)
% 
% The brightest regions of a particle are towards it center...
% 
% K < 0 means that the sphere spread function(ssf) for the particle is 
% inverted so that the center of a particle is the darkest regions of the
% particle.  This is NOT possible in Inorm!
% 
% K = 0 means that the ssf of a particle is flat (i.e., the spherical
% profile is a horizontal line). This also is impossible! 
 
pksToEliminate = find(K(:)<=0 | isnan(K)); %there aren't any



% criterion2: Particle size is known (i.e., 5um +/-0.44um silica,
% manufactured by Bang Laboratories). All particles identified by 'pksNew'
% should have sizes ~Dactual +/- 10% of Dactual 
% 
% IMPORTANT NOTE: Dactual or D are really just parameters of ssf() NOT
% particle size.  You found empirically that particle size is usually some
% combination of D and w
% 
% But some particles are slightly less than 10% larger or smaller than
% Dactual. SO YOU MUST PLOT THE PARTICLE SIZES AS DONE BELOW, MANUALLY 
% LOOK AT THESE PARTICLES, AND ACCORDINGLY ADJUST lBound and uBound TO
% EITHER KEEP OR DELETE THEM. THIS MANUAL ADJUSTMENT MUST BE DONE WITH EACH
% PASS OF YOUR PARTICLE TRACKING ROUTINE!

lBound = 0.7*Dactual;
uBound = 1.3*Dactual;
lowerBound = lBound*ones([size(D,1),1],'single');
upperBound = uBound*ones([size(D,1),1],'single');

figure(1)
hold off
plot(D,'b-*')
hold on
plot(lowerBound,'r--')
plot(upperBound,'r--')
title(sprintf('D = 5.06 +/- 0.44 um according to Bang Labs',Dactual),'fontsize',20)
xlabel(sprintf('actually plotted:\n (%.03f*Dactual), r-- upper \n(%.03f*Dactual), r-- lower \n Dactual = %.02f pixels (size of ideal maskExp)',(uBound/Dactual),(lBound/Dactual),Dactual),'fontsize',20)

pksToEliminate = find(D>uBound | D<lBound); %49 particles

for np = 1:length(pksToEliminate)
    [v,center] = getCutoutSectVertices(pksNewer(pksToEliminate(np),:),D(pksToEliminate(np)),size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    figure(2)
    hold off
    h(1) = plot(rNew,avgProfile(startIdx:stopIdx),'b*');
    hold on
    h(2) = plot(rNew,ssf([K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np))],rNew),'k-');
    legend(h,'InormRed','mask')
    title(sprintf('pksNewer(%d of %d,:) = [%.04f, %.04f, %.04f] \n K = %.04f, D = %.04f, w = %.04f, offset = %.04f \n bad D because D < %.04f OR D > %.04f',pksToEliminate(np),size(pksNewer,1),pksNewer(pksToEliminate(np),1),pksNewer(pksToEliminate(np),2),pksNewer(pksToEliminate(np),3),K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),lowerBound(1),upperBound(1)),'fontsize',20)
    xlabel(sprintf('Rsquared = %.04f',Rsquared(pksToEliminate(np))),'fontsize',20)

    [voronoiVol, Rmask] = peakPlacementB(center,size(cutoutSect),D(pksToEliminate(np)),w(pksToEliminate(np)));
    mask = getCalcImg(Rmask,K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),voronoiVol);
    
    sphericalRegion = find(voronoiVol == 1);
    diff = cutoutSect;
    diff(sphericalRegion) = abs(cutoutSect(sphericalRegion)-mask(sphericalRegion));       
    
    b = 1:numel(mask); 
    b(sphericalRegion) = []; %indices to identify region of the image where there isn't a particle
    mask(b) = cutoutSect(b);
    
    figure(3)
    hold off
    simage([cutoutSect(:,:,round(center(3))), max(cutoutSect(:))*ones([size(cutoutSect,1) 1],'single'), mask(:,:,round(center(3))), max(cutoutSect(:))*ones([size(cutoutSect,1) 1],'single'), diff(:,:,round(center(3)))]);
    colorbar
    hold on
    plot(center(2),center(1),'k*','markersize',20)
    title(sprintf('InormRed \t mask \t abs(mask - InormRed)'),'fontsize',20)
    xlabel(sprintf('pksNewer(%d of %d,:) = [%.04f, %.04f, %.04f] \n K = %.04f, D = %.04f, w = %.04f, offset = %.04f \n bad D because D < %.04f OR D > %.04f',pksToEliminate(np),size(pksNewer,1),pksNewer(pksToEliminate(np),1),pksNewer(pksToEliminate(np),2),pksNewer(pksToEliminate(np),3),K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),lowerBound(1),upperBound(1)),'fontsize',20)
    pause
end


% submit one job at a time to be computed by all available slave nodes
for j = 1:length(job)
    submit(job(j));
    disp(sprintf('job%02d, submitted',j))
end 

% 1st job to start, started at: Tue Aug 04 23:54:10 EDT 2015
% last job to finish, finished at: Wed Aug 05 10:48:57 EDT 2015
% total time: 10hrs54min47seconds for 421 particles

pksBest = zeros(size(pksNewer),'single');
Rsquared = zeros([11,size(pksNewer,1)],'single'); 
for j = 1:length(job)
    state = job(j).State;
    
    startTask = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        try 
            results = fetchOutputs(job(j));
            disp(sprintf('\n------------job%02d, COMPLETED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksBest(startTask+np-1,:) = results{np,1}; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = results{np,2};
                disp(sprintf('pksNewer(%d,:) becomes pksBest(%d,:)',startTask+np-1,startTask+np-1));
            end
        catch err 
            disp(sprintf('\n------------job%02d, FAILED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksBest(startTask+np-1,:) = NaN; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = NaN;
                disp(sprintf('pksNewer(%d,:) = NaN',startTask+np-1));
            end
        end%end of try
    end
end

temp = isnan(pksBest(:,1));
npRedo = find(temp(:) == 1); %There aren't any!

% pksBest =(1.0e+03)*[
% 
%    0.0539843   0.1150039   0.2449922
%    0.1574941   0.1100275   0.2299922
%    0.3559882   0.1049804   0.2609961
%    0.3690823   0.2040039   0.2469843
%    0.4379216   0.1506980   0.1940000
%    0.7540627   0.1419961   0.2090000
%    0.8700000   0.0464824   0.2390196
%    0.7283490   0.2322980   0.2420000
%    1.0010000   0.1929529   0.1390000
%    0.9970000   0.0993765   0.1889961
%    0.0150000   0.3439568   0.1350000
%    0.4080000   0.4479922   0.2538667
%    0.4545451   0.3150039   0.2490000
%    0.5959177   0.5090000   0.2479961
%    0.6805020   0.3570157   0.2110000
%    0.8050078   0.2948235   0.2130118
%    0.9450588   0.5380510   0.2302980
%    0.9380078   0.3150039   0.2160118
%    0.0190000   0.6970118   0.2182392
%    0.1868588   0.6380706   0.2270000
%    0.3609804   0.5308000   0.2379961
%    0.2390000   0.5929529   0.3259961
%    0.6200000   0.6640000   0.2825843
%    0.8039882   0.5789961   0.2039569
%    0.7590000   0.7480432   0.2289725
%    0.0280000   0.8109177   0.2100118
%    0.3465882   0.9950000   0.1370000
%    0.2490039   0.8299922   0.1710784
%    0.2210000   0.9409961   0.2480275
%    0.4869922   0.8930000   0.2199882
%    0.6230784   0.6900039   0.1830275
%    0.7390000   0.8440314   0.2569177
%    0.8489843   0.9760039   0.2439725
%    0.4589843   0.9570157   0.2859765
%    0.5510000   0.9703530   0.2250078
%    0.4557333   0.0759686   0.2800000
%    0.6349961   0.0019961   0.2409255
%    0.9579882   0.1540000   0.2530118
%    0.0509843   0.4870039   0.2190275
%    0.7632549   0.7620392   0.3080118
%    0.8429922   0.9850000   0.1400078
%    0.0389961   0.0599294   0.3709451
%    0.1048549   0.2470039   0.3579961
%    0.1199922   0.0060000   0.3940510
%    0.3119882   0.0349922   0.3790392
%    0.4160000   0.1660039   0.3120039
%    0.4682510   0.0019961   0.3280078
%    0.6861921   0.0910000   0.3690078
%    0.4860000   0.0519373   0.4620039
%    0.5750039   0.0129725   0.4560157
%    0.7490000   0.0239961   0.4008549
%    0.6469961   0.0990196   0.4569686
%    0.7729961   0.0589961   0.2629529
%    0.7969922   0.0920078   0.3568784
%    0.7206706   0.1410275   0.2900157
%    0.8529804   0.0820000   0.4560275
%    0.9399255   0.2430078   0.4191020
%    0.9279961   0.1370000   0.4899451
%    0.1931569   0.2609922   0.3139804
%    0.1479882   0.5017686   0.3730039
%    0.1419608   0.3759177   0.4840118
%    0.0889961   0.4460157   0.4500196
%    0.4660196   0.5060039   0.3861647
%    0.3889882   0.3540000   0.3107255
%    0.3095882   0.4240078   0.4620000
%    0.4349961   0.2560157   0.3270275
%    0.7200000   0.3769647   0.3537843
%    0.6090000   0.4566471   0.4440000
%    0.8500078   0.3280078   0.3099568
%    0.8119922   0.4150078   0.3268863
%    1.0230000   0.3569647   0.3539686
%    0.9630157   0.4279529   0.3850000
%    0.0599137   0.7150039   0.3559882
%    0.1610431   0.6460196   0.3160000
%    0.2480667   0.5672353   0.4819961
%    0.1278745   0.6950000   0.4766196
%    0.3750039   0.6020196   0.4820980
%    0.4710039   0.6064980   0.4250118
%    0.6356706   0.6914510   0.3760000
%    0.7110118   0.6680000   0.3140196
%    0.6849647   0.4380000   0.5050039
%    0.6260549   0.6969961   0.4699961
%    0.8079882   0.4730078   0.4085020
%    0.7290196   0.5480236   0.4500000
%    0.9089961   0.6830157   0.4849882
%    1.0040039   0.6140000   0.2691843
%    0.9189961   0.5380118   0.5259804
%    0.9980000   0.6580784   0.4506941
%    0.1109843   0.8329804   0.3560039
%    0.2669765   0.9460157   0.3669804
%    0.1210000   0.7790039   0.5339647
%    0.1179922   0.8680039   0.4648510
%    0.1899882   0.9100000   0.5309529
%    0.3780000   0.8377686   0.3580196
%    0.3689961   0.7034039   0.4729922
%    0.4490235   0.8470000   0.3149961
%    0.4629568   0.9390078   0.3680039
%    0.6908823   0.9950000   0.3709882
%    0.7160000   0.9300275   0.4590902
%    0.9579961   0.8140078   0.3569882
%    0.8169804   0.8889961   0.3019961
%    0.8319686   0.7870000   0.5119843
%    0.9169568   0.8030078   0.4350039
%    0.8709961   0.8650157   0.3690196
%    1.0020000   0.8451216   0.4349804
%    0.9416981   0.8600000   0.5119373
%    0.1759765   0.9930039   0.3759804
%    0.1349922   0.9650039   0.4570039
%    0.2340078   0.9660118   0.4539725
%    0.6120745   0.9470157   0.3410000
%    0.8729490   0.9550039   0.3370000
%    0.9060000   0.9520039   0.4479059
%    0.1630000   0.1340118   0.6349294
%    0.0639922   0.1070157   0.6521451
%    0.3693726   0.1070000   0.5470000
%    0.2560000   0.1581294   0.5978667
%    0.4489804   0.2549961   0.6289961
%    0.2532353   0.0299961   0.6660157
%    0.4599490   0.1010078   0.6919961
%    0.3370157   0.1290000   0.6609882
%    0.4247372   0.2030000   0.7670627
%    0.5490627   0.1020000   0.5110000
%    0.7390157   0.0852157   0.4829922
%    0.5448079   0.1370471   0.6000039
%    0.5614353   0.1165020   0.6946627
%    0.6950196   0.0919961   0.6839882
%    0.8290196   0.0149961   0.5629177
%    0.8109725   0.1089961   0.5509764
%    0.8649922   0.2120000   0.5710353
%    0.8879059   0.0739961   0.7360000
%    0.8830000   0.2710039   0.6894078
%    1.0020078   0.0910157   0.5991843
%    0.9457961   0.2159804   0.6461294
%    0.9507647   0.2524627   0.7810078
%    0.1490000   0.3120000   0.5719647
%    0.2360000   0.3220000   0.4879843
%    0.0390000   0.3159765   0.7399961
%    0.2099059   0.3180118   0.7380000
%    0.4287843   0.3610000   0.6029961
%    0.3480000   0.4450039   0.5713490
%    0.2849922   0.3919961   0.7355176
%    0.4999882   0.2650078   0.7439922
%    0.3599725   0.4020000   0.6628981
%    0.3510000   0.4529608   0.7440157
%    0.4500432   0.4254941   0.6749882
%    0.6098941   0.3530000   0.6150627
%    0.5179843   0.4599804   0.4689765
%    0.5735412   0.4570588   0.5480039
%    0.6430118   0.3180118   0.7589098
%    0.7309725   0.2429922   0.6889922
%    0.6940000   0.4959529   0.7010392
%    0.7900196   0.4995608   0.5043569
%    0.8599137   0.4710157   0.5769961
%    0.9239843   0.2930000   0.5909922
%    0.9940000   0.3670000   0.5659490
%    0.1782314   0.4750039   0.4560039
%    0.1909569   0.4080078   0.7019961
%    0.1000000   0.7229922   0.7070236
%    0.2810000   0.5670902   0.6195686
%    0.3666510   0.5550039   0.6210000
%    0.3189961   0.5321843   0.7043726
%    0.3630039   0.6170039   0.6880823
%    0.4439843   0.6880078   0.6889922
%    0.4789882   0.5485647   0.5949882
%    0.6201530   0.6430078   0.6100275
%    0.5630667   0.6679255   0.5439764
%    0.5273137   0.6630000   0.6430588
%    0.8440275   0.6670000   0.5594588
%    0.9179137   0.6740000   0.7310196
%    0.9630157   0.6370000   0.5597098
%    0.9800353   0.5330000   0.6970823
%    0.9880000   0.7069922   0.6634785
%    0.1320000   0.8170000   0.7160980
%    0.0379961   0.8730157   0.7160000
%    0.1981059   0.9629922   0.6130000
%    0.3903372   0.6755137   0.6010000
%    0.2730235   0.9600078   0.5420078
%    0.4300000   0.8010392   0.7650118
%    0.5436196   0.9560000   0.6250000
%    0.7249686   0.9300000   0.5590314
%    0.5990000   0.7238823   0.6760549
%    0.5469961   0.9060196   0.7024785
%    0.8009608   0.8650118   0.5550706
%    0.8109882   0.7660000   0.6029804
%    0.8090000   0.9370510   0.7099961
%    0.9469843   0.9359961   0.5769961
%    0.1090000   0.9463961   0.5860078
%    0.0860039   0.0760431   0.8320000
%    0.2740000   0.0780000   0.8260078
%    0.4809882   0.0319882   0.7500196
%    0.6250196   0.0890000   0.7600000
%    0.7439961   0.1480157   0.7280000
%          NaN         NaN         NaN
%          NaN         NaN         NaN
%          NaN         NaN         NaN
%          NaN         NaN         NaN
%          NaN         NaN         NaN
%          NaN         NaN         NaN
%    0.7680039   0.4990000   0.7659529
%    0.9360314   0.3480000   0.7550078
%    0.1460510   0.4420039   0.7760314
%    0.1980000   0.7509412   0.6940039
%    0.2670000   0.6160196   0.7058588
%    0.4070039   0.4737882   0.8177568
%    0.4321608   0.7069882   0.8320000
%    0.5050196   0.7400000   0.7371608
%    0.5670000   0.6689922   0.7710157
%    0.2120078   0.7647961   0.7860000
%    0.6310196   0.7459725   0.7660000
%    0.6060039   0.8220588   0.7039961
%    0.8880039   0.8059804   0.7675372
%    0.9760000   0.7300078   0.8230432
%    0.1114824   0.9144941   0.7620078
%    0.1620039   0.9400236   0.8300039
%    0.4830039   0.9599686   0.8320078
%    0.7299961   0.9190039   0.8100000
%    0.8000118   0.9870039   0.8260000
%    0.3829882   0.5490784   0.4070000
%    0.2550275   0.5264235   0.3939961
%    0.3582392   0.3752118   0.3900078
%    0.4299216   0.3066941   0.4120000
%    0.5170510   0.2930039   0.3630000
%    0.6919843   0.1920314   0.3850196
%    0.9220078   0.0269882   0.4320118
%    0.3918510   0.0150000   0.4369961
%    0.0389961   0.2454471   0.4230000
%    0.1063020   0.3150000   0.4228039
%    0.0890078   0.6330039   0.4089608
%    0.4389882   0.4274196   0.4190039
%    0.6316078   0.5670078   0.4530039
%    0.6210157   0.8890000   0.4320039
%    0.5519961   0.9500000   0.4139608
%    0.1759922   0.0559882   0.1400157
%    0.0631412   0.0209961   0.1460078
%    0.1720314   0.1780000   0.1390000
%    0.2986000   0.2799804   0.1440000
%    0.3529686   0.3850000   0.1430118
%    0.5148784   0.3678039   0.1430000
%    0.5940000   0.2285020   0.1500000
%    0.6107608   0.0549843   0.1430039
%    0.3347490   0.6800118   0.1480039
%    0.3820353   0.5770000   0.1519765
%    0.5279922   0.6750196   0.1460275
%    0.6019843   0.7810000   0.1460549
%    1.0040118   0.7860471   0.1420118
%    0.9810745   0.9530000   0.1460000
%    0.2050000   0.9050157   0.1401333
%    0.4300000   0.9010236   0.1408235
%    0.5700000   0.8909804   0.1489882
%    0.1271137   0.3400039   0.3299765
%    0.2420235   0.6940078   0.3459922
%    0.3779961   0.9680039   0.3374745
%    0.1859529   0.9019922   0.3959882
%    0.2670000   0.8180275   0.3769647
%    0.3349922   0.7499686   0.3769843
%    0.4170941   0.6939764   0.3660000
%    0.5400236   0.4330078   0.3785294
%    0.5350000   0.0400039   0.3909098
%    0.9229961   0.3389373   0.4040000
%    0.9589922   0.9100000   0.3736588
%    1.0110000   0.0910078   0.3136039
%    0.7749882   0.5000078   0.3189922
%    0.8809137   0.0339961   0.3549843
%    0.5600000   0.0699961   0.2999922
%    0.4049568   0.5220078   0.3170039
%    0.4519961   0.4220000   0.3320157
%    0.5308196   0.5373373   0.3200000
%    0.0110000   0.8039764   0.3766039
%    0.4140000   0.0669922   0.3673569
%    0.7790000   0.9450039   0.3659765
%    0.5320118   1.0140078   0.3214039
%    0.6990000   0.8977177   0.3409961
%    0.9629647   0.6547294   0.3509882
%    0.6130353   0.7890000   0.4450000
%    0.4200432   0.7729961   0.4349804
%    0.4370000   0.8670000   0.4369686
%    0.5242078   0.1350000   0.4200000
%    0.2860039   0.6270000   0.4065216
%    0.1109529   0.5400039   0.4450039
%    0.1386314   0.2250000   0.4410432
%    1.0080039   0.9572157   0.4560314
%    1.0020236   0.5639882   0.4790353
%    0.7769764   0.9950000   0.4959725
%    0.4100706   1.0080039   0.4849294
%    0.5420000   0.8459216   0.4790039
%    0.6499647   0.9740118   0.5140000
%    1.0140000   0.0970471   0.4990039
%    0.9630118   0.4100000   0.4709765
%    0.0109843   0.6650549   0.5179961
%    0.2056431   0.1470078   0.5109961
%    0.4549686   0.1420078   0.4840157
%    0.6099922   0.3779961   0.4989686
%    0.8589647   0.3860078   0.4800000
%    0.8790000   0.2949177   0.5070745
%    0.7200000   0.8429764   0.4970275
%    0.2700000   0.7140667   0.4886588
%    0.4019922   0.2400196   0.5529608
%    0.4960039   0.2485608   0.5450000
%    0.0450000   0.8890196   0.5437333
%    0.2390000   0.8180000   0.5155726
%    0.5080000   0.7790000   0.5444588
%    0.5560078   0.5709843   0.5299725
%    0.7059922   0.5389804   0.5380314
%    0.6760000   0.6250039   0.5339961
%    0.3171765   0.5360078   0.5380000
%    0.9970000   0.4820000   0.5340000
%    0.9250118   0.7630078   0.5404078
%    0.6254863   0.1892784   0.5520039
%    0.8999686   1.0160000   0.5399608
%    0.4429961   0.4731961   0.5368196
%    0.8250000   0.5740000   0.5583176
%    0.0669922   0.2150392   0.5909843
%    0.4220078   0.1619922   0.6010078
%    0.6039882   0.0220000   0.5896000
%    0.9079764   0.0630196   0.6039961
%    0.9870549   0.5495137   0.5940039
%    0.6860000   0.8233726   0.5970000
%    0.3229804   0.7390118   0.5980078
%    0.5639804   0.2450000   0.6141961
%    0.8538000   0.9205098   0.6170000
%    0.1644549   0.2360000   0.6370314
%    0.1639961   0.4810275   0.6409843
%    0.6269961   1.0110000   0.6501137
%    0.8940392   0.8189922   0.6170118
%    0.9000471   1.0030000   0.6334196
%    1.0080000   0.8730471   0.6299843
%    0.0049961   0.2790039   0.6547372
%    0.5050118   0.3200039   0.6539922
%    0.9349804   0.0069922   0.6700392
%    0.9840000   0.2949608   0.6829804
%    0.8870000   0.5070196   0.6609882
%    0.5139725   0.8159529   0.6679961
%    0.4649686   0.5910510   0.6700000
%    0.5529686   0.4386784   0.6489922
%    0.8014353   0.6749373   0.6480000
%    0.8330000   0.3760196   0.8130353
%    0.3075059   0.0450000   0.7290157
%    0.3750275   0.1140039   0.7430118
%    0.3340000   0.2060353   0.7427608
%    0.3009490   0.2939608   0.7210000
%    0.2109804   0.0640000   0.7487372
%    0.1139922   0.0430157   0.7450000
%    0.0639961   0.1520000   0.7449686
%    0.0480157   0.6390039   0.7470510
%    0.0130000   0.0719961   0.7595412
%    0.8880039   0.5780392   0.7390039
%    0.7310078   0.9980000   0.7440510
%    0.6819961   0.0099922   0.7730039
%    0.8930078   0.9720275   0.7270000
%    1.0160118   0.3390000   0.8030275
%    0.4999922   0.3002588   0.8220078
%    0.3080314   0.3354471   0.8050000
%    0.2140000   0.3750196   0.8191725
%    0.5649804   1.0150118   0.8260078
%    0.3870118   0.0179882   0.7940039
%    0.8049961   0.0299882   0.7800196
%    0.8149961   0.1290353   0.7900157
%    0.8567922   0.2880000   0.7740039
%    0.9940000   0.1710000   0.7860667
%    0.9640039   0.4330275   0.7850236
%    0.9522941   0.6200000   0.8070000
%    0.8370275   0.8809882   0.7850078
%    0.7639843   0.8209961   0.7690000
%    0.6060745   0.8870000   0.7860039
%    0.5089255   0.8800118   0.7910078
%    0.3110667   0.6483294   0.7870000
%    0.0470000   0.4040118   0.7730039
%    0.1229922   0.3480039   0.7890236
%    0.0390000   0.7349608   0.8190078
%    0.4199647   0.3059961   0.7790000
%    0.9092589   0.1290314   0.8080000
%    1.0030078   0.9061177   0.8120000
%    0.1649882   0.0069490   0.8250118
%    0.6809686   0.8320549   0.8180039
%    0.4909255   0.5343647   0.8130000
%    0.2120196   0.5621725   0.8310039
%    0.1010235   0.6020039   0.8240157
%    0.0810000   0.7650196   0.1320000
%    0.1090000   0.7130157   0.1790157
%    1.0010000   0.6859882   0.1909843
%    0.9880000   0.7600745   0.2420000
%    0.9160471   0.8929922   0.1909961
%    0.9570000   0.9560667   0.2459843
%    0.3018902   0.4459843   0.2210118
%    0.2419961   0.5000000   0.2510000
%    0.1509961   0.2123177   0.2410000
%    0.1250039   0.1830078   0.2819216
%    0.1339922   0.4930000   0.2589922
%    0.1139765   0.5590000   0.3140823
%    0.7660000   0.2189843   0.3200078
%    0.7080000   0.2667098   0.3290000
%    0.9310039   0.5039216   0.3329686
%    0.9920000   0.5412745   0.3669568
%    0.3361647   0.1240000   0.3569922
%    0.3259255   0.2030039   0.3820000
%    0.6739843   0.2999922   0.4750000
%    0.7359686   0.3610000   0.4820235
%    0.3350000   0.8660000   0.4430432
%    0.3349961   0.7970000   0.4772431
%    0.3470823   0.9200000   0.5840118
%    0.3869922   0.9900000   0.6080078
%    0.1450392   0.6920039   0.5819725
%    0.2350000   0.6920078   0.5851647
%    0.1869647   0.7810118   0.6030000
%    0.2340000   0.8500000   0.6340078
%    0.7581098   0.0509961   0.6020000
%    0.7661530   0.0030000   0.6419725
%    0.6849961   0.5770000   0.6398196
%    0.6230314   0.5630000   0.6730000
%    0.7250000   0.9540196   0.6629686
%    0.6529961   0.9220000   0.6919529
%    0.5389961   0.5270039   0.6548588
%    0.5390000   0.5599922   0.7287686
%    0.6090236   0.5000196   0.7560039
%    0.6180000   0.5540000   0.7840784
%    0.8600275   0.7160118   0.8130039
%    0.8830000   0.7040157   0.8220000
%    0.8250000   0.6359922   0.8060157
%    0.8000000   0.5610000   0.8200196
%    1.0220039   0.5363019   0.8230000
%    0.9510157   0.5190039   0.8230078];

clear job
clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

for j = 1:numel(npRedo)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue
    
    particleRadius = ((D(npRedo(j))/2)+(2*w(npRedo(j)))); %this is the apperant particle radius
    [v,center] = getCutoutSectVertices(pksNewer(npRedo(j),:),particleRadius,size(InormRed)); %vertices of rectangular portion of 'Inorm' that is (2*Dactual)^3 and centered on validPks(npRedo(j),:) 
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    centerIdx = sub2ind(size(cutoutSect),center(1),center(2),center(3));
    createTask(job(j),@getPositionSubPixel,2,{v(1,:),centerIdx,K(npRedo(j)),D(npRedo(j)),w(npRedo(j)),offset(npRedo(j)),cutoutSect,rangeColormap});
    disp(sprintf('locating pksNewer(%d,:) at SUB-pixel resolution...',npRedo(j)));
end

% submit one job at a time to be computed by all available slave nodes
for j = 1:length(job)
    submit(job(j));
    disp(sprintf('job%02d, submitted',j))
end 

for j = 1:length(job)
    state = job(j).State;
    
    if(state(1:3) == 'fin')%'fin' for finished 
        try 
            results = fetchOutputs(job(j));
            disp(sprintf('\n------------job%02d, COMPLETED with the following tasks------------',j));

            pksBest(npRedo(j),:) = results{1}; %minus 1 becuase MATLAB indexing starts at 1
            Rsquared(:,npRedo(j)) = results{2};
            disp(sprintf('\t pksNewer(%d,:) becomes pksBest(%d,:)',npRedo(j),npRedo(j)));
        catch err 
            disp(sprintf('\n------------job%02d, FAILED with the following tasks------------',j));
            
            pksBest(npRedo(j),:) = NaN; %minus 1 becuase MATLAB indexing starts at 1
            Rsquared(:,npRedo(j)) = NaN;
            disp(sprintf('\t pksBest(%d,:) = NaN',npRedo(j)));
        end
    end
end

pksBest =(1.0e+03)*[

   0.0539843   0.1150039   0.2449922
   0.1574941   0.1100275   0.2299922
   0.3559882   0.1049804   0.2609961
   0.3690823   0.2040039   0.2469843
   0.4379216   0.1506980   0.1940000
   0.7540627   0.1419961   0.2090000
   0.8700000   0.0464824   0.2390196
   0.7283490   0.2322980   0.2420000
   1.0010000   0.1929529   0.1390000
   0.9970000   0.0993765   0.1889961
   0.0150000   0.3439568   0.1350000
   0.4080000   0.4479922   0.2538667
   0.4545451   0.3150039   0.2490000
   0.5959177   0.5090000   0.2479961
   0.6805020   0.3570157   0.2110000
   0.8050078   0.2948235   0.2130118
   0.9450588   0.5380510   0.2302980
   0.9380078   0.3150039   0.2160118
   0.0190000   0.6970118   0.2182392
   0.1868588   0.6380706   0.2270000
   0.3609804   0.5308000   0.2379961
   0.2390000   0.5929529   0.3259961
   0.6200000   0.6640000   0.2825843
   0.8039882   0.5789961   0.2039569
   0.7590000   0.7480432   0.2289725
   0.0280000   0.8109177   0.2100118
   0.3465882   0.9950000   0.1370000
   0.2490039   0.8299922   0.1710784
   0.2210000   0.9409961   0.2480275
   0.4869922   0.8930000   0.2199882
   0.6230784   0.6900039   0.1830275
   0.7390000   0.8440314   0.2569177
   0.8489843   0.9760039   0.2439725
   0.4589843   0.9570157   0.2859765
   0.5510000   0.9703530   0.2250078
   0.4557333   0.0759686   0.2800000
   0.6349961   0.0019961   0.2409255
   0.9579882   0.1540000   0.2530118
   0.0509843   0.4870039   0.2190275
   0.7632549   0.7620392   0.3080118
   0.8429922   0.9850000   0.1400078
   0.0389961   0.0599294   0.3709451
   0.1048549   0.2470039   0.3579961
   0.1199922   0.0060000   0.3940510
   0.3119882   0.0349922   0.3790392
   0.4160000   0.1660039   0.3120039
   0.4682510   0.0019961   0.3280078
   0.6861921   0.0910000   0.3690078
   0.4860000   0.0519373   0.4620039
   0.5750039   0.0129725   0.4560157
   0.7490000   0.0239961   0.4008549
   0.6469961   0.0990196   0.4569686
   0.7729961   0.0589961   0.2629529
   0.7969922   0.0920078   0.3568784
   0.7206706   0.1410275   0.2900157
   0.8529804   0.0820000   0.4560275
   0.9399255   0.2430078   0.4191020
   0.9279961   0.1370000   0.4899451
   0.1931569   0.2609922   0.3139804
   0.1479882   0.5017686   0.3730039
   0.1419608   0.3759177   0.4840118
   0.0889961   0.4460157   0.4500196
   0.4660196   0.5060039   0.3861647
   0.3889882   0.3540000   0.3107255
   0.3095882   0.4240078   0.4620000
   0.4349961   0.2560157   0.3270275
   0.7200000   0.3769647   0.3537843
   0.6090000   0.4566471   0.4440000
   0.8500078   0.3280078   0.3099568
   0.8119922   0.4150078   0.3268863
   1.0230000   0.3569647   0.3539686
   0.9630157   0.4279529   0.3850000
   0.0599137   0.7150039   0.3559882
   0.1610431   0.6460196   0.3160000
   0.2480667   0.5672353   0.4819961
   0.1278745   0.6950000   0.4766196
   0.3750039   0.6020196   0.4820980
   0.4710039   0.6064980   0.4250118
   0.6356706   0.6914510   0.3760000
   0.7110118   0.6680000   0.3140196
   0.6849647   0.4380000   0.5050039
   0.6260549   0.6969961   0.4699961
   0.8079882   0.4730078   0.4085020
   0.7290196   0.5480236   0.4500000
   0.9089961   0.6830157   0.4849882
   1.0040039   0.6140000   0.2691843
   0.9189961   0.5380118   0.5259804
   0.9980000   0.6580784   0.4506941
   0.1109843   0.8329804   0.3560039
   0.2669765   0.9460157   0.3669804
   0.1210000   0.7790039   0.5339647
   0.1179922   0.8680039   0.4648510
   0.1899882   0.9100000   0.5309529
   0.3780000   0.8377686   0.3580196
   0.3689961   0.7034039   0.4729922
   0.4490235   0.8470000   0.3149961
   0.4629568   0.9390078   0.3680039
   0.6908823   0.9950000   0.3709882
   0.7160000   0.9300275   0.4590902
   0.9579961   0.8140078   0.3569882
   0.8169804   0.8889961   0.3019961
   0.8319686   0.7870000   0.5119843
   0.9169568   0.8030078   0.4350039
   0.8709961   0.8650157   0.3690196
   1.0020000   0.8451216   0.4349804
   0.9416981   0.8600000   0.5119373
   0.1759765   0.9930039   0.3759804
   0.1349922   0.9650039   0.4570039
   0.2340078   0.9660118   0.4539725
   0.6120745   0.9470157   0.3410000
   0.8729490   0.9550039   0.3370000
   0.9060000   0.9520039   0.4479059
   0.1630000   0.1340118   0.6349294
   0.0639922   0.1070157   0.6521451
   0.3693726   0.1070000   0.5470000
   0.2560000   0.1581294   0.5978667
   0.4489804   0.2549961   0.6289961
   0.2532353   0.0299961   0.6660157
   0.4599490   0.1010078   0.6919961
   0.3370157   0.1290000   0.6609882
   0.4247372   0.2030000   0.7670627
   0.5490627   0.1020000   0.5110000
   0.7390157   0.0852157   0.4829922
   0.5448079   0.1370471   0.6000039
   0.5614353   0.1165020   0.6946627
   0.6950196   0.0919961   0.6839882
   0.8290196   0.0149961   0.5629177
   0.8109725   0.1089961   0.5509764
   0.8649922   0.2120000   0.5710353
   0.8879059   0.0739961   0.7360000
   0.8830000   0.2710039   0.6894078
   1.0020078   0.0910157   0.5991843
   0.9457961   0.2159804   0.6461294
   0.9507647   0.2524627   0.7810078
   0.1490000   0.3120000   0.5719647
   0.2360000   0.3220000   0.4879843
   0.0390000   0.3159765   0.7399961
   0.2099059   0.3180118   0.7380000
   0.4287843   0.3610000   0.6029961
   0.3480000   0.4450039   0.5713490
   0.2849922   0.3919961   0.7355176
   0.4999882   0.2650078   0.7439922
   0.3599725   0.4020000   0.6628981
   0.3510000   0.4529608   0.7440157
   0.4500432   0.4254941   0.6749882
   0.6098941   0.3530000   0.6150627
   0.5179843   0.4599804   0.4689765
   0.5735412   0.4570588   0.5480039
   0.6430118   0.3180118   0.7589098
   0.7309725   0.2429922   0.6889922
   0.6940000   0.4959529   0.7010392
   0.7900196   0.4995608   0.5043569
   0.8599137   0.4710157   0.5769961
   0.9239843   0.2930000   0.5909922
   0.9940000   0.3670000   0.5659490
   0.1782314   0.4750039   0.4560039
   0.1909569   0.4080078   0.7019961
   0.1000000   0.7229922   0.7070236
   0.2810000   0.5670902   0.6195686
   0.3666510   0.5550039   0.6210000
   0.3189961   0.5321843   0.7043726
   0.3630039   0.6170039   0.6880823
   0.4439843   0.6880078   0.6889922
   0.4789882   0.5485647   0.5949882
   0.6201530   0.6430078   0.6100275
   0.5630667   0.6679255   0.5439764
   0.5273137   0.6630000   0.6430588
   0.8440275   0.6670000   0.5594588
   0.9179137   0.6740000   0.7310196
   0.9630157   0.6370000   0.5597098
   0.9800353   0.5330000   0.6970823
   0.9880000   0.7069922   0.6634785
   0.1320000   0.8170000   0.7160980
   0.0379961   0.8730157   0.7160000
   0.1981059   0.9629922   0.6130000
   0.3903372   0.6755137   0.6010000
   0.2730235   0.9600078   0.5420078
   0.4300000   0.8010392   0.7650118
   0.5436196   0.9560000   0.6250000
   0.7249686   0.9300000   0.5590314
   0.5990000   0.7238823   0.6760549
   0.5469961   0.9060196   0.7024785
   0.8009608   0.8650118   0.5550706
   0.8109882   0.7660000   0.6029804
   0.8090000   0.9370510   0.7099961
   0.9469843   0.9359961   0.5769961
   0.1090000   0.9463961   0.5860078
   0.0860039   0.0760431   0.8320000
   0.2740000   0.0780000   0.8260078
   0.4809882   0.0319882   0.7500196
   0.6250196   0.0890000   0.7600000
   0.7439961   0.1480157   0.7280000
   0.7699804   0.2310157   0.7630000
   0.9099882   0.1799608   0.7299961
   0.9920000   0.0769843   0.8120275
   0.5739961   0.3750000   0.8090314
   0.6699922   0.2299059   0.7930000
   0.5260000   0.4610078   0.7370000
   0.7680039   0.4990000   0.7659529
   0.9360314   0.3480000   0.7550078
   0.1460510   0.4420039   0.7760314
   0.1980000   0.7509412   0.6940039
   0.2670000   0.6160196   0.7058588
   0.4070039   0.4737882   0.8177568
   0.4321608   0.7069882   0.8320000
   0.5050196   0.7400000   0.7371608
   0.5670000   0.6689922   0.7710157
   0.2120078   0.7647961   0.7860000
   0.6310196   0.7459725   0.7660000
   0.6060039   0.8220588   0.7039961
   0.8880039   0.8059804   0.7675372
   0.9760000   0.7300078   0.8230432
   0.1114824   0.9144941   0.7620078
   0.1620039   0.9400236   0.8300039
   0.4830039   0.9599686   0.8320078
   0.7299961   0.9190039   0.8100000
   0.8000118   0.9870039   0.8260000
   0.3829882   0.5490784   0.4070000
   0.2550275   0.5264235   0.3939961
   0.3582392   0.3752118   0.3900078
   0.4299216   0.3066941   0.4120000
   0.5170510   0.2930039   0.3630000
   0.6919843   0.1920314   0.3850196
   0.9220078   0.0269882   0.4320118
   0.3918510   0.0150000   0.4369961
   0.0389961   0.2454471   0.4230000
   0.1063020   0.3150000   0.4228039
   0.0890078   0.6330039   0.4089608
   0.4389882   0.4274196   0.4190039
   0.6316078   0.5670078   0.4530039
   0.6210157   0.8890000   0.4320039
   0.5519961   0.9500000   0.4139608
   0.1759922   0.0559882   0.1400157
   0.0631412   0.0209961   0.1460078
   0.1720314   0.1780000   0.1390000
   0.2986000   0.2799804   0.1440000
   0.3529686   0.3850000   0.1430118
   0.5148784   0.3678039   0.1430000
   0.5940000   0.2285020   0.1500000
   0.6107608   0.0549843   0.1430039
   0.3347490   0.6800118   0.1480039
   0.3820353   0.5770000   0.1519765
   0.5279922   0.6750196   0.1460275
   0.6019843   0.7810000   0.1460549
   1.0040118   0.7860471   0.1420118
   0.9810745   0.9530000   0.1460000
   0.2050000   0.9050157   0.1401333
   0.4300000   0.9010236   0.1408235
   0.5700000   0.8909804   0.1489882
   0.1271137   0.3400039   0.3299765
   0.2420235   0.6940078   0.3459922
   0.3779961   0.9680039   0.3374745
   0.1859529   0.9019922   0.3959882
   0.2670000   0.8180275   0.3769647
   0.3349922   0.7499686   0.3769843
   0.4170941   0.6939764   0.3660000
   0.5400236   0.4330078   0.3785294
   0.5350000   0.0400039   0.3909098
   0.9229961   0.3389373   0.4040000
   0.9589922   0.9100000   0.3736588
   1.0110000   0.0910078   0.3136039
   0.7749882   0.5000078   0.3189922
   0.8809137   0.0339961   0.3549843
   0.5600000   0.0699961   0.2999922
   0.4049568   0.5220078   0.3170039
   0.4519961   0.4220000   0.3320157
   0.5308196   0.5373373   0.3200000
   0.0110000   0.8039764   0.3766039
   0.4140000   0.0669922   0.3673569
   0.7790000   0.9450039   0.3659765
   0.5320118   1.0140078   0.3214039
   0.6990000   0.8977177   0.3409961
   0.9629647   0.6547294   0.3509882
   0.6130353   0.7890000   0.4450000
   0.4200432   0.7729961   0.4349804
   0.4370000   0.8670000   0.4369686
   0.5242078   0.1350000   0.4200000
   0.2860039   0.6270000   0.4065216
   0.1109529   0.5400039   0.4450039
   0.1386314   0.2250000   0.4410432
   1.0080039   0.9572157   0.4560314
   1.0020236   0.5639882   0.4790353
   0.7769764   0.9950000   0.4959725
   0.4100706   1.0080039   0.4849294
   0.5420000   0.8459216   0.4790039
   0.6499647   0.9740118   0.5140000
   1.0140000   0.0970471   0.4990039
   0.9630118   0.4100000   0.4709765
   0.0109843   0.6650549   0.5179961
   0.2056431   0.1470078   0.5109961
   0.4549686   0.1420078   0.4840157
   0.6099922   0.3779961   0.4989686
   0.8589647   0.3860078   0.4800000
   0.8790000   0.2949177   0.5070745
   0.7200000   0.8429764   0.4970275
   0.2700000   0.7140667   0.4886588
   0.4019922   0.2400196   0.5529608
   0.4960039   0.2485608   0.5450000
   0.0450000   0.8890196   0.5437333
   0.2390000   0.8180000   0.5155726
   0.5080000   0.7790000   0.5444588
   0.5560078   0.5709843   0.5299725
   0.7059922   0.5389804   0.5380314
   0.6760000   0.6250039   0.5339961
   0.3171765   0.5360078   0.5380000
   0.9970000   0.4820000   0.5340000
   0.9250118   0.7630078   0.5404078
   0.6254863   0.1892784   0.5520039
   0.8999686   1.0160000   0.5399608
   0.4429961   0.4731961   0.5368196
   0.8250000   0.5740000   0.5583176
   0.0669922   0.2150392   0.5909843
   0.4220078   0.1619922   0.6010078
   0.6039882   0.0220000   0.5896000
   0.9079764   0.0630196   0.6039961
   0.9870549   0.5495137   0.5940039
   0.6860000   0.8233726   0.5970000
   0.3229804   0.7390118   0.5980078
   0.5639804   0.2450000   0.6141961
   0.8538000   0.9205098   0.6170000
   0.1644549   0.2360000   0.6370314
   0.1639961   0.4810275   0.6409843
   0.6269961   1.0110000   0.6501137
   0.8940392   0.8189922   0.6170118
   0.9000471   1.0030000   0.6334196
   1.0080000   0.8730471   0.6299843
   0.0049961   0.2790039   0.6547372
   0.5050118   0.3200039   0.6539922
   0.9349804   0.0069922   0.6700392
   0.9840000   0.2949608   0.6829804
   0.8870000   0.5070196   0.6609882
   0.5139725   0.8159529   0.6679961
   0.4649686   0.5910510   0.6700000
   0.5529686   0.4386784   0.6489922
   0.8014353   0.6749373   0.6480000
   0.8330000   0.3760196   0.8130353
   0.3075059   0.0450000   0.7290157
   0.3750275   0.1140039   0.7430118
   0.3340000   0.2060353   0.7427608
   0.3009490   0.2939608   0.7210000
   0.2109804   0.0640000   0.7487372
   0.1139922   0.0430157   0.7450000
   0.0639961   0.1520000   0.7449686
   0.0480157   0.6390039   0.7470510
   0.0130000   0.0719961   0.7595412
   0.8880039   0.5780392   0.7390039
   0.7310078   0.9980000   0.7440510
   0.6819961   0.0099922   0.7730039
   0.8930078   0.9720275   0.7270000
   1.0160118   0.3390000   0.8030275
   0.4999922   0.3002588   0.8220078
   0.3080314   0.3354471   0.8050000
   0.2140000   0.3750196   0.8191725
   0.5649804   1.0150118   0.8260078
   0.3870118   0.0179882   0.7940039
   0.8049961   0.0299882   0.7800196
   0.8149961   0.1290353   0.7900157
   0.8567922   0.2880000   0.7740039
   0.9940000   0.1710000   0.7860667
   0.9640039   0.4330275   0.7850236
   0.9522941   0.6200000   0.8070000
   0.8370275   0.8809882   0.7850078
   0.7639843   0.8209961   0.7690000
   0.6060745   0.8870000   0.7860039
   0.5089255   0.8800118   0.7910078
   0.3110667   0.6483294   0.7870000
   0.0470000   0.4040118   0.7730039
   0.1229922   0.3480039   0.7890236
   0.0390000   0.7349608   0.8190078
   0.4199647   0.3059961   0.7790000
   0.9092589   0.1290314   0.8080000
   1.0030078   0.9061177   0.8120000
   0.1649882   0.0069490   0.8250118
   0.6809686   0.8320549   0.8180039
   0.4909255   0.5343647   0.8130000
   0.2120196   0.5621725   0.8310039
   0.1010235   0.6020039   0.8240157
   0.0810000   0.7650196   0.1320000
   0.1090000   0.7130157   0.1790157
   1.0010000   0.6859882   0.1909843
   0.9880000   0.7600745   0.2420000
   0.9160471   0.8929922   0.1909961
   0.9570000   0.9560667   0.2459843
   0.3018902   0.4459843   0.2210118
   0.2419961   0.5000000   0.2510000
   0.1509961   0.2123177   0.2410000
   0.1250039   0.1830078   0.2819216
   0.1339922   0.4930000   0.2589922
   0.1139765   0.5590000   0.3140823
   0.7660000   0.2189843   0.3200078
   0.7080000   0.2667098   0.3290000
   0.9310039   0.5039216   0.3329686
   0.9920000   0.5412745   0.3669568
   0.3361647   0.1240000   0.3569922
   0.3259255   0.2030039   0.3820000
   0.6739843   0.2999922   0.4750000
   0.7359686   0.3610000   0.4820235
   0.3350000   0.8660000   0.4430432
   0.3349961   0.7970000   0.4772431
   0.3470823   0.9200000   0.5840118
   0.3869922   0.9900000   0.6080078
   0.1450392   0.6920039   0.5819725
   0.2350000   0.6920078   0.5851647
   0.1869647   0.7810118   0.6030000
   0.2340000   0.8500000   0.6340078
   0.7581098   0.0509961   0.6020000
   0.7661530   0.0030000   0.6419725
   0.6849961   0.5770000   0.6398196
   0.6230314   0.5630000   0.6730000
   0.7250000   0.9540196   0.6629686
   0.6529961   0.9220000   0.6919529
   0.5389961   0.5270039   0.6548588
   0.5390000   0.5599922   0.7287686
   0.6090236   0.5000196   0.7560039
   0.6180000   0.5540000   0.7840784
   0.8600275   0.7160118   0.8130039
   0.8830000   0.7040157   0.8220000
   0.8250000   0.6359922   0.8060157
   0.8000000   0.5610000   0.8200196
   1.0220039   0.5363019   0.8230000
   0.9510157   0.5190039   0.8230078];

npDimer = 378:size(pksBest,1);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single1stPass12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for np = 1:(npDimer(1)-1)
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest(np,1),pksBest(np,2),pksBest(np,3),K(np),D(np),w(np),offset(np));
end
fclose(fid);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for np = npDimer(1):2:npDimer(numel(npDimer))
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest(np,1),pksBest(np,2),pksBest(np,3),K(np),D(np),w(np),offset(np));
    fprintf(fid,'%d %f %f %f %f %f %f %f\n\n',np+1,pksBest(np+1,1),pksBest(np+1,2),pksBest(np+1,3),K(np+1),D(np+1),w(np+1),offset(np+1));
end
fclose(fid);


% ------(start: remove particles successfully located from Inorm)---------

% [Kactual, Dactual, wActual, offsetActual]
% 
% ans =
% 
%    0.124835500000000  70.564804100000003  12.486556100000000   0.586136000000000


% [min(K(:)), min(D(:)), min(w(:)), min(offset(:))]
% 
% ans =
% 
%    0.082475341856480  53.766948699951172   9.410371780395508   0.019292334094644
   
Kactual = min(K(:));
Dactual = min(D(:));
wActual = min(w(:));
offsetActual = min(offset(:));

clearvars -except InormRed Kactual Dactual wActual offsetActual pksBest K D w offset

cd('~/poincareProgs/')

% run results6may2015_2_2626nm_8xS5 %load pksBest, K, D, w, and offset from all iterations of your particle tracker

tic;
[voronoiVol, Rmask] = peakPlacementB(pksBest,size(InormRed),D,w); %this MUST be done in series because you're editing a single large data set
elapsedTime1 = toc; %2.433391413000000e+03seconds =  40.556523549999994min for 651particles

tic;
Inorm2 = InormRed;
for np = 100:size(pksBest,1) 
    %particles are not drawn in order i.e., pksBest(1,:), pksBest(2,:),
    %etc. instead particles are drawn in the order provided by voronoiVol 
    
    sphericalRegion = find(voronoiVol == np);
  
    Inorm2(sphericalRegion) = InormRed(sphericalRegion)-ssf([K(np), D(np), w(np), offset(np)],Rmask(sphericalRegion));
    disp(sprintf('removed particle %d of %d\n',np, size(pksBest,1)))
end
elapsedTime2 = toc; %1.691247193000000e+03seconds =  28.187453216666665min

figure 
for slice = 5:size(Inorm2,3)
    hold off
    subplot(1,2,1)
    simage(Inorm2(:,:,slice));
    pksToPlot = find(round(pksBest(:,3)) == slice);
    if(~isempty(pksToPlot))
        hold on
        plot(pksBest(pksToPlot,2),pksBest(pksToPlot,1),'k*','markersize',20)
    end
    title(sprintf('Inorm2(:,:,%d)',slice),'fontsize',20);
    
    
    hold off
    subplot(1,2,2)
    simage(InormRed(:,:,slice));
    if(~isempty(pksToPlot))
        hold on
        plot(pksBest(pksToPlot,2),pksBest(pksToPlot,1),'k*','markersize',20)
    end
    title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20);
      
    pause
end

clear InormRed
InormRed = Inorm2;
% ------(END: remove particles successfully located from Inorm)---------

% --------------------(end: particle trakcing 1st Pass)-----------------


% --------------------(start: particle trakcing 2nd Pass)-----------------

% ------------(start: convolution of Inorm2 & maskExp)--------------------

% --------------------(start: thresholding)-------------------------
clearvars -except Inorm* Kactual Dactual wActual offsetActual

ss=2*fix((Dactual/2)+(2*wActual))-1;  % size of ideal particle image NOT the size of the particle in the image
os=(ss-1)/2; % (size-1)/2 of ideal particle image
center = [os+1 os+1 os+1]; %center of ideal particle image
[voronoiVol, Rmask] = peakPlacementB(center,[ss ss ss],Dactual,wActual);
maskExp = getCalcImg(Rmask,Kactual,Dactual,wActual,offsetActual,voronoiVol);%the 'actual' maskImg to be use in particleTracking algorithm

N = size(InormRed) + size(maskExp) - 1;
a = single(fftn(InormRed,N));
b = single(fftn(maskExp,N));

prod = @(x,y) x.*y;

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

numTasks = 100; %one slice == one task
tasksPerJob = 1:numTasks:N(3);
for j = 1:length(tasksPerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue

    startSlice = tasksPerJob(j);
    stopSlice = startSlice+numTasks-1; %subtract to account for MATLAB indexing displacement
   
    if(stopSlice > N(3)) 
        stopSlice = N(3);
    end
    
    for slice = startSlice:stopSlice %slice == task #
        createTask(job(j),prod,1,{a(:,:,slice),b(:,:,slice)});
        disp(sprintf('job%02d,slice%03d',j,slice))
    end
end

% submit one job at a time to be computed by all available slave nodes
clear a b %first clear RAM of 'a' and 'b' to make room for computation of 'c'
for j = 1:length(job)
    submit(job(j)); %~3min for 61 slices
    disp(sprintf('job%02d, submitted',j))
end 

% populate 'c' with the convolution matrix data and delete jobs from memory
cFourier = zeros(N,'single');
for j = 1:length(job)
    state = job(j).State;
    
    startSlice = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        results = fetchOutputs(job(j));
        for slice = 1:size(results,1)
            cFourier(:,:,startSlice+slice-1) = results{slice}; %minus 1 becuase MATLAB indexing starts at 1
            disp(sprintf('c(:,:,%d)',startSlice+slice-1))
        end
    end
end

c = ifftn(cFourier,N);

clear cFourier

RowColStart = ((size(maskExp,1)-1)/2)+1;%plus 1 to correct for matlab indexing/counting which starts at 1 NOT 0 
RowColStop = N(1)-(RowColStart-1); %minus 1 to correct for matlab indexing/counting which starts at 1 NOT 0 

SliceStart = RowColStart;
SliceStop = N(3)-(RowColStart-1);

cScaled = c(RowColStart:RowColStop,RowColStart:RowColStop,SliceStart:SliceStop);

% ------------(stop: convolution of Inorm & maskExp)--------------------


% --------------------(start: thresholding)-------------------------

clearvars -except cScaled Inorm* maskExp Kactual Dactual wActual offsetActual 

sizeSect = round(size(maskExp,1)); 
rangeRow = 1:sizeSect:size(InormRed,1); %size(InormRed,1) = 1024pixels
rangeCol = 1:sizeSect:size(InormRed,2); %size(InormRed,1) = 1024pixels
rangeSlice = 1:sizeSect:size(InormRed,3); %size(InormRed,3) = 795pixels

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

counter = 0;
for slice = 1:length(rangeSlice)
    job(slice) = createJob(clusterD);
    disp(sprintf('\n------------job%02d, created with the following tasks------------',slice));
    for col = 1:length(rangeCol)  
        for row = 1:length(rangeRow)
            stopRow = rangeRow(row)+sizeSect-1;
            stopCol = rangeCol(col)+sizeSect-1;
            stopSlice = rangeSlice(slice)+sizeSect-1;
            
            if(stopRow > size(InormRed,1))
                stopRow = size(InormRed,1);
            end
            
            if(stopCol > size(InormRed,2))
                stopCol = size(InormRed,2);
            end
            
            if(stopSlice > size(InormRed,3))
                stopSlice = size(InormRed,3);
            end
            counter = counter + 1;
            disp(sprintf('thresholding...section %d, InormRed(%d:%d,%d:%d,%d:%d)',counter,rangeRow(row),stopRow,rangeCol(col),stopCol,rangeSlice(slice),stopSlice))
            cutoutSect = cScaled(rangeRow(row):stopRow,rangeCol(col):stopCol,rangeSlice(slice):stopSlice);
            createTask(job(slice),@threshold,1,{cutoutSect,0.999,size(InormRed),[rangeRow(row), rangeCol(col), rangeSlice(slice)]});
        end
    end
end
clear maskExp

% submit one job at a time to computed by all available slave nodes
for num = 1:length(job)
    submit(job(num));
    disp(sprintf('job%02d, submitted',num));
end

% populate tPks...
numJob = 1;
numSect = 1;
pks = fetchOutputs(job(numJob));
thresholdPks = pks{numSect};
disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
for numSect = 2:length(pks)
    thresholdPks = vertcat(thresholdPks,pks{numSect});
    disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
end

for numJob = 1:length(job)
    pks = fetchOutputs(job(numJob));
    for numSect = 1:length(pks)
        thresholdPks = vertcat(thresholdPks,pks{numSect});
        disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
    end
end
thresholdPks = single(thresholdPks);
[tPks(:,1), tPks(:,2), tPks(:,3)] = ind2sub(size(InormRed),thresholdPks);
% --------------------(end: thresholding)-------------------------




% ----(start: eliminate mulitple peaks identifying the same particle)----
clearvars -except tPks sizeSect Inorm* Kactual Dactual wActual offsetActual 

% sizeSect = size(maskExp,1);
numBoxes = 0.8; %OVERESTIMATE THE NUMBER OF PARTICLES BY MAKING THE number of boxes with particles smaller (i.e., sizeSect) smaller so that there are more boxes!  round(size(maskExp,1)/2) makes each box half the size of the mask image and thereby assumes that there are 2 particles (mask images) per box

rangeRow = round(1:(numBoxes*sizeSect):size(InormRed,1)); %size(InormRed,1) = 1024pixels
rangeCol = round(1:(numBoxes*sizeSect):size(InormRed,2)); %size(InormRed,1) = 1024pixels
rangeSlice = round(1:(numBoxes*sizeSect):size(InormRed,3)); %size(InormRed,3) = 795pixels

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

counter = 0;

for slice = 1:length(rangeSlice)
    job(slice) = createJob(clusterD); %1024x1024x795pixels == 9x9x7sections for size(maskExp) = [119,119,119]. 1section == 1 task. 9x9sections or 81 sections per job 
    disp(sprintf('\n------------job%02d, created with the following tasks------------',slice));
    for col = 1:length(rangeCol)  
        for row = 1:length(rangeRow)
            stopRow = rangeRow(row)+(sizeSect*numBoxes)-1;
            stopCol = rangeCol(col)+(sizeSect*numBoxes)-1;
            stopSlice = rangeSlice(slice)+(sizeSect*numBoxes)-1;
            
            if(stopRow > size(InormRed,1))
                stopRow = size(InormRed,1);
            end
            
            if(stopCol > size(InormRed,2))
                stopCol = size(InormRed,2);
            end
            
            if(stopSlice > size(InormRed,3))
                stopSlice = size(InormRed,3);
            end
            
            corrSlice = find(tPks(:,3) >= rangeSlice(slice) & tPks(:,3) <= stopSlice);
            corrCol = find(tPks(:,2) >= rangeCol(col) & tPks(:,2) <= stopCol);
            corrRow = find(tPks(:,1) >= rangeRow(row) & tPks(:,1) <= stopRow);
        
            idxPks = intersect(intersect(corrSlice,corrCol),corrRow);%thresholding peaks located inside cScaled(rangeRow(row):stopRow,rangeCol(col):stopCol,rangeSlice(slice):stopSlice)
            
            counter = counter + 1;
            disp(sprintf('locating peaks %0.2fpixels apart inside...LARGE section %d, InormRed(%d:%d,%d:%d,%d:%d)',Dactual,counter,rangeRow(row),stopRow,rangeCol(col),stopCol,rangeSlice(slice),stopSlice));
            createTask(job(slice),@getNeighborsValid,1,{tPks(idxPks,:),Dactual});
        end
    end
end

% submit one job at a time to computed by all available slave nodes
for num = 1:length(job)
    submit(job(num));
    disp(sprintf('job%02d, submitted',num));
end

% populate validPks...
numJob = 1;
numSect = 1;
pks = fetchOutputs(job(numJob));
validPks = pks{numSect};
disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
for numSect = 2:length(pks)
    validPks = vertcat(validPks,pks{numSect});
    disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
end

for numJob = 1:length(job)
    pks = fetchOutputs(job(numJob));
    for numSect = 1:length(pks)
        validPks = vertcat(validPks,pks{numSect});
        disp(sprintf('job %d, task %d (aka section %d) has %d peaks',numJob,numSect,(length(rangeRow)*length(rangeCol)*(numJob-1))+numSect,length(pks{numSect})));
    end
end
validPks = single(validPks);



% validPks =[
% 
%          269          65          45
%          537          45          40
%          713          45          45
%          801          45          45
%          857          45          45
%          981          45          45
%           69          89          45
%           77          88          45
%          265          90          45
%          289          90          45
%          397          81          45
%          449          88          45
%          537          89          45
%          713          83          45
%          841          90          45
%          893          90          45
%          981          90          33
%          353         179          45
%          357         179          45
%          561         178          45
%          593         177          45
%          709         178          45
%          713         178          45
%          781         177          45
%          785         178          45
%           45         251          45
%           89         266          45
%          189         267          42
%          445         239          45
%          497         254          45
%          505         251          45
%          705         240          45
%          713         243          45
%          789         228          45
%          893         267          43
%          981         267          30
%           45         350          44
%          177         346          45
%          217         356          45
%          357         356          44
%          621         335          45
%          877         325          45
%          981         356          45
%           45         363          45
%          221         376          45
%          349         408          45
%          357         406          44
%          525         368          42
%          617         357          45
%          881         357          45
%          105         445          45
%           89         498          45
%          237         479          45
%          345         495          45
%          357         497          45
%          801         445          45
%          981         444          45
%          973         498          45
%           89         502          45
%          349         502          45
%          357         501          45
%          449         512          45
%          625         510          45
%          853         518          45
%          857         514          45
%          965         503          45
%           89         626          45
%          153         606          45
%          389         624          44
%          429         624          45
%          761         620          45
%          789         623          45
%          957         592          45
%          269         674          45
%          533         680          39
%          625         680          45
%          689         712          45
%          801         712          45
%          889         702          45
%          177         735          45
%          269         713          45
%          697         759          45
%          713         760          45
%          889         714          45
%          981         713          37
%           45         800          36
%           45         854          45
%          173         802          44
%          353         819          45
%          357         821          45
%          445         799          45
%          497         825          45
%          501         801          45
%          713         851          45
%           45         866          45
%          269         893          42
%          445         924          36
%          713         856          45
%          757         891          45
%          805         874          45
%          893         862          44
%          981         890          39
%          113         979          45
%          157         979          45
%          213         979          42
%          269         980          30
%          445         931          35
%          497         979          45
%          501         979          45
%          709         980          45
%          713         980          45
%          769         980          45
%          925         979          41
%          929         979          41
%          269          65          45
%          537          45          40
%          713          45          45
%          801          45          45
%          857          45          45
%          981          45          45
%           69          89          45
%           77          88          45
%          265          90          45
%          289          90          45
%          397          81          45
%          449          88          45
%          537          89          45
%          713          83          45
%          841          90          45
%          893          90          45
%          981          90          33
%          353         179          45
%          357         179          45
%          561         178          45
%          593         177          45
%          709         178          45
%          713         178          45
%          781         177          45
%          785         178          45
%           45         251          45
%           89         266          45
%          189         267          42
%          445         239          45
%          497         254          45
%          505         251          45
%          705         240          45
%          713         243          45
%          789         228          45
%          893         267          43
%          981         267          30
%           45         350          44
%          177         346          45
%          217         356          45
%          357         356          44
%          621         335          45
%          877         325          45
%          981         356          45
%           45         363          45
%          221         376          45
%          349         408          45
%          357         406          44
%          525         368          42
%          617         357          45
%          881         357          45
%          105         445          45
%           89         498          45
%          237         479          45
%          345         495          45
%          357         497          45
%          801         445          45
%          981         444          45
%          973         498          45
%           89         502          45
%          349         502          45
%          357         501          45
%          449         512          45
%          625         510          45
%          853         518          45
%          857         514          45
%          965         503          45
%           89         626          45
%          153         606          45
%          389         624          44
%          429         624          45
%          761         620          45
%          789         623          45
%          957         592          45
%          269         674          45
%          533         680          39
%          625         680          45
%          689         712          45
%          801         712          45
%          889         702          45
%          177         735          45
%          269         713          45
%          697         759          45
%          713         760          45
%          889         714          45
%          981         713          37
%           45         800          36
%           45         854          45
%          173         802          44
%          353         819          45
%          357         821          45
%          445         799          45
%          497         825          45
%          501         801          45
%          713         851          45
%           45         866          45
%          269         893          42
%          445         924          36
%          713         856          45
%          757         891          45
%          805         874          45
%          893         862          44
%          981         890          39
%          113         979          45
%          157         979          45
%          213         979          42
%          269         980          30
%          445         931          35
%          497         979          45
%          501         979          45
%          709         980          45
%          713         980          45
%          769         980          45
%          925         979          41
%          929         979          41
%          289          62          90
%          801          45          90
%          849          45         137
%          889          45          90
%          977          33          90
%           65          89          90
%           89          89          90
%          353          90          90
%          401          80          90
%          449          85          90
%          537          89          90
%          713          89          90
%          889          90          90
%          977          90          90
%          537         164          90
%          713         178          90
%          785         178          90
%           89         267          90
%          449         253          90
%          537         264          90
%          705         239          90
%          713         268          90
%          801         221          90
%          889         267          90
%          977         267          90
%          209         356          90
%          353         356          90
%          625         334          90
%          873         326          90
%          977         356          90
%           89         426          90
%          625         357          90
%          889         357          90
%           89         434          90
%           89         498          90
%          233         445          90
%          337         490          90
%          361         489          90
%          441         445          90
%          801         445          90
%          969         497          90
%           89         499          90
%          449         535          90
%          537         558          90
%          625         513          90
%          881         535          90
%          969         499          90
%           89         624          90
%          161         610         129
%          273         624          90
%          721         623          90
%          953         590          90
%          433         712          90
%          625         662          90
%          705         712          90
%          801         710          90
%          889         695          90
%           33         713          90
%          177         727          90
%          705         755          90
%          713         756          90
%          825         740          90
%          905         713          90
%          977         713          90
%          185         802          90
%          353         800          90
%          449         786          90
%          497         817          90
%          537         801          90
%          889         802          90
%           49         883          90
%           89         886          90
%          713         863          90
%          753         869         140
%          801         872          90
%          889         891          90
%          977         874          90
%           49         980          90
%           89         979          90
%          177         991          90
%          449         979          90
%          513         979          90
%          801         980          90
%          921         980          90
%          977         980          90
%          273          68         179
%          289          69         179
%          401          45         179
%          449          45         179
%          529          24         206
%          705          43         179
%          801          45         179
%           97         126         179
%          273          82         178
%          289          72         179
%          385          90         179
%          897          90         179
%           49         178         179
%          273         184         178
%          289         179         178
%          513         178         179
%          625         160         202
%          657         168         179
%          993         178         186
%           49         267         178
%           97         267         179
%          401         278         179
%          481         251         146
%          801         221         179
%          897         225         206
%           65         287         178
%           97         304         178
%          193         356         179
%          769         356         179
%          849         356         179
%          865         356         179
%          177         369         179
%          273         357         213
%          337         364         179
%          897         417         171
%          977         427         202
%           97         446         179
%          177         446         179
%          801         432         179
%          977         430         178
%          273         567         173
%          449         528         179
%          529         534         179
%          993         565         179
%          129         613         179
%          529         573         179
%          657         598         182
%          721         623         179
%          977         624         179
%          209         712         179
%          721         643         179
%          801         713         179
%          881         713         179
%          481         802         179
%          977         800         204
%          209         909         165
%          289         914         178
%          705         890         179
%          753         870         143
%          801         887         179
%          273         980         178
%          481         979         179
%          625         979         179
%          705         941         178
%          721         942         178
%          801         980         189
%           33          24         255
%           97          26         237
%          193          38         268
%          353          12         268
%          993          20         267
%          257          89         268
%          257         143         268
%          289         179         268
%          417         162         273
%          513         178         268
%          641         178         268
%          801         166         267
%           33         267         257
%           97         267         268
%          289         283         242
%          353         284         268
%          561         267         255
%          577         267         257
%          641         267         268
%          993         254         267
%           49         294         228
%           97         356         268
%          289         286         242
%          353         286         268
%          769         356         268
%           49         357         230
%          273         357         215
%          545         377         266
%          577         377         267
%          897         410         268
%          993         415         268
%          705         441         251
%          865         446         268
%          993         446         268
%            1         536         273
%          481         534         268
%          993         535         276
%           33         623         267
%           97         623         268
%          257         623         268
%          353         641         267
%          449         601         245
%          897         623         267
%           97         712         269
%          257         712         268
%          353         642         267
%          513         648         265
%          641         712         268
%          897         672         267
%           97         771         267
%          193         777         267
%          289         713         255
%          353         713         258
%          545         713         258
%          609         771         267
%          129         784         267
%          193         802         267
%          257         802         271
%          609         802         267
%          897         802         278
%           97         914         264
%          257         891         268
%          321         890         227
%          609         891         268
%          897         891         267
%          993         863         267
%           97         980         267
%          321         979         267
%          353         933         253
%          609         987         268
%          705         980         215
%          673         979         268
%          257        1005         267
%          993        1013         267
%          481          22         320
%          609          18         356
%          705          16         324
%          193         140         356
%          897          89         298
%           33         178         342
%          193         143         356
%          801         158         286
%          897         167         356
%          993         177         356
%          353         267         297
%          993         268         304
%          257         356         356
%          897         356         356
%          193         409         356
%          257         401         346
%          609         386         296
%           33         475         342
%           97         446         333
%          257         446         329
%          353         445         356
%          289         446         328
%          609         446         340
%          641         569         356
%           33         623         292
%          545         623         356
%          801         596         356
%          865         623         356
%          545         655         347
%          865         668         356
%          993         712         356
%          193         783         286
%          865         713         356
%          993         738         356
%          193         784         287
%          609         801         329
%          705         799         356
%          545         890         305
%          993         871         287
%           33         979         356
%          545         980         356
%          865         980         345
%          993         979         332
%          257        1004         290
%          609          18         359
%          193          89         357
%          993          89         400
%           33         179         357
%          193         152         394
%          225         171         422
%          449         178         420
%          865         176         396
%          993         143         396
%          641         267         395
%          833         268         427
%          193         356         410
%          257         356         381
%          609         356         357
%          801         356         406
%          897         342         387
%          193         399         393
%          257         396         357
%          609         362         357
%          641         417         357
%          801         363         406
%          897         357         385
%          993         357         357
%           33         487         357
%          193         446         374
%          257         446         357
%          353         445         365
%          641         445         357
%          673         487         357
%           33         535         372
%          609         534         371
%          897         534         421
%          929         514         426
%           33         624         357
%          353         612         357
%          545         623         357
%          609         571         366
%          801         623         404
%          865         623         391
%          545         655         357
%          705         712         417
%          801         712         409
%          897         712         357
%          993         712         371
%          161         765         399
%          257         731         417
%          545         713         364
%          705         769         402
%          801         735         408
%          897         720         357
%          993         740         391
%          161         802         425
%          449         802         420
%          545         843         383
%          609         801         357
%          705         802         392
%          769         804         406
%          801         802         406
%            1         867         357
%          545         891         357
%          897         890         357
%           33         979         372
%           97         932         357
%          225         966         412
%          353         979         415
%          641         979         427
%          865         980         357
%          993         979         357
%          257        1011         357
%           33          42         445
%          225          30         445
%          673          15         445
%          833          10         445
%           33         142         445
%          353          89         445
%          353         139         466
%          385          89         445
%          545          89         445
%           33         143         445
%           33         179         495
%          321         151         486
%          417         200         432
%          449         178         446
%          609         178         445
%          801         178         445
%          897         178         446
%          993         178         446
%          225         220         445
%          257         227         496
%          321         267         445
%          289         224         497
%          545         260         445
%          577         251         445
%          833         268         428
%          993         268         445
%           33         351         445
%          353         288         445
%           33         357         445
%          513         357         445
%          417         488         470
%          705         445         446
%          545         534         445
%          865         534         446
%          993         504         445
%          929         510         446
%          193         615         432
%          545         624         446
%          705         623         446
%          801         623         446
%          897         581         446
%          705         712         446
%          769         703         446
%          801         663         446
%          993         712         446
%           65         782         445
%           97         782         444
%          193         779         446
%          449         769         429
%          769         713         446
%          801         713         446
%          993         752         446
%           65         802         445
%           33         789         487
%           97         802         445
%          193         802         428
%          353         805         445
%          801         854         446
%          801         891         445
%           33         979         446
%          353         979         446
%          641         979         428
%          769         997         446
%          801         997         446
%          513         999         445
%          609        1007         430
%          769        1000         433
%          801        1001         433
%          257          45         534
%          385          24         534
%          449          45         534
%          513          20         534
%          705          24         523
%          897          42         518
%          961          29         525
%           65          89         534
%          257          90         519
%          705         171         534
%          961         178         538
%           65         267         535
%          193         243         513
%          257         268         534
%          641         268         534
%          993         221         502
%          961         268         535
%           65         343         526
%          385         343         537
%          449         356         505
%          513         356         535
%           65         357         534
%          385         357         534
%          449         358         504
%          513         364         534
%          897         420         534
%          961         389         532
%          193         431         534
%          257         476         534
%           65         534         534
%          193         534         534
%           65         624         534
%          193         591         534
%          257         623         534
%          449         623         534
%          449         674         534
%          641         712         564
%          193         726         535
%          385         759         534
%          641         783         534
%          705         713         569
%          641         784         534
%          897         802         534
%            1         891         530
%          385         866         534
%          449         866         534
%          513         890         534
%          641         891         534
%            1         980         534
%          193         979         534
%          449         980         548
%          513         955         514
%          961         979         530
%            1         998         534
%          193        1005         534
%           65          45         580
%          129          34         614
%          193          31         624
%          385          30         624
%          449          45         598
%          513          45         601
%          449          90         593
%          833         130         623
%            1         174         640
%          641         178         623
%          321         267         624
%          385         267         623
%          449         258         613
%          641         262         623
%          769         267         624
%          833         267         623
%           65         356         624
%          257         356         624
%          321         341         575
%          321         314         624
%          385         312         626
%          705         333         623
%          833         356         623
%          897         305         611
%           65         405         574
%          705         425         623
%          833         397         623
%          897         419         623
%           65         446         582
%          257         445         575
%          257         459         627
%          385         475         623
%          449         480         623
%          705         432         618
%          769         445         623
%          833         446         623
%          961         445         623
%          641         531         592
%          769         535         623
%          833         535         623
%           65         623         622
%          193         597         578
%          193         624         625
%          321         623         624
%          769         627         600
%          897         623         623
%          961         623         623
%           65         647         621
%          449         712         624
%          705         658         623
%          897         696         623
%            1         713         574
%           65         782         598
%            1         782         623
%          449         759         608
%          705         713         571
%          897         713         623
%           65         784         599
%          193         851         585
%          385         801         623
%          449         799         596
%          513         873         571
%          641         913         623
%          769         858         623
%            1         980         623
%          257         979         623
%          705         980         596
%          769         979         623
%          833         979         623
%          961         979         623
%          257        1002         623
%          769         999         623
%          833         998         630
%          513          26         653
%          705          71         704
%            1          90         712
%          705          73         704
%          833          89         646
%          961          89         712
%            1         174         642
%            1         180         707
%          193         170         712
%          257         166         712
%          513         178         712
%          641         178         657
%          833         178         655
%          961         179         709
%          641         259         677
%          833         267         644
%          449         353         712
%          705         348         666
%          833         355         666
%          897         356         666
%          961         356         712
%          449         357         712
%          577         387         712
%          641         390         712
%          705         357         671
%          833         396         666
%          897         402         664
%          513         445         712
%          577         458         712
%          769         446         643
%          833         446         652
%          961         432         674
%           65         534         677
%          449         533         711
%          513         546         712
%           65         623         642
%            1         624         712
%          129         623         660
%          641         624         712
%          705         623         712
%          833         573         659
%          897         623         644
%          321         712         712
%          385         712         712
%          897         697         649
%          577         721         704
%          705         713         702
%          769         713         712
%          897         713         656
%          257         801         712
%          321         801         712
%          705         853         684
%          961         801         712
%          257         890         712
%          385         891         712
%          641         919         674
%          897         880         698
%          961         890         712
%           65         979         691
%            1         981         697
%          193         979         699
%          257         980         712
%          385         980         712
%          449         979         712
%          641         980         712
%          961         980         700
%          513        1004         712
%            1          13         713
%          321          45         741
%          961          26         752
%            1         105         713
%          193          89         774
%          385         113         764
%          705          79         713
%          961         107         713
%            1         180         713
%          257         150         732
%          833         178         713
%          961         179         713
%          641         243         713
%          577         268         713
%          833         267         713
%           65         312         713
%          129         317         713
%          449         356         763
%          705         356         713
%          833         327         713
%          897         356         741
%          961         356         731
%          129         357         713
%          449         385         778
%          577         390         717
%          705         377         713
%          897         413         713
%          961         420         713
%            1         445         713
%           65         445         713
%          513         460         726
%          897         446         727
%           65         534         714
%          449         534         713
%          513         547         734
%            1         624         729
%           65         624         713
%          129         613         713
%          193         623         775
%          449         624         780
%          641         624         718
%          833         616         713
%          961         623         725
%            1         657         744
%          321         712         728
%          833         683         723
%           65         720         726
%          321         750         732
%          385         722         724
%          833         722         713
%          257         801         725
%          961         801         721
%          193         890         713
%          257         871         744
%          513         886         766
%          641         890         713
%          705         861         713
%          897         876         713
%          961         887         732
%           65         991         713
%          193         979         713
%          257         980         713
%          449         979         778
%          641         938         713
%          513        1004         751
%           65           5         820
%          385          45         789
%          833          39         789
%          961          16         802
%          193         142         789
%          193          89         802
%          385         139         802
%          513          89         801
%          705         120         798
%          897         133         802
%           65         178         801
%          193         145         789
%          257         196         802
%          385         179         801
%          513         179         789
%          577         179         789
%          641         179         803
%          833         178         816
%          961         178         802
%           65         246         801
%            1         251         818
%          257         268         789
%          577         268         802
%          769         267         815
%          833         215         815
%          897         222         818
%          449         356         802
%          705         356         802
%          961         345         802
%           65         395         814
%          129         368         802
%          449         387         789
%          705         390         802
%          961         357         802
%          257         445         801
%          257         498         791
%          321         445         801
%          513         454         802
%          577         445         802
%          705         446         789
%          833         445         802
%          897         453         802
%           65         533         802
%          257         499         789
%          321         534         801
%          385         534         798
%          833         535         804
%          961         527         813
%            1         624         802
%          193         623         802
%          385         573         789
%          449         607         789
%          577         571         802
%          705         621         789
%          833         624         809
%            1         658         802
%          193         655         793
%          321         712         802
%          577         712         802
%          961         685         802
%          385         759         802
%          513         764         801
%          705         713         801
%          769         713         789
%           65         801         802
%            1         802         802
%          129         801         802
%          257         801         802
%          321         801         802
%          513         802         812
%          769         802         802
%          897         801         820
%          961         801         802
%           65         890         801
%            1         890         816
%          257         859         802
%          321         868         789
%          385         890         789
%          641         868         802
%          897         890         802
%           65         979         789
%            1         979         789
%          257         985         802
%          321         980         789
%          385         980         791
%          385         927         802
%          449         997         785
%          641         979         802
%          705         979         812
%          897         979         801
%          961         980         789
%           65        1005         802
%          449        1000         802
%          641        1006         801
%          705        1003         801];

npDelete = find(validPks(:,3) < 80);

validPks(npDelete,:) = [];

% --(start: move peaks to minimize error and locate at pixel resolution)--

% Unsuccessful at locating validPks found from thresholding using
% 1.5*Dactual.  ALL pksNew from validPks were WRONG!
% 
% Therefore, relocated validPks by hand and saved as
% validPks2ndPass_row_col_slice.txt
% 
% pksNew for validPks2ndPass_row_col_slice.txt using 1.5*Dactual were again
% WRONG!


clearvars -except Inorm* Kactual Dactual wActual offsetActual

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('validPks2ndPass_row_col_slice.txt','r');
temp = fscanf(fid,'%f %f %f\n',[3,Inf]);
fclose(fid);

validPks(:,1) = temp(1,:)';
validPks(:,2) = temp(2,:)';
validPks(:,3) = temp(3,:)';

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers2ndPass_row_col_slice.txt','r');
dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksDimers(:,1) = dimers(1,:)';
pksDimers(:,2) = dimers(2,:)';
pksDimers(:,3) = dimers(3,:)';

npDimer = (size(validPks,1)+1):(size(pksDimers,1)+size(validPks,1));
validPks = vertcat(validPks,pksDimers);

pksNew = validPks;

% % print single particles in pksNew
% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('pksNew2ndPass_id_row_col_slice.txt','w');
% for np = 1:size(pksNew,1)
%     fprintf(fid,'%d %d %d %d\n',np,pksNew(np,1),pksNew(np,2),pksNew(np,3));
% end
% fclose(fid);

clear *pks*

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single2ndPass_id_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
fclose(fid);

order = temp(1,:)';
pksNew(order,1) = temp(2,:)';
pksNew(order,2) = temp(3,:)';
pksNew(order,3) = temp(4,:)';

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers2ndPass_row_col_slice.txt','r');
dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksDimers(:,1) = dimers(1,:)';
pksDimers(:,2) = dimers(2,:)';
pksDimers(:,3) = dimers(3,:)';

npDimer = (size(pksNew,1)+1):(size(pksDimers,1)+size(pksNew,1));
pksNew = vertcat(pksNew,pksDimers);

npGood = 1:size(pksNew,1);
npDelete = find(isnan(pksNew(:,1)) == 1);

npGood(npDelete) = [];

pksNew(npDelete,:) = [];

% [dist, pksOrder] = orderPeaks(pksNew);
% 
% % view multimers:
% for j = 96:numel(npGood)
%     np = npGood(j);
%     npp = pksOrder(1:3,np);
%     
% %     %for 2mers:
% %     npp = union(pksOrder(1:3,np),[np:(np+1)]);
%     
% %     npp = pksOrder(1:3,np);
%     
% %     %for 5mers:
% %     npp = union(pksOrder(1:6,np),[np:(np+4)]);
% 
% %     %for 6mers
% %     npp = union(pksOrder(1:6,np),[np:(np+5)]);
%     
%     clear label 
%     disp(sprintf('---------------------------x-----------------------------'));
%     for jj = 1:numel(npp)
%         label{jj} = num2str(sprintf('%d',npp(jj)));
%         disp(sprintf('\tpksNew(%d,:) = [%.03f, %.03f, %.03f], distance(%d,%d) = %.04f',npp(jj),pksNew(npp(jj),1),pksNew(npp(jj),2),pksNew(npp(jj),3),npp(jj),np,dist(npp(jj),np)));
%     end
%     disp(sprintf('---------------------------x-----------------------------\n'));
%    
%     
%     startSlice = round(min(pksNew(np,3))-10);
%     stopSlice = round(max(pksNew(np,3))+10);
% 
%     if(startSlice < 1)
%         startSlice = 1;
%     end
% 
%     if(stopSlice > size(InormRed,3))
%         stopSlice = size(InormRed,3);
%     end
% 
%     nppClose = npp;
%     nppClose(nppClose == np) = [];
%     
%     for slice = startSlice:stopSlice
%         figure(1)
%         subplot(1,2,1)
%         hold off
%         simage(InormRed(:,:,slice));
%         hold on
%         plot(pksNew(npp,2),pksNew(npp,1),'k*','markersize',20)
%         plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
%         text(pksNew(npp,2),pksNew(npp,1),label,'fontsize',12)
%         title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
%         xlabel(sprintf('pksNew(%d of %d,:) = [%.04f, %.04f, %.04f] \n %d particles plotted',np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3),numel(npp)),'fontsize',20)
% 
%         subplot(1,2,2)
%         hold off
%         simage(InormGreen(:,:,slice));
%         hold on
%         plot(pksNew(npp,2),pksNew(npp,1),'k*','markersize',20)
%         plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
%         text(pksNew(npp,2),pksNew(npp,1),label,'fontsize',12)
%         title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
%         xlabel(sprintf('pksNew(%d of %d,:) = [%.04f, %.04f, %.04f] \npksNew(%d of %d) = [%.04f, %.04f, %.04f]',nppClose(1),size(pksNew,1),pksNew(nppClose(1),1),pksNew(nppClose(1),2),pksNew(nppClose(1),3),nppClose(2),size(pksNew,1),pksNew(nppClose(2),1),pksNew(nppClose(2),2),pksNew(nppClose(2),3)),'fontsize',20)
% 
%         pause
%     end
% end

% -------(start: get all parameters K,D,w,offset for each particle)--------
 
% Kactual = 0.1248355; 
% Dactual = 70.5648041;  
% wActual = 12.4865561;   
% offsetActual = 0.5861360;

Kactual = 0.082475341856480;  
Dactual = 53.766948699951172;
wActual = 9.410371780395508;
offsetActual = 0.019292334094644;

clearvars -except Inorm* Kactual Dactual wActual offsetActual

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single2ndPass_id_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
fclose(fid);

order = temp(1,:)';
pksNew(order,1) = temp(2,:)';
pksNew(order,2) = temp(3,:)';
pksNew(order,3) = temp(4,:)';

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers2ndPass_row_col_slice.txt','r');
dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksDimers(:,1) = dimers(1,:)';
pksDimers(:,2) = dimers(2,:)';
pksDimers(:,3) = dimers(3,:)';

npDimer = (size(pksNew,1)+1):(size(pksDimers,1)+size(pksNew,1));
pksNew = vertcat(pksNew,pksDimers);

npGood = 1:size(pksNew,1);
npDelete = find(isnan(pksNew(:,1)) == 1);

npGood(npDelete) = [];

pksNew(npDelete,:) = [];

numTry = 10; %# of tries to calculate param(np,:)

param = zeros([size(pksNew,1), 4]);
Rsquared = zeros([size(pksNew,1), 1],'single');

tic;
for np = 1:size(pksNew,1)
    [v,center] = getCutoutSectVertices(pksNew(np,:),0.75*Dactual,size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    %bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    %rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    [param(np,:), Rsquared(np)] = getParameters(avgProfile(startIdx:stopIdx),Dactual,wActual,numTry);
    
    if(~isnan(param(np,1)))
        disp(sprintf('SUCCESSFUL computation of parameters for pksNew(%d,:)',np));
    end
    
    if(isnan(param(np,1)))
        disp(sprintf('FAILED computation of parameters for pksNew(%d,:)',np));
    end
end
elapsedTime = toc; %16.5133669seconds for 285 particles

K = param(:,1);
D = param(:,2);
w = param(:,3);
offset = param(:,4);

% param =(1.0e+02)*[
% 
%    0.001626689136028   0.700081710815430   0.130961837768555   0.005748860836029
%    0.004213369488716   0.996036376953125   0.217023582458496   0.003216860592365
%    0.001373026967049   0.702410354614258   0.180403270721436   0.005888956785202
%    0.001379778385162   0.678561096191406   0.156764383316040   0.005711719989777
%    0.001069074049592   0.511003913879394   0.125844030380249   0.005694840550423
%    0.001155533343554   0.616300582885742   0.228901710510254   0.006325969696045
%    0.001620403528214   0.779639053344727   0.120543031692505   0.005939540863037
%    0.003275525569916   0.980249938964844   0.239234390258789   0.004348288476467
%    0.001366440504789   0.651432800292969   0.158616085052490   0.006138064861298
%    0.000561881996691   0.486667098999023   0.099430541992187   0.006916378140450
%    0.001095413714647   0.672352523803711   0.099170017242432   0.005872829556465
%    0.002364481091499   0.772096405029297   0.178349151611328   0.004530478417873
%    0.002185411900282   0.742180023193359   0.172064189910889   0.004815630614758
%    0.006511397957802   1.140401916503906   0.373016357421875   0.000811772122979
%    0.002093155980110   0.750231933593750   0.134220943450928   0.004499850869179
%    0.002113492786884   0.729236373901367   0.173541584014893   0.004519376754761
%    0.001694871932268   0.690377426147461   0.129800052642822   0.005144584178925
%    0.005744611024857   1.156977844238281   0.297900009155273   0.001173442378640
%    0.001809971630573   0.690082168579102   0.171362323760986   0.004847627878189
%    0.002244126200676   0.702211380004883   0.185413589477539   0.004900319278240
%    0.001877402514219   0.667368240356445   0.147517976760864   0.005075790286064
%    0.003459694087505   0.883341674804688   0.287954540252686   0.003308722376823
%    0.002882576584816   0.818693161010742   0.236235466003418   0.004363935887814
%    0.001924621164799   0.622936363220215   0.191988048553467   0.005083823204041
%    0.006041778326035   1.045429458618164   0.181191635131836   0.000934812799096
%    0.004718746244907   0.971140136718750   0.382003440856934   0.002871006727219
%    0.001919619888067   0.647881927490234   0.133355522155762   0.004808759689331
%    0.001624137759209   0.626550521850586   0.141737451553345   0.005088241696358
%    0.005026643872261   1.065180892944336   0.366894989013672   0.001977983117104
%    0.001409883052111   0.476794281005859   0.126459751129150   0.004962096810341
%    0.005860819220543   1.156775207519531   0.362996101379395   0.001343727409840
%    0.001595026105642   0.650788879394531   0.137286481857300   0.005274205207825
%    0.001241949126124   0.590686378479004   0.127111473083496   0.005484029650688
%    0.005335355401039   1.089024658203125   0.274056034088135   0.001630156934261
%    0.004791601598263   1.102586212158203   0.222138309478760   0.001937154829502
%    0.001446709632874   0.643937988281250   0.133395891189575   0.005426278114319
%    0.001420859247446   0.634015998840332   0.140148706436157   0.005381944179535
%    0.002054664045572   0.766882019042969   0.176652355194092   0.004691389203072
%    0.001377572417259   0.573375816345215   0.124520902633667   0.005742106437683
%    0.002445891201496   0.768236541748047   0.172218017578125   0.004643688797951
%    0.002547087073326   0.824421005249023   0.266195907592773   0.004559686779976
%    0.001326843053102   0.654190368652344   0.120268554687500   0.005522387027740
%    0.002065649330616   0.658573455810547   0.142349386215210   0.005135874152184
%    0.004562323391438   1.036813507080078   0.483943214416504   0.003123162090778
%    0.001974675506353   0.735239791870117   0.204902553558350   0.004951849281788
%    0.004871541857719   0.976655426025391   0.252081317901611   0.002327967584133
%    0.001192963272333   0.626847915649414   0.067487297058105   0.005549242496490
%    0.002618892192841   0.763054733276367   0.182832736968994   0.004381961822510
%    0.004168935716152   1.002888259887695   0.225747623443604   0.002653227746487
%    0.002429427951574   0.798900375366211   0.199062385559082   0.004336560964584
%    0.001223559826612   0.601633415222168   0.156105251312256   0.005613782405853
%    0.001580297499895   0.703565216064453   0.121904039382935   0.005642821192741
%    0.001984894275665   0.702149887084961   0.147528181076050   0.004810744822025
%    0.003420042991638   0.876248855590820   0.219333629608154   0.003530635535717
%    0.002870749235153   0.911317291259766   0.228953056335449   0.004004518687725
%    0.001529164761305   0.587966423034668   0.132678098678589   0.005622845888138
%    0.001601684093475   0.701760482788086   0.164174976348877   0.005391453504562
%    0.002112801074982   0.753033523559570   0.190997867584229   0.004844561815262
%    0.003355203568935   0.942418823242188   0.233060150146484   0.003795498311520
%    0.005419316887856   1.118215484619141   0.295230693817139   0.001914349347353
%    0.002620644867420   0.836034317016602   0.165793685913086   0.004093385636806
%    0.003536171913147   0.920922851562500   0.217005519866943   0.003646059334278
%    0.001401395350695   0.536404838562012   0.116088199615479   0.005439224243164
%    0.000939638763666   0.608305740356445   0.091707134246826   0.005720565319061
%    0.001361616849899   0.589052581787109   0.113090114593506   0.005774816274643
%    0.001193960309029   0.527999191284180   0.104667787551880   0.006354433298111
%    0.005626214146614   1.140025024414062   0.280113658905029   0.001741332858801
%    0.002205626368523   0.735335998535156   0.177536983489990   0.004850660264492
%    0.001757955402136   0.599343109130859   0.203318920135498   0.005675869584084
%    0.004269691407681   0.943299865722656   0.346356620788574   0.003406108617783
%    0.003735060691833   1.008851547241211   0.299804859161377   0.003649284839630
%    0.000773053839803   0.312806930541992   0.134984989166260   0.006900937557220
%    0.002478671371937   0.822429580688477   0.179206504821777   0.004524542391300
%    0.001363940984011   0.616501083374023   0.131515703201294   0.006111965775490
%    0.006402861475945   1.243406219482422   0.411938438415527   0.001291826218367
%    0.002724638283253   0.824125976562500   0.256528434753418   0.004684564769268
%    0.001476149111986   0.602739830017090   0.127141647338867   0.005502166152000
%    0.002739206850529   0.729943542480469   0.218589859008789   0.003996579349041
%    0.003141257762909   0.919690322875977   0.243488998413086   0.003758335411549
%    0.000949956849217   0.569686889648438   0.064372725486755   0.005861158967018
%    0.002501482665539   0.787898406982422   0.211962528228760   0.004817315042019
%    0.006941682100296   1.393177032470703   0.464124870300293   0.000802412182093
%    0.001181585043669   0.620869178771973   0.106589927673340   0.005617925524712
%    0.001641951948404   0.619386978149414   0.153954563140869   0.005740472078323
%    0.003535416424274   1.052474670410156   0.344355392456055   0.003868314027786
%    0.004412575960159   0.945795974731445   0.194506320953369   0.002775108814240
%    0.001303784400225   0.627402801513672   0.115606327056885   0.005934365391731
%    0.001428082883358   0.625668449401856   0.082013301849365   0.005278236865997
%    0.003689153194427   0.976111221313477   0.298747272491455   0.003990989625454
%    0.001621634960175   0.603891563415527   0.169248027801514   0.005941306352615
%    0.001059360951185   0.584055328369141   0.093173522949219   0.006229863166809
%    0.001203651204705   0.565502395629883   0.155331296920776   0.006060971021652
%    0.006803213953972   1.394926147460938   0.429790878295898   0.001176348924637
%    0.000811292678118   0.597857818603516   0.071349234580994   0.006653156280518
%    0.002663824856281   0.856084365844727   0.208784294128418   0.004867908358574
%    0.003670360147953   1.061364822387695   0.343649406433105   0.003970317840576
%    0.005566983819008   1.232793426513672   0.311877422332764   0.001989878863096
%    0.003338814973831   0.867361602783203   0.316131801605225   0.004332161247730
%    0.004064394235611   1.039696578979492   0.308820858001709   0.003435157537460
%    0.007371908426285   1.476032257080078   0.533445129394531   0.000345843285322
%    0.001223269775510   0.727508468627930   0.147539672851562   0.006088330149651
%    0.004511981308460   1.096388778686523   0.329707336425781   0.003085968494415
%    0.002795204222202   0.847239761352539   0.186861610412598   0.004641666114330
%    0.005936530828476   1.214724349975586   0.256533069610596   0.001483623832464
%    0.002679903209209   0.916892700195313   0.280026016235352   0.004851260483265
%    0.000753264874220   0.543055610656738   0.138795099258423   0.006875401735306
%    0.006581504344940   1.175639648437500   0.286277542114258   0.001095846295357
%    0.001464165002108   0.629545860290527   0.234728240966797   0.005809217691422
%    0.001376052945852   0.661222991943359   0.166008319854736   0.005790936350822
%    0.000693211257458   0.543311309814453   0.118732271194458   0.006560686826706
%    0.003898775577545   1.071934432983398   0.269695091247559   0.003291625976562
%    0.005992536544800   1.267260971069336   0.460817413330078   0.001936322301626
%    0.000630677640438   0.484582519531250   0.099317150115967   0.006687180995941
%    0.000668008029461   0.606775398254395   0.084868202209473   0.006785557866096
%    0.000913734063506   0.631240806579590   0.093477802276611   0.006559562683105
%    0.003420960605145   0.976324234008789   0.281819038391113   0.004085080623627
%    0.004019453525543   1.218638763427734   0.420235633850098   0.003479030728340
%    0.001056955233216   0.614048118591309   0.127768344879150   0.005967635512352
%    0.003192428648472   1.012978134155273   0.217884540557861   0.004343386590481
%    0.004981375634670   1.140769042968750   0.241959037780762   0.002657686471939
%    0.001215327829123   0.657607727050781   0.118116569519043   0.005550574064255
%    0.003601420521736   0.968219299316406   0.349507980346680   0.004588645398617
%    0.003972346186638   1.107997817993164   0.424730644226074   0.003688510060310
%    0.001284062117338   0.710308303833008   0.103462629318237   0.005851390361786
%    0.001448235362768   0.647645568847656   0.309835948944092   0.006113651990891
%    0.004963103234768   1.107601394653320   0.520706977844238   0.002139875888824
%    0.002422955036163   0.769995956420898   0.210508308410645   0.003974082171917
%    0.002372161149979   0.779569625854492   0.161610946655273   0.004317984580994
%    0.004921101629734   1.211427001953125   0.388109207153320   0.002074864357710];

% ----------------------(eliminate multiple peaks)--------------------

% criterion1: K must be POSITIVE (K>0)! Eliminate all pksNew with
% corresponding K values that are <=0 OR K == NaN (note: K == NaN means
% that nlinfit was not able to solve for parameters K, D, w, offset)
% 
% The brightest regions of a particle are towards it center...
% 
% K < 0 means that the sphere spread function(ssf) for the particle is 
% inverted so that the center of a particle is the darkest regions of the
% particle.  This is NOT possible in Inorm!
% 
% K = 0 means that the ssf of a particle is flat (i.e., the spherical
% profile is a horizontal line). This also is impossible! 
 
pksToEliminate = find(K(:)<=0 | isnan(K)); %there aren't any



% criterion2: Particle size is known (i.e., 5um +/-0.44um silica,
% manufactured by Bang Laboratories). All particles identified by 'pksNew'
% should have sizes ~Dactual +/- 10% of Dactual 
% 
% IMPORTANT NOTE: Dactual or D are really just parameters of ssf() NOT
% particle size.  You found empirically that particle size is usually some
% combination of D and w
% 
% But some particles are slightly less than 10% larger or smaller than
% Dactual. SO YOU MUST PLOT THE PARTICLE SIZES AS DONE BELOW, MANUALLY 
% LOOK AT THESE PARTICLES, AND ACCORDINGLY ADJUST lBound and uBound TO
% EITHER KEEP OR DELETE THEM. THIS MANUAL ADJUSTMENT MUST BE DONE WITH EACH
% PASS OF YOUR PARTICLE TRACKING ROUTINE!

lBound = 0.7*Dactual;
uBound = 1.3*Dactual;
lowerBound = lBound*ones([size(D,1),1],'single');
upperBound = uBound*ones([size(D,1),1],'single');

figure(1)
hold off
plot(D,'b-*')
hold on
plot(lowerBound,'r--')
plot(upperBound,'r--')
title(sprintf('D = 5.06 +/- 0.44 um according to Bang Labs',Dactual),'fontsize',20)
xlabel(sprintf('actually plotted:\n (%.03f*Dactual), r-- upper \n(%.03f*Dactual), r-- lower \n Dactual = %.02f pixels (size of ideal maskExp)',(uBound/Dactual),(lBound/Dactual),Dactual),'fontsize',20)

pksToEliminate = find(D>uBound | D<lBound); %49 particles

for np = 1:length(pksToEliminate)
    [v,center] = getCutoutSectVertices(pksNew(pksToEliminate(np),:),D(pksToEliminate(np)),size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    figure(2)
    hold off
    h(1) = plot(rNew,avgProfile(startIdx:stopIdx),'b*');
    hold on
    h(2) = plot(rNew,ssf([K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np))],rNew),'k-');
    legend(h,'InormRed','mask')
    title(sprintf('pksNew(%d of %d,:) = [%.04f, %.04f, %.04f] \n K = %.04f, D = %.04f, w = %.04f, offset = %.04f \n bad D because D < %.04f OR D > %.04f',pksToEliminate(np),size(pksNew,1),pksNew(pksToEliminate(np),1),pksNew(pksToEliminate(np),2),pksNew(pksToEliminate(np),3),K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),lowerBound(1),upperBound(1)),'fontsize',20)
    xlabel(sprintf('Rsquared = %.04f',Rsquared(pksToEliminate(np))),'fontsize',20)

    [voronoiVol, Rmask] = peakPlacementB(center,size(cutoutSect),D(pksToEliminate(np)),w(pksToEliminate(np)));
    mask = getCalcImg(Rmask,K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),voronoiVol);
    
    sphericalRegion = find(voronoiVol == 1);
    diff = cutoutSect;
    diff(sphericalRegion) = abs(cutoutSect(sphericalRegion)-mask(sphericalRegion));       
    
    b = 1:numel(mask); 
    b(sphericalRegion) = []; %indices to identify region of the image where there isn't a particle
    mask(b) = cutoutSect(b);
    
    figure(3)
    hold off
    simage([cutoutSect(:,:,round(center(3))), max(cutoutSect(:))*ones([size(cutoutSect,1) 1],'single'), mask(:,:,round(center(3))), max(cutoutSect(:))*ones([size(cutoutSect,1) 1],'single'), diff(:,:,round(center(3)))]);
    colorbar
    hold on
    plot(center(2),center(1),'k*','markersize',20)
    title(sprintf('InormRed \t mask \t abs(mask - InormRed)'),'fontsize',20)
    xlabel(sprintf('pksNew(%d of %d,:) = [%.04f, %.04f, %.04f] \n K = %.04f, D = %.04f, w = %.04f, offset = %.04f \n bad D because D < %.04f OR D > %.04f',pksToEliminate(np),size(pksNew,1),pksNew(pksToEliminate(np),1),pksNew(pksToEliminate(np),2),pksNew(pksToEliminate(np),3),K(pksToEliminate(np)),D(pksToEliminate(np)),w(pksToEliminate(np)),offset(pksToEliminate(np)),lowerBound(1),upperBound(1)),'fontsize',20)
    pause
end


% criterion3: Now identfity multiple peaks locating the same particle. 
% (i.e., peaks or centers of maskExp's separated by less than Dactaul) 
% This can happen because of Monte Carlo minimization of Rsquared. 
% The random movement of maskExp (i.e., peak location) may move a peak on 
% the edge of one particle in Inorm onto a neighboring particle. 
% Keep the particle with the lowest value of Rsquared. (note: you cannot
% just use getNeighborsValid.m for this because this doesn't compare
% Rsquared values of each peak it just keeps the first peak and eliminates
% all subsequent peaks/mask images identifying the same particle)
[dist, pksOrder] = orderPeaks(pksNew);

dummyVal = 0; %this should be orders of magnitude larger than any real Rsquared value
RsquaredCompare = dummyVal*ones([size(pksNew,1), size(pksNew,1)],'single');

% How to read RsquaredCompare...
% 
% provided RsquaredCompare ~= dummyVal
% 
% RsquaredCompare(1,np) = value of Rsquared for particle 1 NOT np
% RsquaredCompare(2,np) = value of Rsquared for particle 2 NOT np
% ...
% RsquaredCompare(np,np) = value of Rsquared for particle np
tic;
for np = 1:size(pksNew,1)
    idxClose = find((dist(:,np) < Dactual*0.4)); %positions in pksNew that are too close to pksNew(np,:) 
    disp(sprintf('-------------\nfinding peaks too close to pksNew(%d,:)...',np))
    for npClose = 1:length(idxClose)%goes through positions too close to pksNew(np,:) and computes Rsquared
        disp(sprintf('pksNew(%d,:) is too close',idxClose(npClose)))
        [v, peakNew] = getCutoutSectVertices(pksNew(idxClose(npClose),:),2*D(idxClose(npClose)),size(InormRed)); %vertices of rectangular portion of 'Inorm' that is (2*Dactual)^3 and centered on validPks(np,:) 
        [naught,RsquaredCompare(idxClose(npClose),np)] = removeParticle(peakNew,K(idxClose(npClose)),D(idxClose(npClose)),w(idxClose(npClose)),offset(idxClose(npClose)),InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3)),idxClose(npClose));
    end
    disp(sprintf('-------------\n'))
end 
elapsedTime = toc; %6.351172820000000e+02 = 10.709min for 1054 particles

RsquaredCompare1 = RsquaredCompare;
for np = 1:size(pksNew,1)
    idxClose = find((dist(:,np) < Dactual*0.4)); %positions in pksNew that are too close to pksNew(np,:) 
    if(numel(idxClose)>1)%>1 because dist(np,np) == 0 so idxClose will always be at least 1 
        disp(sprintf('particle %d has %d peaks: ',np,numel(idxClose)))
        disp(idxClose)
        minRsquared = min(RsquaredCompare(idxClose,np));
        numMinRsquared = find(RsquaredCompare(idxClose,np) == minRsquared);
        if (length(numMinRsquared) > 1)%in the event that 2 or more peaks are EXACTLY the same or have exactly the same minimum Rsquared value then just select one peak delete all others
            idxClose(numMinRsquared(1)) = [];%one of the multipe peaks having the same min. value of Rsquared is deleted from idxClose
            npDelete = 1:numel(idxClose);%ALL other peak in idxClose are tagged for deletion
        end
        if (length(numMinRsquared) == 1) %only one peak has the minimum Rsquared value
            npDelete = find(RsquaredCompare(idxClose,np) > min(RsquaredCompare(idxClose,np)));
        end
        RsquaredCompare(idxClose(npDelete),np) = NaN;
    end
end


pksToEliminate = find(isnan(diag(RsquaredCompare)));

pksToEliminate(pksToEliminate >= npDimer(1)) = []; %There aren't any!

% % pksToEliminate2 is a list that provides the shortest method of viewing
% % all peaks scheduled for deletion
% pksToEliminate2 = pksToEliminate;
% 
% for jj = 1:numel(pksToEliminate2)
%     np = pksToEliminate2(jj);
%     
%     if(~isnan(np))
%         idxClose = find((dist(:,np) < Dactual*0.4)); %positions in pksNew that are too close to pksNew(np,:) 
%         idxClose(idxClose == np) = [];
% 
%         for ii = 1:numel(idxClose)
%             pksToEliminate2(pksToEliminate2 == idxClose(ii)) = nan;
%         end
%     end
% end
% 
% temp = isnan(pksToEliminate2);
% tempDelete = find(temp == 1);
% pksToEliminate2(tempDelete) = [];


% --(start: move peaks to minimize error and locate at pixel resolution)--

clearvars -except K D w offset pksNew Inorm* Kactual Dactual wActual offsetActual 

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

numTasks = 1; %one particle == one task
tasksPerJob = 1:numTasks:size(pksNew,1);
for j = 1:length(tasksPerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue
    disp(sprintf('\n------------job%02d, created with the following tasks------------',j));

    startTask = tasksPerJob(j);
    stopTask = startTask+numTasks-1; %subtract to account for MATLAB indexing displacement
   
    if(stopTask > size(pksNew,1)) 
        stopTask = size(pksNew,1);
    end
    
    for np = startTask:stopTask %slice == task #
        [v,center] = getCutoutSectVertices(pksNew(np,:),0.4*D(np),size(InormRed)); %vertices of rectangular portion of 'InormRed' that is (2*Dactual)^3 and centered on pksNew(np,:) 
        createTask(job(j),@getPosition,2,{v(1,:),center,K(np),D(np),w(np),offset(np),InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3))});
        disp(sprintf('locating pksNew(%d,:) at pixel resolution...',np));
    end
end

% submit one job at a time to be computed by all available slave nodes
for j = 1:length(job)
    submit(job(j));
    disp(sprintf('job%02d, submitted',j))
end 
% Actual time for 0.6*D(np):
%   first job started at: Mon Aug 03 19:36:33 EDT 2015
%   last job finished at: Mon Aug 03 22:08:14 EDT 2015
%   total time: 2hr31min41sec for locating 420 particles


% Actual time for 0.4*D(np):
%   first job started at: Tue Aug 04 15:11:51 EDT 2015
%   last job finished at: Tue Aug 04 16:21:28 EDT 2015
%   total time: 1hr09min37sec for locating 421 particles

pksNewer = zeros(size(pksNew),'single');
Rsquared = zeros([8,size(pksNew,1)],'single'); 
for j = 1:length(job)
    state = job(j).State;
    
    startTask = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        try 
            results = fetchOutputs(job(j));
            disp(sprintf('\n------------job%02d, COMPLETED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksNewer(startTask+np-1,:) = results{np,1}; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = results{np,2};
                disp(sprintf('pksNew(%d,:) becomes pksNewer(%d,:)',startTask+np-1,startTask+np-1));
            end
        catch err 
            disp(sprintf('\n------------job%02d, FAILED with the following tasks------------',j));
            for np = 1:numel(job(j).Tasks)%fetchOutputs(job(j)) doesn't work so you have to read the number of tasks using numel(job(j).Tasks)
                pksNewer(startTask+np-1,:) = NaN; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = NaN;
                disp(sprintf('pksNewer(%d,:) = NaN',startTask+np-1));
            end
        end%end of try
    end
end

temp = isnan(pksNewer(:,1));
npRedo = find(temp(:) == 1); %There aren't any!

% pksNewer =[
% 
%          597         417         284
%          719         375         358
%           10         801         553
%          403         765         536
%          660         263         622
%         1017         967         656
%          885         277         687
%          446         406         764
%          411         517         703
%         1018         162         699
%          744         582         719
%          185         276         171
%          251        1007         174
%          746         861         146
%           84         411         173
%          497         970         152
%          396         281         174
%          906         403         177
%          270          61         146
%          395          81         145
%           91         121         146
%          225         371         141
%          236         476         148
%          108         501         145
%          153         599         146
%          338         497         147
%          207         742         143
%          492         829         144
%          697         754         143
%          821         741         143
%          920         693         150
%          945         592         150
%          963         497         148
%          877         315         145
%          794         215         142
%          698         236         149
%          627         331         145
%          493         251         149
%          843          43         153
%           57         296         211
%          260         102         230
%          284         196         190
%          288         278         239
%          180         426         212
%          280         564         184
%          436         503         182
%          546         557         188
%          606         412         194
%          605         140         201
%          885         235         196
%          989         416         223
%          811         893         208
%          726         939         199
%          301         925         188
%          522         151         255
%          536         275         276
%          637         275         255
%          815         146         273
%          763         361         274
%          698         428         288
%          662         593         233
%          469         558         255
%          923         648         237
%          623        1017         265
%          374         938         247
%          296         995         257
%          118         778         287
%           69         646         285
%           22         262         296
%          274         201         289
%          335         269         314
%         1024         517         285
%          676         749         279
%          530         900         305
%           47         480         312
%          153         754         359
%          545         658         365
%          624         821         355
%          487         778         361
%          255         383         391
%          728         732         384
%          887         652         412
%          826         731         406
%          789         824         403
%          986         746         409
%          876         169         409
%          153         404         398
%          198         628         392
%           86         784         436
%          527         860         393
%          706         450         411
%           12         340         443
%          936         507         433
%          822         883         460
%          185         768         441
%          947         238         519
%          532         953         532
%          413         868         526
%           36         342         541
%          207         241         549
%          516          22         558
%          714         135         564
%          529         373         564
%          213         503         550
%          602         755         539
%          186        1015         534
%          521         862         583
%         1017         719         578
%          194         607         582
%          930         423         618
%          275         455         659
%          158         649         666
%          558           9         688
%          986         441         697
%          883         875         697
%          245         162         691
%          505         192         675
%          690         334         670
%          369         569         769
%          466         647         763
%          743         411         730
%          530         108         787
%          164         133         797
%          214         661         818
%          594         157         826
%          846         516         147
%          812         449         149
%          766         660         163
%          825         675         198];


% check that pksNewer is correct. Before accepting pksBest
cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single2ndPass_id_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
fclose(fid);

order = temp(1,:)';
pksNew(order,1) = temp(2,:)';
pksNew(order,2) = temp(3,:)';
pksNew(order,3) = temp(4,:)';

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers2ndPass_row_col_slice.txt','r');
dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksDimers(:,1) = dimers(1,:)';
pksDimers(:,2) = dimers(2,:)';
pksDimers(:,3) = dimers(3,:)';

npDimer = (size(pksNew,1)+1):(size(pksDimers,1)+size(pksNew,1));
pksNew = vertcat(pksNew,pksDimers);

npGood = 1:size(pksNew,1);
npDelete = find(isnan(pksNew(:,1)) == 1);

npGood(npDelete) = [];

pksNew(npDelete,:) = [];

 
for np = 121:size(pksNew,1)
    startSlice = nanmin([pksNewer(np,3), pksNew(np,3)])-15;
    stopSlice = nanmax([pksNewer(np,3), pksNew(np,3)])+15;

    if(startSlice < 1)
        startSlice = 1;
    end

    if(stopSlice > size(InormRed,3))
        stopSlice = size(InormRed,3);
    end
    
    for slice = startSlice:stopSlice
        figure(1)
        subplot(1,2,1)
        hold off
        simage(InormRed(:,:,slice));
        hold on
        plot(pksNewer(np,2),pksNewer(np,1),'k*','markersize',20)
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNewer(%d,:) = [%d, %d, %d] (k*)\npksNew(%d of %d,:) = [%d, %d, %d](g*)',np,pksNewer(np,1),pksNewer(np,2),pksNewer(np,3),np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)

        subplot(1,2,2)
        hold off
        simage(InormGreen(:,:,slice));
        hold on
        plot(pksNewer(np,2),pksNewer(np,1),'k*','markersize',20)
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNew(%d of %d,:) = [%d, %d, %d](k*)',np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)
        pause
    end
end

% -----(start: locate position of particle at sub-pixel resolution)-------

clearvars -except K D w offset pksNewer Inorm* Kactual Dactual wActual offsetActual 

clear pksNewer

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single_pksNewer2ndPass_id_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
fclose(fid);

order = temp(1,:)';
pksNewer(order,1) = temp(2,:)';
pksNewer(order,2) = temp(3,:)';
pksNewer(order,3) = temp(4,:)';

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers2ndPass_row_col_slice.txt','r');
dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
fclose(fid);

pksDimers(:,1) = dimers(1,:)';
pksDimers(:,2) = dimers(2,:)';
pksDimers(:,3) = dimers(3,:)';

npDimer = (size(pksNewer,1)+1):(size(pksDimers,1)+size(pksNewer,1));
pksNewer = vertcat(pksNewer,pksDimers);

% for np = 1:size(pksNew,1)
%     startSlice = nanmin([pksNewer(np,3), pksNew(np,3)])-7;
%     stopSlice = nanmax([pksNewer(np,3), pksNew(np,3)])+7;
% 
%     if(startSlice < 1)
%         startSlice = 1;
%     end
% 
%     if(stopSlice > size(InormRed,3))
%         stopSlice = size(InormRed,3);
%     end
%     
%     for slice = startSlice:stopSlice
%         figure(1)
%         subplot(1,2,1)
%         hold off
%         simage(InormRed(:,:,slice));
%         hold on
%         plot(pksNewer(np,2),pksNewer(np,1),'k*','markersize',20)
%         plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
%         title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
%         xlabel(sprintf('pksNewer(%d,:) = [%d, %d, %d] (k*)\npksNew(%d of %d,:) = [%d, %d, %d](g*)',np,pksNewer(np,1),pksNewer(np,2),pksNewer(np,3),np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)
% 
%         subplot(1,2,2)
%         hold off
%         simage(InormGreen(:,:,slice));
%         hold on
%         plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
%         plot(pksNewer(np,2),pksNewer(np,1),'k*','markersize',20)
%         title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
%         xlabel(sprintf('pksNew(%d of %d,:) = [%d, %d, %d](k*)',np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)
%         pause
%     end
% end

numTry = 10; %# of tries to calculate param(np,:)

param = zeros([size(pksNew,1), 4]);
Rsquared = zeros([size(pksNew,1), 1],'single');

tic;
for np = 1:size(pksNew,1)
    [v,center] = getCutoutSectVertices(pksNew(np,:),0.75*Dactual,size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    %bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    %rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    [param(np,:), Rsquared(np)] = getParameters(avgProfile(startIdx:stopIdx),Dactual,wActual,numTry);
    
    if(~isnan(param(np,1)))
        disp(sprintf('SUCCESSFUL computation of parameters for pksNew(%d,:)',np));
    end
    
    if(isnan(param(np,1)))
        disp(sprintf('FAILED computation of parameters for pksNew(%d,:)',np));
    end
end
elapsedTime = toc; %16.5133669seconds for 285 particles

K = param(:,1);
D = param(:,2);
w = param(:,3);
offset = param(:,4);

% param =(1.0e+02)*[
% 
%    0.001626689136028   0.700081710815430   0.130961837768555   0.005748860836029
%    0.004213369488716   0.996036376953125   0.217023582458496   0.003216860592365
%    0.001373026967049   0.702410354614258   0.180403270721436   0.005888956785202
%    0.001379778385162   0.678561096191406   0.156764383316040   0.005711719989777
%    0.001069074049592   0.511003913879394   0.125844030380249   0.005694840550423
%    0.001155533343554   0.616300582885742   0.228901710510254   0.006325969696045
%    0.001620403528214   0.779639053344727   0.120543031692505   0.005939540863037
%    0.003275525569916   0.980249938964844   0.239234390258789   0.004348288476467
%    0.001366440504789   0.651432800292969   0.158616085052490   0.006138064861298
%    0.000561881996691   0.486667098999023   0.099430541992187   0.006916378140450
%    0.001095413714647   0.672352523803711   0.099170017242432   0.005872829556465
%    0.002364481091499   0.772096405029297   0.178349151611328   0.004530478417873
%    0.002185411900282   0.742180023193359   0.172064189910889   0.004815630614758
%    0.006511397957802   1.140401916503906   0.373016357421875   0.000811772122979
%    0.002093155980110   0.750231933593750   0.134220943450928   0.004499850869179
%    0.002113492786884   0.729236373901367   0.173541584014893   0.004519376754761
%    0.001694871932268   0.690377426147461   0.129800052642822   0.005144584178925
%    0.005744611024857   1.156977844238281   0.297900009155273   0.001173442378640
%    0.001809971630573   0.690082168579102   0.171362323760986   0.004847627878189
%    0.002244126200676   0.702211380004883   0.185413589477539   0.004900319278240
%    0.001877402514219   0.667368240356445   0.147517976760864   0.005075790286064
%    0.003459694087505   0.883341674804688   0.287954540252686   0.003308722376823
%    0.002882576584816   0.818693161010742   0.236235466003418   0.004363935887814
%    0.001924621164799   0.622936363220215   0.191988048553467   0.005083823204041
%    0.006041778326035   1.045429458618164   0.181191635131836   0.000934812799096
%    0.004718746244907   0.971140136718750   0.382003440856934   0.002871006727219
%    0.001919619888067   0.647881927490234   0.133355522155762   0.004808759689331
%    0.001624137759209   0.626550521850586   0.141737451553345   0.005088241696358
%    0.005026643872261   1.065180892944336   0.366894989013672   0.001977983117104
%    0.001409883052111   0.476794281005859   0.126459751129150   0.004962096810341
%    0.005860819220543   1.156775207519531   0.362996101379395   0.001343727409840
%    0.001595026105642   0.650788879394531   0.137286481857300   0.005274205207825
%    0.001241949126124   0.590686378479004   0.127111473083496   0.005484029650688
%    0.005335355401039   1.089024658203125   0.274056034088135   0.001630156934261
%    0.004791601598263   1.102586212158203   0.222138309478760   0.001937154829502
%    0.001446709632874   0.643937988281250   0.133395891189575   0.005426278114319
%    0.001420859247446   0.634015998840332   0.140148706436157   0.005381944179535
%    0.002054664045572   0.766882019042969   0.176652355194092   0.004691389203072
%    0.001377572417259   0.573375816345215   0.124520902633667   0.005742106437683
%    0.002445891201496   0.768236541748047   0.172218017578125   0.004643688797951
%    0.002547087073326   0.824421005249023   0.266195907592773   0.004559686779976
%    0.001326843053102   0.654190368652344   0.120268554687500   0.005522387027740
%    0.002065649330616   0.658573455810547   0.142349386215210   0.005135874152184
%    0.004562323391438   1.036813507080078   0.483943214416504   0.003123162090778
%    0.001974675506353   0.735239791870117   0.204902553558350   0.004951849281788
%    0.004871541857719   0.976655426025391   0.252081317901611   0.002327967584133
%    0.001192963272333   0.626847915649414   0.067487297058105   0.005549242496490
%    0.002618892192841   0.763054733276367   0.182832736968994   0.004381961822510
%    0.004168935716152   1.002888259887695   0.225747623443604   0.002653227746487
%    0.002429427951574   0.798900375366211   0.199062385559082   0.004336560964584
%    0.001223559826612   0.601633415222168   0.156105251312256   0.005613782405853
%    0.001580297499895   0.703565216064453   0.121904039382935   0.005642821192741
%    0.001984894275665   0.702149887084961   0.147528181076050   0.004810744822025
%    0.003420042991638   0.876248855590820   0.219333629608154   0.003530635535717
%    0.002870749235153   0.911317291259766   0.228953056335449   0.004004518687725
%    0.001529164761305   0.587966423034668   0.132678098678589   0.005622845888138
%    0.001601684093475   0.701760482788086   0.164174976348877   0.005391453504562
%    0.002112801074982   0.753033523559570   0.190997867584229   0.004844561815262
%    0.003355203568935   0.942418823242188   0.233060150146484   0.003795498311520
%    0.005419316887856   1.118215484619141   0.295230693817139   0.001914349347353
%    0.002620644867420   0.836034317016602   0.165793685913086   0.004093385636806
%    0.003536171913147   0.920922851562500   0.217005519866943   0.003646059334278
%    0.001401395350695   0.536404838562012   0.116088199615479   0.005439224243164
%    0.000939638763666   0.608305740356445   0.091707134246826   0.005720565319061
%    0.001361616849899   0.589052581787109   0.113090114593506   0.005774816274643
%    0.001193960309029   0.527999191284180   0.104667787551880   0.006354433298111
%    0.005626214146614   1.140025024414062   0.280113658905029   0.001741332858801
%    0.002205626368523   0.735335998535156   0.177536983489990   0.004850660264492
%    0.001757955402136   0.599343109130859   0.203318920135498   0.005675869584084
%    0.004269691407681   0.943299865722656   0.346356620788574   0.003406108617783
%    0.003735060691833   1.008851547241211   0.299804859161377   0.003649284839630
%    0.000773053839803   0.312806930541992   0.134984989166260   0.006900937557220
%    0.002478671371937   0.822429580688477   0.179206504821777   0.004524542391300
%    0.001363940984011   0.616501083374023   0.131515703201294   0.006111965775490
%    0.006402861475945   1.243406219482422   0.411938438415527   0.001291826218367
%    0.002724638283253   0.824125976562500   0.256528434753418   0.004684564769268
%    0.001476149111986   0.602739830017090   0.127141647338867   0.005502166152000
%    0.002739206850529   0.729943542480469   0.218589859008789   0.003996579349041
%    0.003141257762909   0.919690322875977   0.243488998413086   0.003758335411549
%    0.000949956849217   0.569686889648438   0.064372725486755   0.005861158967018
%    0.002501482665539   0.787898406982422   0.211962528228760   0.004817315042019
%    0.006941682100296   1.393177032470703   0.464124870300293   0.000802412182093
%    0.001181585043669   0.620869178771973   0.106589927673340   0.005617925524712
%    0.001641951948404   0.619386978149414   0.153954563140869   0.005740472078323
%    0.003535416424274   1.052474670410156   0.344355392456055   0.003868314027786
%    0.004412575960159   0.945795974731445   0.194506320953369   0.002775108814240
%    0.001303784400225   0.627402801513672   0.115606327056885   0.005934365391731
%    0.001428082883358   0.625668449401856   0.082013301849365   0.005278236865997
%    0.003689153194427   0.976111221313477   0.298747272491455   0.003990989625454
%    0.001621634960175   0.603891563415527   0.169248027801514   0.005941306352615
%    0.001059360951185   0.584055328369141   0.093173522949219   0.006229863166809
%    0.001203651204705   0.565502395629883   0.155331296920776   0.006060971021652
%    0.006803213953972   1.394926147460938   0.429790878295898   0.001176348924637
%    0.000811292678118   0.597857818603516   0.071349234580994   0.006653156280518
%    0.002663824856281   0.856084365844727   0.208784294128418   0.004867908358574
%    0.003670360147953   1.061364822387695   0.343649406433105   0.003970317840576
%    0.005566983819008   1.232793426513672   0.311877422332764   0.001989878863096
%    0.003338814973831   0.867361602783203   0.316131801605225   0.004332161247730
%    0.004064394235611   1.039696578979492   0.308820858001709   0.003435157537460
%    0.007371908426285   1.476032257080078   0.533445129394531   0.000345843285322
%    0.001223269775510   0.727508468627930   0.147539672851562   0.006088330149651
%    0.004511981308460   1.096388778686523   0.329707336425781   0.003085968494415
%    0.002795204222202   0.847239761352539   0.186861610412598   0.004641666114330
%    0.005936530828476   1.214724349975586   0.256533069610596   0.001483623832464
%    0.002679903209209   0.916892700195313   0.280026016235352   0.004851260483265
%    0.000753264874220   0.543055610656738   0.138795099258423   0.006875401735306
%    0.006581504344940   1.175639648437500   0.286277542114258   0.001095846295357
%    0.001464165002108   0.629545860290527   0.234728240966797   0.005809217691422
%    0.001376052945852   0.661222991943359   0.166008319854736   0.005790936350822
%    0.000693211257458   0.543311309814453   0.118732271194458   0.006560686826706
%    0.003898775577545   1.071934432983398   0.269695091247559   0.003291625976562
%    0.005992536544800   1.267260971069336   0.460817413330078   0.001936322301626
%    0.000630677640438   0.484582519531250   0.099317150115967   0.006687180995941
%    0.000668008029461   0.606775398254395   0.084868202209473   0.006785557866096
%    0.000913734063506   0.631240806579590   0.093477802276611   0.006559562683105
%    0.003420960605145   0.976324234008789   0.281819038391113   0.004085080623627
%    0.004019453525543   1.218638763427734   0.420235633850098   0.003479030728340
%    0.001056955233216   0.614048118591309   0.127768344879150   0.005967635512352
%    0.003192428648472   1.012978134155273   0.217884540557861   0.004343386590481
%    0.004981375634670   1.140769042968750   0.241959037780762   0.002657686471939
%    0.001215327829123   0.657607727050781   0.118116569519043   0.005550574064255
%    0.003601420521736   0.968219299316406   0.349507980346680   0.004588645398617
%    0.003972346186638   1.107997817993164   0.424730644226074   0.003688510060310
%    0.001284062117338   0.710308303833008   0.103462629318237   0.005851390361786
%    0.001448235362768   0.647645568847656   0.309835948944092   0.006113651990891
%    0.004963103234768   1.107601394653320   0.520706977844238   0.002139875888824
%    0.002422955036163   0.769995956420898   0.210508308410645   0.003974082171917
%    0.002372161149979   0.779569625854492   0.161610946655273   0.004317984580994
%    0.004921101629734   1.211427001953125   0.388109207153320   0.002074864357710];

% -----(start: locate position of particle at sub-pixel resolution)-------

clearvars -except K D w offset pksNewer Inorm* Kactual Dactual wActual offsetActual 

% clear pksNewer
% 
% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('single_pksNewer2ndPass_id_row_col_slice.txt','r');
% temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
% fclose(fid);
% 
% order = temp(1,:)';
% pksNewer(order,1) = temp(2,:)';
% pksNewer(order,2) = temp(3,:)';
% pksNewer(order,3) = temp(4,:)';
% 
% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('2mers2ndPass_row_col_slice.txt','r');
% dimers = fscanf(fid,'%f %f %f\n',[3, Inf]);
% fclose(fid);
% 
% pksDimers(:,1) = dimers(1,:)';
% pksDimers(:,2) = dimers(2,:)';
% pksDimers(:,3) = dimers(3,:)';
% 
% npDimer = (size(pksNewer,1)+1):(size(pksDimers,1)+size(pksNewer,1));
% pksNewer = vertcat(pksNewer,pksDimers);

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

rangeColormap = 255;

numTasks = 1; %one particle == one task
tasksPerJob = 1:numTasks:size(pksNewer,1);
for j = 1:length(tasksPerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue
    disp(sprintf('\n------------job%02d, created with the following tasks------------',j));

    startTask = tasksPerJob(j);
    stopTask = startTask+numTasks-1; %subtract to account for MATLAB indexing displacement
   
    if(stopTask > size(pksNewer,1)) 
        stopTask = size(pksNewer,1);
    end
    
    for np = startTask:stopTask %slice == task #
        particleRadius = ((D(np)/2)+(2*w(np))); %this is the apperant particle radius
        [v,center] = getCutoutSectVertices(pksNewer(np,:),particleRadius,size(InormRed)); %vertices of rectangular portion of 'Inorm' that is (2*Dactual)^3 and centered on validPks(np,:) 
        cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
        centerIdx = sub2ind(size(cutoutSect),center(1),center(2),center(3));
        createTask(job(j),@getPositionSubPixel,2,{v(1,:),centerIdx,K(np),D(np),w(np),offset(np),cutoutSect,rangeColormap});
        disp(sprintf('locating pksNewer(%d,:) at SUB-pixel resolution...',np));
    end
end

% submit one job at a time to be computed by all available slave nodes
for j = 1:length(job)
    submit(job(j));
    disp(sprintf('job%02d, submitted',j))
end 

% 1st job to start, started at: Tue Aug 11 23:12:10 EDT 2015
% last job to finish, finished at: Wed Aug 12 04:05:09 EDT 2015
% total time: 4hrs52min59seconds for 129 particles

pksBest = zeros(size(pksNewer),'single');
Rsquared = zeros([11,size(pksNewer,1)],'single'); 
for j = 1:length(job)
    state = job(j).State;
    
    startTask = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        try 
            results = fetchOutputs(job(j));
            disp(sprintf('\n------------job%02d, COMPLETED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksBest(startTask+np-1,:) = results{np,1}; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = results{np,2};
                disp(sprintf('pksNewer(%d,:) becomes pksBest(%d,:)',startTask+np-1,startTask+np-1));
            end
        catch err 
            disp(sprintf('\n------------job%02d, FAILED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksBest(startTask+np-1,:) = NaN; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = NaN;
                disp(sprintf('pksNewer(%d,:) = NaN',startTask+np-1));
            end
        end%end of try
    end
end

temp = isnan(pksBest(:,1));
npRedo = find(temp(:) == 1); %There aren't any!

% pksBest =(1.0e+03)*[
% 
%    0.5960432   0.4110000   0.2861373
%    0.7199843   0.3759922   0.3589765
%    0.0180000   0.7920000   0.5522549
%    0.4040000   0.7645529   0.5350000
%    0.6607647   0.2632823   0.6210039
%    1.0160000   0.9660392   0.6569725
%    0.8840078   0.2760275   0.6860510
%    0.4497608   0.3990157   0.7640039
%    0.4111216   0.5162353   0.7020078
%    1.0170078   0.1610039   0.6980157
%    0.7444078   0.5830000   0.7180000
%    0.1840902   0.2750118   0.1640039
%    0.2519961   1.0060000   0.1749686
%    0.7469059   0.8550157   0.1450078
%    0.0849608   0.4050000   0.1739490
%    0.4950823   0.9680078   0.1510000
%    0.3969294   0.2804784   0.1730000
%    0.9050000   0.4040000   0.1762823
%    0.2729568   0.0530000   0.1400078
%    0.3952628   0.0800000   0.1440000
%    0.0900824   0.1219961   0.1450000
%    0.2320000   0.3611020   0.1400078
%    0.2359961   0.4720118   0.1380235
%    0.1070980   0.5020000   0.1390000
%    0.1549686   0.5920039   0.1435843
%    0.3420000   0.4921059   0.1350000
%    0.2079176   0.7419686   0.1420000
%    0.4979725   0.8180000   0.1430235
%    0.6979961   0.7531255   0.1390039
%    0.8219961   0.7419764   0.1420235
%    0.9190039   0.6934510   0.1420000
%    0.9460039   0.5919725   0.1420000
%    0.9620078   0.4964784   0.1420000
%    0.8760236   0.3141333   0.1440000
%    0.7930000   0.2159922   0.1412549
%    0.6970039   0.2357961   0.1480039
%    0.6260823   0.3300157   0.1440000
%    0.4936274   0.2460157   0.1480078
%    0.8422628   0.0399882   0.1410000
%    0.0580157   0.2920000   0.2112000
%    0.2609686   0.1010000   0.2399961
%    0.2830353   0.1900039   0.1890431
%    0.2889804   0.2770039   0.2459961
%    0.1780157   0.4220039   0.2139804
%    0.2789490   0.5640078   0.1800039
%    0.4350078   0.5039922   0.1830000
%    0.5457529   0.5574902   0.1870000
%    0.6069843   0.4110118   0.1949961
%    0.6040078   0.1390157   0.2019961
%    0.8840000   0.2359961   0.1950745
%    0.9880078   0.4100275   0.2220000
%    0.8120000   0.8920196   0.2070314
%    0.7263765   0.9380118   0.1980000
%    0.3000196   0.9230000   0.1889961
%    0.5210039   0.1500000   0.2547804
%    0.5370000   0.2740196   0.2769216
%    0.6360236   0.2755569   0.2560000
%    0.8180000   0.1450039   0.2733686
%    0.7666470   0.3560235   0.2749922
%    0.6980275   0.4270118   0.2889922
%    0.6620039   0.5860039   0.2339882
%    0.4680078   0.5570157   0.2559843
%    0.9220039   0.6470000   0.2362549
%    0.6299922   1.0190078   0.2641059
%    0.3732235   0.9330000   0.2479098
%    0.2969568   0.9940000   0.2580000
%    0.1170078   0.7770157   0.2879961
%    0.0698980   0.6469059   0.2860000
%    0.0229961   0.2610000   0.2950000
%    0.2749725   0.1970000   0.2920667
%    0.3359843   0.2680118   0.3130039
%    1.0230078   0.5160000   0.2842667
%    0.6769961   0.7480078   0.2800000
%    0.5349961   0.8990000   0.3032706
%    0.0480000   0.4760432   0.3149961
%    0.1530000   0.7502235   0.3579922
%    0.5420000   0.6550157   0.3640275
%    0.6244274   0.8310000   0.3560000
%    0.4910196   0.7740000   0.3619961
%    0.2544157   0.3820157   0.3900039
%    0.7270157   0.7270000   0.3930078
%    0.8824667   0.6450353   0.4080000
%    0.8269843   0.7300118   0.4050432
%    0.7893569   0.8230039   0.4039961
%    0.9850000   0.7451765   0.4099922
%    0.8769804   0.1680039   0.4099961
%    0.1569804   0.4030549   0.3910000
%    0.1972902   0.6270118   0.3910039
%    0.0840000   0.7770039   0.4355725
%    0.5280000   0.8609804   0.3940000
%    0.7100000   0.4509647   0.4119725
%    0.0129804   0.3410000   0.4420157
%    0.9300039   0.5079765   0.4350000
%    0.8210039   0.8820039   0.4592039
%    0.1890078   0.7670157   0.4419961
%    0.9500157   0.2330000   0.5189843
%    0.5338902   0.9470039   0.5310039
%    0.4139804   0.8640039   0.5321765
%    0.0399961   0.3410000   0.5419725
%    0.2060314   0.2419922   0.5492314
%    0.5150353   0.0230000   0.5589922
%    0.7130078   0.1340078   0.5649568
%    0.5299922   0.3739765   0.5649961
%    0.2139961   0.5039961   0.5504431
%    0.6050000   0.7540078   0.5430392
%    0.1851686   1.0180078   0.5419922
%    0.5219843   0.8629686   0.5840000
%    1.0160000   0.7180902   0.5770157
%    0.1930039   0.6060078   0.5810078
%    0.9290236   0.4220000   0.6170078
%    0.2759686   0.4559961   0.6599922
%    0.1589765   0.6499961   0.6669882
%    0.5589843   0.0099922   0.6889764
%    0.9850118   0.4419882   0.6960039
%    0.8839882   0.8740118   0.6960039
%    0.2441882   0.1610000   0.6919686
%    0.5059098   0.1929922   0.6760000
%    0.6950000   0.3344784   0.6709255
%    0.3698667   0.5680157   0.7680078
%    0.4669882   0.6460314   0.7620039
%    0.7437177   0.4100078   0.7290000
%    0.5291451   0.1077137   0.7860000
%    0.1630039   0.1334510   0.7960000
%    0.2189882   0.6550157   0.8176470
%    0.5930314   0.1580000   0.8250118
%    0.8559804   0.5109843   0.1380000
%    0.8159882   0.4439216   0.1410039
%    0.7649725   0.6550000   0.1559373
%    0.8339843   0.6680078   0.1989725];

npDimer = 126:size(pksNewer,1);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single2ndPass12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for np = 1:(npDimer(1)-1)
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest(np,1),pksBest(np,2),pksBest(np,3),K(np),D(np),w(np),offset(np));
end
fclose(fid);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers2ndPass12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for np = npDimer(1):2:npDimer(numel(npDimer))
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest(np,1),pksBest(np,2),pksBest(np,3),K(np),D(np),w(np),offset(np));
    fprintf(fid,'%d %f %f %f %f %f %f %f\n\n',np+1,pksBest(np+1,1),pksBest(np+1,2),pksBest(np+1,3),K(np+1),D(np+1),w(np+1),offset(np+1));
end
fclose(fid);

% ------------(start: check particles found in last pass)---------------

clearvars -except Inorm* 

cd('~/poincareProgs/')
run results12may2015_6may2015_3_1920A1000r_13xS3

[dist, pksOrder] = orderPeaks(pksBest);

% view multimers:
for np = 514:size(pksBest,1)
    npp = pksOrder(1:5,np);

    [naught, idx] = sort(pksBest(npp,3));
    npp = npp(idx);
    
%     %for 2mers:
%     npp = union(pksOrder(1:3,np),[np:(np+1)]);
    
%     npp = pksOrder(1:3,np);
    
%     %for 5mers:
%     npp = union(pksOrder(1:6,np),[np:(np+4)]);

%     %for 6mers
%     npp = union(pksOrder(1:6,np),[np:(np+5)]);
    
    clear label 
    disp(sprintf('---------------------------x-----------------------------'));
    for jj = 1:numel(npp)
        label{jj} = num2str(sprintf('%d',npp(jj)));
        disp(sprintf('\tpksBest(%d,:) = [%.03f, %.03f, %.03f], distance(%d,%d) = %.04f',npp(jj),pksBest(npp(jj),1),pksBest(npp(jj),2),pksBest(npp(jj),3),npp(jj),np,dist(npp(jj),np)));
    end
    disp(sprintf('---------------------------x-----------------------------\n'));
   
    
    startSlice = round(min(pksBest(npp,3))-10);
    stopSlice = round(max(pksBest(npp,3))+10);

    if(startSlice < 1)
        startSlice = 1;
    end

    if(stopSlice > size(InormRed,3))
        stopSlice = size(InormRed,3);
    end
    
    for slice = startSlice:stopSlice
        figure(1)
        subplot(1,2,1)
        hold off
        simage(InormRed(:,:,slice));
        hold on
        plot(pksBest(npp,2),pksBest(npp,1),'k*','markersize',20)
        plot(pksBest(np,2),pksBest(np,1),'g*','markersize',20)
        text(pksBest(npp,2),pksBest(npp,1),label,'fontsize',12)
        title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksBest(%d of %d,:) = [%.04f, %.04f, %.04f] \n %d particles plotted',np,size(pksBest,1),pksBest(np,1),pksBest(np,2),pksBest(np,3),numel(npp)),'fontsize',20)

        subplot(1,2,2)
        hold off
        simage(InormGreen(:,:,slice));
        hold on
        plot(pksBest(npp,2),pksBest(npp,1),'k*','markersize',20)
        plot(pksBest(np,2),pksBest(np,1),'g*','markersize',20)
        text(pksBest(npp,2),pksBest(npp,1),label,'fontsize',12)
        title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksBest(%d,3) = %.03f \npksBest(%d,3) = %.03f \npksBest(%d,3) = %.03f',npp(1),pksBest(npp(1),3),npp(2),pksBest(npp(2),3),npp(3),pksBest(npp(3),3)),'fontsize',20)

        pause
    end
end

npDelete = [2,7];

pksBest2(npDelete,:) = [];
K2(npDelete) = [];
D2(npDelete) = [];
w2(npDelete) = [];
offset2(npDelete) = [];

% correct lists of particles found in 2ndPass:
npDimer = 124:size(pksBest2,1);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single2ndPass12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for np = 1:(npDimer(1)-1)
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest2(np,1),pksBest2(np,2),pksBest2(np,3),K2(np),D2(np),w2(np),offset2(np));
end
fclose(fid);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('2mers2ndPass12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for np = npDimer(1):2:npDimer(numel(npDimer))
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest2(np,1),pksBest2(np,2),pksBest2(np,3),K2(np),D2(np),w2(np),offset2(np));
    fprintf(fid,'%d %f %f %f %f %f %f %f\n\n',np+1,pksBest2(np+1,1),pksBest2(np+1,2),pksBest2(np+1,3),K2(np+1),D2(np+1),w2(np+1),offset2(np+1));
end
fclose(fid);

% --------------------(end: particle trakcing 2nd Pass)-----------------


% --------------------(start: particle trakcing 3rd Pass)-----------------
% locate all remaining particles by hand

clearvars -except Inorm* 

cd('~/poincareProgs/')
run results12may2015_6may2015_3_1920A1000r_13xS3

pksBest = double(pksBest);
for slice = 100:size(InormRed,3)
    npTemp1 = find(round(pksBest(:,3)) >= (slice-40));
    npTemp2 = find(round(pksBest(:,3)) <= (slice + 40));
    npp = intersect(npTemp1,npTemp2);
    
    clear label 
    disp(sprintf('---------------------------x-----------------------------'));
    for jj = 1:numel(npp)
        label{jj} = num2str(sprintf('%d',npp(jj)));
        disp(sprintf('\tpksBest(%d,:) = [%d, %d, %d]',npp(jj),pksBest(npp(jj),1),pksBest(npp(jj),2),pksBest(npp(jj),3)));
    end
    disp(sprintf('---------------------------x-----------------------------\n'));  
    
    figure(1)
    subplot(1,2,1)
    hold off
    simage(InormRed(:,:,slice));
    hold on
    plot(pksBest(npp,2),pksBest(npp,1),'k*','markersize',20)
    text(pksBest(npp,2),pksBest(npp,1),label,'fontsize',12)
    title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20);

    subplot(1,2,2)
    hold off
    simage(InormGreen(:,:,slice));
    hold on
    plot(pksBest(npp,2),pksBest(npp,1),'k*','markersize',20)
    text(pksBest(npp,2),pksBest(npp,1),label,'fontsize',12)
    title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)

    pause
end

% -----------------------------x-----------------------------

clearvars -except Inorm* Kactual Dactual wActual offsetActual 

% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('pksNew3rdPass_row_col_slice.txt','r');
% temp = fscanf(fid,'%f %f %f\n',[3,Inf]);
% fclose(fid);
% 
% pksNew(:,1) = temp(1,:)';
% pksNew(:,2) = temp(2,:)';
% pksNew(:,3) = temp(3,:)';

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('pksNew4thPass_row_col_slice.txt','r');
temp = fscanf(fid,'%f %f %f\n',[3,Inf]);
fclose(fid);

pksNew(:,1) = temp(1,:)';
pksNew(:,2) = temp(2,:)';
pksNew(:,3) = temp(3,:)';

for np = 1:size(pksNew,1)
    startSlice = nanmin(pksNew(np,3))-40;
    stopSlice = nanmax(pksNew(np,3))+40;

    if(startSlice < 1)
        startSlice = 1;
    end

    if(stopSlice > size(InormRed,3))
        stopSlice = size(InormRed,3);
    end
    
    for slice = startSlice:stopSlice
        figure(1)
        subplot(1,2,1)
        hold off
        simage(InormRed(:,:,slice));
        hold on
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)

        subplot(1,2,2)
        hold off
        simage(InormGreen(:,:,slice));
        hold on
        plot(pksNew(np,2),pksNew(np,1),'k*','markersize',20)
        title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNew(%d of %d,:) = [%d, %d, %d](k*)',np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)
        pause
    end
end


% --(start: move peaks to minimize error and locate at pixel resolution)--

numTry = 10; %# of tries to calculate param(np,:)

param = zeros([size(pksNew,1), 4]);
Rsquared = zeros([size(pksNew,1), 1],'single');

tic;
for np = 1:size(pksNew,1)
    [v,center] = getCutoutSectVertices(pksNew(np,:),0.75*Dactual,size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    %bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    %rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    [param(np,:), Rsquared(np)] = getParameters(avgProfile(startIdx:stopIdx),Dactual,wActual,numTry);
    
    if(~isnan(param(np,1)))
        disp(sprintf('SUCCESSFUL computation of parameters for pksNew(%d,:)',np));
    end
    
    if(isnan(param(np,1)))
        disp(sprintf('FAILED computation of parameters for pksNew(%d,:)',np));
    end
end
elapsedTime = toc; %16.5133669seconds for 285 particles

K = param(:,1);
D = param(:,2);
w = param(:,3);
offset = param(:,4);

% param =(1.0e+02)*[
% 
%    0.003227757513523   0.697536849975586   0.282194442749023   0.003994338810444
%    0.001934662163258   0.794641647338867   0.185456256866455   0.005564644336700
%    0.001633525788784   0.607884101867676   0.133369512557983   0.005036738514900
%    0.001557430922985   0.655199966430664   0.178628005981445   0.005827140212059
%    0.001722699850798   0.763781738281250   0.148881130218506   0.005776244997978
%    0.001435083597898   0.675904617309570   0.201933364868164   0.005952192544937
%    0.001242846921086   0.619646263122559   0.111098365783691   0.006195434331894
%    0.001728632301092   0.699696731567383   0.198300876617432   0.005647221207619
%    0.003362057805061   1.013091430664063   0.312162303924561   0.003696921169758
%    0.003529913723469   1.194706878662109   0.343828315734863   0.003838486075401
%    0.002187668383121   0.758569030761719   0.167971515655518   0.005163049101830
%    0.001362594664097   0.577598953247070   0.231816158294678   0.005646393299103
%    0.001782626211643   0.801628112792969   0.186470146179199   0.005436823368073
%    0.000949202924967   0.453152847290039   0.121362752914429   0.005418688654900
%    0.001930826157331   0.814857711791992   0.150160083770752   0.005722764134407
%    0.001991224884987   0.760727844238281   0.160855884552002   0.005531769394875
%    0.001976958364248   0.796563034057617   0.231466464996338   0.005355238914490
%    0.002342429310083   0.745821762084961   0.161337871551514   0.005342081785202
%    0.001929880380630   0.715039215087891   0.167042083740234   0.004809198081493
%    0.001680015027523   0.765760116577148   0.130840768814087   0.005303588509560
%    0.001113213151693   0.682846679687500   0.130632925033569   0.006153529286385
%    0.002568353414536   0.754957885742187   0.254111194610596   0.005207121968269
%    0.001236899867654   0.662300720214844   0.165663166046143   0.006506183147430];


% param =[
% 
%    0.173397436738014  72.840103149414062  17.492570877075195   0.562047004699707
%    0.075591310858727  44.855827331542969   5.307141304016113   0.548380494117737
%    0.057756066322327  66.363334655761719   5.256161212921143   0.670129954814911];

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

numTasks = 1; %one particle == one task
tasksPerJob = 1:numTasks:size(pksNew,1);
for j = 1:length(tasksPerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue
    disp(sprintf('\n------------job%02d, created with the following tasks------------',j));

    startTask = tasksPerJob(j);
    stopTask = startTask+numTasks-1; %subtract to account for MATLAB indexing displacement
   
    if(stopTask > size(pksNew,1)) 
        stopTask = size(pksNew,1);
    end
    
    for np = startTask:stopTask %slice == task #
        [v,center] = getCutoutSectVertices(pksNew(np,:),0.4*D(np),size(InormRed)); %vertices of rectangular portion of 'InormRed' that is (2*Dactual)^3 and centered on pksNew(np,:) 
        createTask(job(j),@getPosition,2,{v(1,:),center,K(np),D(np),w(np),offset(np),InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3))});
        disp(sprintf('locating pksNew(%d,:) at pixel resolution...',np));
    end
end

% submit one job at a time to be computed by all available slave nodes
for j = 1:length(job)
    submit(job(j));
    disp(sprintf('job%02d, submitted',j))
end 
% Actual time for 0.6*D(np):
%   first job started at: Mon Aug 03 19:36:33 EDT 2015
%   last job finished at: Mon Aug 03 22:08:14 EDT 2015
%   total time: 2hr31min41sec for locating 420 particles


% Actual time for 0.4*D(np):
%   first job started at: Tue Aug 04 15:11:51 EDT 2015
%   last job finished at: Tue Aug 04 16:21:28 EDT 2015
%   total time: 1hr09min37sec for locating 421 particles

pksNewer = zeros(size(pksNew),'single');
Rsquared = zeros([8,size(pksNew,1)],'single'); 
for j = 1:length(job)
    state = job(j).State;
    
    startTask = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        try 
            results = fetchOutputs(job(j));
            disp(sprintf('\n------------job%02d, COMPLETED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksNewer(startTask+np-1,:) = results{np,1}; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = results{np,2};
                disp(sprintf('pksNew(%d,:) becomes pksNewer(%d,:)',startTask+np-1,startTask+np-1));
            end
        catch err 
            disp(sprintf('\n------------job%02d, FAILED with the following tasks------------',j));
            for np = 1:numel(job(j).Tasks)%fetchOutputs(job(j)) doesn't work so you have to read the number of tasks using numel(job(j).Tasks)
                pksNewer(startTask+np-1,:) = NaN; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = NaN;
                disp(sprintf('pksNewer(%d,:) = NaN',startTask+np-1));
            end
        end%end of try
    end
end

temp = isnan(pksNewer(:,1));
npRedo = find(temp(:) == 1); %There aren't any!

% pksNewer =[
% 
%          624         515         134
%          414         326         680
%          171         182         353
%           42          29         457
%          127         326         679
%           71         269         827
%           80         177         823
%          108         815         804
%            1         975         561
%         1023         123         398
%          805         378         399
%          634         115         623
%          642         176         674
%          663         415         687
%          903         431         708
%          578         384         712
%          704         125         809
%          912         882         288
%          743         591         282
%          622         905         606
%          896         703         633
%          728         851         669
%         1015         733         746];

% pksNewer =[
% 
%     28    15   283
%     84    73   303
%    342     2   329];
   
for np = 1:size(pksNew,1)
    startSlice = nanmin([pksNewer(np,3), pksNew(np,3)])-45;
    stopSlice = nanmax([pksNewer(np,3), pksNew(np,3)])+45;

    if(startSlice < 1)
        startSlice = 1;
    end

    if(stopSlice > size(InormRed,3))
        stopSlice = size(InormRed,3);
    end
    
    for slice = startSlice:stopSlice
        figure(1)
        subplot(1,2,1)
        hold off
        simage(InormRed(:,:,slice));
        hold on
        plot(pksNewer(np,2),pksNewer(np,1),'k*','markersize',20)
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNewer(%d,:) = [%d, %d, %d] (k*)\npksNew(%d of %d,:) = [%d, %d, %d](g*)',np,pksNewer(np,1),pksNewer(np,2),pksNewer(np,3),np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)

        subplot(1,2,2)
        hold off
        simage(InormGreen(:,:,slice));
        hold on
        plot(pksNewer(np,2),pksNewer(np,1),'k*','markersize',20)
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNew(%d of %d,:) = [%d, %d, %d](k*)',np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)
        pause
    end
end

% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('single_pksNewer4thPass_id_row_col_slice.txt','w');
% for np = 1:size(pksNewer,1)
%     fprintf(fid,'%d %d %d %d\n',np,pksNewer(np,1),pksNewer(np,2),pksNewer(np,3));
% end
% fclose(fid);

% -----(start: locate position of particle at sub-pixel resolution)-------

clearvars -except pksNew Inorm* Kactual Dactual wActual offsetActual 

% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('single_pksNewer3rdPass_id_row_col_slice.txt','r');
% temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
% fclose(fid);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single_pksNewer4thPass_id_row_col_slice.txt','r');
temp = fscanf(fid,'%d %f %f %f\n',[4,Inf]);
fclose(fid);

order = temp(1,:)';
pksNewer(order,1) = temp(2,:)';
pksNewer(order,2) = temp(3,:)';
pksNewer(order,3) = temp(4,:)';


for np = 1:size(pksNew,1)
    startSlice = nanmin([pksNewer(np,3), pksNew(np,3)])-45;
    stopSlice = nanmax([pksNewer(np,3), pksNew(np,3)])+45;

    if(startSlice < 1)
        startSlice = 1;
    end

    if(stopSlice > size(InormRed,3))
        stopSlice = size(InormRed,3);
    end
    
    for slice = startSlice:stopSlice
        figure(1)
        subplot(1,2,1)
        hold off
        simage(InormRed(:,:,slice));
        hold on
        plot(pksNewer(np,2),pksNewer(np,1),'k*','markersize',20)
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNewer(%d,:) = [%d, %d, %d] (k*)\npksNew(%d of %d,:) = [%d, %d, %d](g*)',np,pksNewer(np,1),pksNewer(np,2),pksNewer(np,3),np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)

        subplot(1,2,2)
        hold off
        simage(InormGreen(:,:,slice));
        hold on
        plot(pksNewer(np,2),pksNewer(np,1),'k*','markersize',20)
        plot(pksNew(np,2),pksNew(np,1),'g*','markersize',20)
        title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksNew(%d of %d,:) = [%d, %d, %d](k*)',np,size(pksNew,1),pksNew(np,1),pksNew(np,2),pksNew(np,3)),'fontsize',20)
        pause
    end
end


clear pksNew

numTry = 10; %# of tries to calculate param(np,:)

param = zeros([size(pksNewer,1), 4]);
Rsquared = zeros([size(pksNewer,1), 1],'single');

tic;
for np = 1:size(pksNewer,1)
    [v,center] = getCutoutSectVertices(pksNewer(np,:),0.75*Dactual,size(InormRed));
    cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
    avgProfile = getSphericalProfile(center,Dactual,wActual,cutoutSect);
    
    b = ceil((length(avgProfile)-1)/2);
   
    idxMin1 = find(avgProfile == min(avgProfile(1:b)));
    dist1 = idxMin1;

    idxMin2 = find(avgProfile == min(avgProfile(b:numel(avgProfile))));
    dist2 = numel(avgProfile)-idxMin2;

    startIdx = floor((dist1+dist2)/2);
    if(startIdx == 0)
        startIdx = 1;
    end
    stopIdx = numel(avgProfile)-startIdx;

    %bNew = ceil((numel(avgProfile(startIdx:stopIdx))-1)/2);
    %rNew = linspace(-bNew,bNew,numel(avgProfile(startIdx:stopIdx)));

    [param(np,:), Rsquared(np)] = getParameters(avgProfile(startIdx:stopIdx),Dactual,wActual,numTry);
    
    if(~isnan(param(np,1)))
        disp(sprintf('SUCCESSFUL computation of parameters for pksNewer(%d,:)',np));
    end
    
    if(isnan(param(np,1)))
        disp(sprintf('FAILED computation of parameters for pksNewer(%d,:)',np));
    end
end
elapsedTime = toc; %16.5133669seconds for 285 particles

K = param(:,1);
D = param(:,2);
w = param(:,3);
offset = param(:,4);

% param =(1.0e+02)*[
% 
%    0.002391418665648   0.636307106018066   0.161573524475098   0.004709480702877
%    0.002314427345991   0.850484771728516   0.219881629943848   0.005327161550522
%    0.002193623483181   0.644549789428711   0.189459705352783   0.004794358015060
%    0.001761718839407   0.651306838989258   0.217511157989502   0.005684220790863
%    0.005630346536636   1.288073272705078   0.354563293457031   0.002120949774981
%    0.001177412644029   0.676541137695313   0.119641342163086   0.006036115288734
%    0.001351763606071   0.612169570922852   0.135741386413574   0.006167746186256
%    0.002078406363726   0.720023193359375   0.217848682403564   0.005466148257256
%    0.003580550849438   1.134198226928711   0.481642684936523   0.004021507501602
%    0.001953471601009   0.828086166381836   0.291907768249512   0.005672151446342
%    0.002175459563732   0.738300704956055   0.163036670684814   0.005252463221550
%    0.001073420271277   0.579455413818359   0.128640747070313   0.005848444104195
%    0.001751979291439   0.795834503173828   0.183925952911377   0.005471484661102
%    0.000877207741141   0.534456977844238   0.093492650985718   0.005403815507889
%    0.002415401488543   0.875793304443359   0.226175212860107   0.005307715535164
%    0.002209520787001   0.776091461181641   0.182998867034912   0.005441777110100
%    0.001637811958790   0.742945098876953   0.181869163513184   0.005664714574814
%    0.002196521610022   0.722824249267578   0.157049303054810   0.005537425279617
%    0.002124168127775   0.629904098510742   0.168019104003906   0.004877961874008
%    0.002070362567902   0.786465835571289   0.192194595336914   0.005195144414902
%    0.002798812091351   0.899332885742188   0.294453277587891   0.004906968772411
%    0.002228068858385   0.788183441162109   0.157313184738159   0.005280736088753
%    0.001753545254469   0.755720520019531   0.278544960021973   0.006068524122238];

% param =[
% 
%    0.114578612148762  53.218952178955078  10.549324035644531   0.621265769004822
%    0.104685977101326  48.793304443359375  10.447098731994629   0.542110979557037
%    0.078063420951366  69.125640869140625   8.880936622619629   0.658804655075073];

clusterD = parcluster(); %uses default profile 
delete(clusterD.Jobs); %delete jobs stored in /home/eru/.matlab/local_cluster_jobs/R2013b to free up hard disk memory

rangeColormap = 255;

numTasks = 1; %one particle == one task
tasksPerJob = 1:numTasks:size(pksNewer,1);
for j = 1:length(tasksPerJob)
    job(j) = createJob(clusterD); %each job is stored in the MATLAB session so I think its good memory allocation to create only one job with many tasks then delete the job after its complete to free of memory on the MJS queue
    disp(sprintf('\n------------job%02d, created with the following tasks------------',j));

    startTask = tasksPerJob(j);
    stopTask = startTask+numTasks-1; %subtract to account for MATLAB indexing displacement
   
    if(stopTask > size(pksNewer,1)) 
        stopTask = size(pksNewer,1);
    end
    
    for np = startTask:stopTask %slice == task #
        particleRadius = ((D(np)/2)+(2*w(np))); %this is the apperant particle radius
        [v,center] = getCutoutSectVertices(pksNewer(np,:),particleRadius,size(InormRed)); %vertices of rectangular portion of 'Inorm' that is (2*Dactual)^3 and centered on validPks(np,:) 
        cutoutSect = InormRed(v(1,1):v(7,1),v(1,2):v(7,2),v(1,3):v(7,3));
        centerIdx = sub2ind(size(cutoutSect),center(1),center(2),center(3));
        createTask(job(j),@getPositionSubPixel,2,{v(1,:),centerIdx,K(np),D(np),w(np),offset(np),cutoutSect,rangeColormap});
        disp(sprintf('locating pksNewer(%d,:) at SUB-pixel resolution...',np));
    end
end

% submit one job at a time to be computed by all available slave nodes
for j = 1:length(job)
    submit(job(j));
    disp(sprintf('job%02d, submitted',j))
end 

% 1st job to start, started at: Wed Aug 12 19:49:27 EDT 2015
% last job to finish, finished at: Wed Aug 12 20:21:06 EDT 2015
% total time: 30min54seconds for 23 particles

pksBest = zeros(size(pksNewer),'single');
Rsquared = zeros([11,size(pksNewer,1)],'single'); 
for j = 1:length(job)
    state = job(j).State;
    
    startTask = tasksPerJob(j);

    if(state(1:3) == 'fin')%'fin' for finished 
        try 
            results = fetchOutputs(job(j));
            disp(sprintf('\n------------job%02d, COMPLETED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksBest(startTask+np-1,:) = results{np,1}; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = results{np,2};
                disp(sprintf('pksNewer(%d,:) becomes pksBest(%d,:)',startTask+np-1,startTask+np-1));
            end
        catch err 
            disp(sprintf('\n------------job%02d, FAILED with the following tasks------------',j));
            for np = 1:size(results,1)
                pksBest(startTask+np-1,:) = NaN; %minus 1 becuase MATLAB indexing starts at 1
                Rsquared(:,startTask+np-1) = NaN;
                disp(sprintf('pksNewer(%d,:) = NaN',startTask+np-1));
            end
        end%end of try
    end
end

temp = isnan(pksBest(:,1));
npRedo = find(temp(:) == 1); %There aren't any!

% pksBest =(1.0e+03)*[
% 
%    0.6233019   0.5159961   0.1330000
%    0.4150000   0.3215804   0.6790039
%    0.1700118   0.1829882   0.3520039
%    0.0439882   0.0300000   0.4673294
%    0.1260118   0.3270000   0.6799843
%    0.0719961   0.2680000   0.8260000
%    0.0809412   0.1775725   0.8220000
%    0.1070078   0.8140471   0.8030000
%    0.0020000   0.9740157   0.5605333
%    1.0220000   0.1239961   0.3970784
%    0.8043333   0.3786745   0.3980039
%    0.6339020   0.1110196   0.6210000
%    0.6410510   0.1750000   0.6733059
%    0.6638902   0.4151686   0.6860039
%    0.8989961   0.4270157   0.7209882
%    0.5774000   0.3849961   0.7110078
%    0.7030314   0.1180000   0.8080000
%    0.9110157   0.8810236   0.2874039
%    0.7420039   0.5916196   0.2818863
%    0.6210078   0.9059882   0.6050000
%    0.8969686   0.7039882   0.6320039
%    0.7284117   0.8500000   0.6680000
%    1.0140000   0.7339647   0.7450039];

% pksBest4 =(1.0e+02)*[
% 
%    0.2898824   0.1599216   2.8347842
%    0.8300784   0.7201569   3.0399609
%    3.6905096   0.0199608   3.0501569];


% ---------------(check pksBest3 located in 3rd Pass)---------

clearvars -except Inorm* Kactual Dactual wActual offsetActual

cd('~/poincareProgs/')
run results12may2015_6may2015_3_1920A1000r_13xS3

pksBest = double(pksBest);

[dist, pksOrder] = orderPeaks(pksBest);

% view multimers:
for np = (size(pksBest,1)-size(pksBest4,1)+1):size(pksBest,1)
    npp = pksOrder(1:5,np);

    [naught, idx] = sort(pksBest(npp,3));
    npp = npp(idx);
    
%     %for 2mers:
%     npp = union(pksOrder(1:3,np),[np:(np+1)]);
    
%     npp = pksOrder(1:3,np);
    
%     %for 5mers:
%     npp = union(pksOrder(1:6,np),[np:(np+4)]);

%     %for 6mers
%     npp = union(pksOrder(1:6,np),[np:(np+5)]);
    
    clear label 
    disp(sprintf('---------------------------x-----------------------------'));
    for jj = 1:numel(npp)
        label{jj} = num2str(sprintf('%d',npp(jj)));
        disp(sprintf('\tpksBest(%d,:) = [%.03f, %.03f, %.03f], distance(%d,%d) = %.04f',npp(jj),pksBest(npp(jj),1),pksBest(npp(jj),2),pksBest(npp(jj),3),npp(jj),np,dist(npp(jj),np)));
    end
    disp(sprintf('---------------------------x-----------------------------\n'));
   
    
    startSlice = round(min(pksBest(npp,3))-10);
    stopSlice = round(max(pksBest(npp,3))+10);

    if(startSlice < 1)
        startSlice = 1;
    end

    if(stopSlice > size(InormRed,3))
        stopSlice = size(InormRed,3);
    end
    
    for slice = startSlice:stopSlice
        figure(1)
        subplot(1,2,1)
        hold off
        simage(InormRed(:,:,slice));
        hold on
        plot(pksBest(npp,2),pksBest(npp,1),'k*','markersize',20)
        plot(pksBest(np,2),pksBest(np,1),'g*','markersize',20)
        text(pksBest(npp,2),pksBest(npp,1),label,'fontsize',12)
        title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksBest(%d of %d,:) = [%.04f, %.04f, %.04f] \n %d particles plotted',np,size(pksBest,1),pksBest(np,1),pksBest(np,2),pksBest(np,3),numel(npp)),'fontsize',20)

        subplot(1,2,2)
        hold off
        simage(InormGreen(:,:,slice));
        hold on
        plot(pksBest(npp,2),pksBest(npp,1),'k*','markersize',20)
        plot(pksBest(np,2),pksBest(np,1),'g*','markersize',20)
        text(pksBest(npp,2),pksBest(npp,1),label,'fontsize',12)
        title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)
        xlabel(sprintf('pksBest(%d,3) = %.03f \npksBest(%d,3) = %.03f \npksBest(%d,3) = %.03f',npp(1),pksBest(npp(1),3),npp(2),pksBest(npp(2),3),npp(3),pksBest(npp(3),3)),'fontsize',20)

        pause
    end
end

% % save pksBest from 3rd pass:
% cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
% fid = fopen('single3rdPass12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
% for np = 1:size(pksBest3,1)
%     fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest3(np,1),pksBest3(np,2),pksBest3(np,3),K3(np),D3(np),w3(np),offset3(np));
% end
% fclose(fid);

% save pksBest from 3rd pass:
cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('single4thPass12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for np = 1:size(pksBest4,1)
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest4(np,1),pksBest4(np,2),pksBest4(np,3),K4(np),D4(np),w4(np),offset4(np));
end
fclose(fid);

% double-check that ALL particles have been found:
for slice = 100:size(InormRed,3)
    npTemp1 = find(round(pksBest(:,3)) >= (slice-40));
    npTemp2 = find(round(pksBest(:,3)) <= (slice + 40));
    npp = intersect(npTemp1,npTemp2);
    
    clear label 
    disp(sprintf('---------------------------x-----------------------------'));
    for jj = 1:numel(npp)
        label{jj} = num2str(sprintf('%d',npp(jj)));
        disp(sprintf('\tpksBest(%d,:) = [%d, %d, %d]',npp(jj),pksBest(npp(jj),1),pksBest(npp(jj),2),pksBest(npp(jj),3)));
    end
    disp(sprintf('---------------------------x-----------------------------\n'));  
    
    figure(1)
    subplot(1,2,1)
    hold off
    simage(InormRed(:,:,slice));
    hold on
    plot(pksBest(npp,2),pksBest(npp,1),'k*','markersize',20)
    text(pksBest(npp,2),pksBest(npp,1),label,'fontsize',12)
    title(sprintf('InormRed(:,:,%d)',slice),'fontsize',20);

    subplot(1,2,2)
    hold off
    simage(InormGreen(:,:,slice));
    hold on
    plot(pksBest(npp,2),pksBest(npp,1),'k*','markersize',20)
    text(pksBest(npp,2),pksBest(npp,1),label,'fontsize',12)
    title(sprintf('InormGreen(:,:,%d)',slice),'fontsize',20)

    pause
end


clearvars -except pksBest K D w offset Inorm* Kactual Dactual wActual offsetActual

% save ALL particles:
npDimer = [378:421,545:548]; %There's ONLY dimers
npSingle = 1:size(pksBest,1);
npSingle(npDimer) = [];

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('ALLsingle12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for jj = 1:numel(npSingle)
    np = npSingle(jj);
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest(np,1),pksBest(np,2),pksBest(np,3),K(np),D(np),w(np),offset(np));
end
fclose(fid);

cd('~/poincareProgs/analysis/analysis12may2015_6may2015_3_1920A1000r_13xS3/')
fid = fopen('ALL2mers12may2015_6may2015_3_1920A1000r_13xS3_id_row_col_slice_K_D_w_offset.txt','w');
for jj = 1:2:numel(npDimer)
    np = npDimer(jj);
    np2 = npDimer(jj+1);
    
    fprintf(fid,'%d %f %f %f %f %f %f %f\n',np,pksBest(np,1),pksBest(np,2),pksBest(np,3),K(np),D(np),w(np),offset(np));
    fprintf(fid,'%d %f %f %f %f %f %f %f\n\n',np2,pksBest(np2,1),pksBest(np2,2),pksBest(np2,3),K(np2),D(np2),w(np2),offset(np2));
end
fclose(fid);
