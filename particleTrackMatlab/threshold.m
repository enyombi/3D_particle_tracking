function pos = threshold(A,fracOfMax,sizeFullScaleImg,displacement)
% Objective: returns the positions pixels located by thresholding. These 
% are the positions in full-sized 3D-image NOT a smaller section of this
% image hence 'displacement'

% A = partitioned section of the convolution image 'cScaled' 
% 
% fracOfMax = fraction of max pixel found in 'A' that will be used as 
%             cutoff/threshold value
% 
% sizeFullScaleImg = the size of cScaled
% 
% displacment = first index of 'A' identified by its original location in
%               'cScaled' (i.e., these are indices identifying the pixel
%               in the top-front-left corner of section from 'cScaled'
% 
% e.g., displacement for Task1 = [1 1 1]
%       displacement for Task2 = [1 167 1]
%       displacement for Task2 = [1 334 1]
% 
% cScaled = 
% 
%    /....                                                 /
%   /_____________________________________________________/
%  /       /       /       /       / .....               /|
% /_______/_______/_______/_______/_____________________/ |
% |       |       |                                     | |
% | Task1 | Task1025                                    | |
% |_______|_______|                                     | |
% |       |       |                                     | |
% | Task2 | Task1026                                    | |
% |_______|_______|                                     | |
% |       |       |                                     | |
% | Task3 | Task1027                                    | |
% |_______|_______|                                     | |
% |                                                     | |
% | .....                                               | |  /
% |_______ _____________________________________________| | / 
% |       |       |                              |      | |/  
% | Task1024 Task2048 ....                       | Task(1024*1024)  
% |_______|_______|______________________________|______|/
% 
% note: sections/Tasks of cScaled are numbered the same that MATLAB
% linearly indexes a matix
cutoffVal = fracOfMax*max(A(:));
idxPos = find(A(:) >= cutoffVal); %position of max pixels in A. 

% Note: find() sequentially goes down each column of 'A' and retuns the
% max. values in order. This means that max pixel value on the 
% edges of a particle are listed before those identifying the center of a
% particle. The next step, i.e., getNeighborsValid().m keeps the first 
% pixel position on the edge and eliminates all subsequent max pixel 
% positions that more accurately identify a particle.  To be rid of the
% nuisance, randomize the max pixel positions...
idxPos = idxPos(randperm(length(idxPos)));

[rowA, colA, sliceA] = ind2sub(size(A),idxPos); %position of max pixels in A

%position of max pixels on cScaled (note:subtract 1 because matlab indexing starts at 1)
rowScaled = rowA+(displacement(1)-1);
colScaled = colA+(displacement(2)-1); 
sliceScaled = sliceA+(displacement(3)-1);

pos = sub2ind(sizeFullScaleImg,rowScaled,colScaled,sliceScaled); %liner-index position of max pixels on cScaled
end %threshold()
















%Matlab's numbering of indices for a matrix
% a = zeros([3 3 3]);

% % a(:,:,1) =
% % 
% %      1     4     7
% %      2     5     8
% %      3     6     9
% % 
% % 
% % a(:,:,2) =
% % 
% %     10    13    16
% %     11    14    17
% %     12    15    18
% % 
% % 
% % a(:,:,3) =
% % 
% %     19    22    25
% %     20    23    26
% %     21    24    27

% ----------------------------------------------------------



% this version of matlab has a glitch that makes it always print out a
% matrix with the same dimensions as its indices
% a = magic(3);
% r = [1 1 1];
% c = [2 3 1];

% a(r,c)
% 
% ans =
% 
%      1     6     8
%      1     6     8
%      1     6     8
%      
% %..but
% [a(1,2) a(1,3) a(1,1)]
% 
% ans =
% 
%      1     6     8

% Test threshold.m in parallel 

% a = magic(3);
% 
% clusterD = parcluster(); %uses default profile 
% job = createJob(clusterD);
% createTask(job,@threshold,1,{magic(3) 0.5});
% submit(job);
% wait(job);
% [peaks val]= fetchOutputs(job); %max pixel value from each cutout section of cScaled
