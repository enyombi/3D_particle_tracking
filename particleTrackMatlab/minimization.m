function [spx_min, spy_min, spz_min, min_value] = minimization(peak,test_dist,real_img,ideal_particle)

%Author: Eru Kyeyune-Nyombi
%18-Feb-2013
%filename: minization

%Primary objective: To locate particle centers and pixel resolution

%Objective: implement minization technique.  Move peaks found by
%thresholding until the difference between the ideal particle mask and the
%real particle centered around these peaks is minimized. Finally, return
%the location (i.e., row, column, and slice) of the pixel where the
%sum of the the differences is minimumized.  This location is the 
%location of the particle center

%peak = peak located by thresholding the original image
%
%test_dist = the distance in pixels away from any peak. The minization
%            routine will be conducted with the ideal particle centered
%            at any of the pixels in the box that is 
%            (2*test_dist)x(2*test_dist)x(2*test_dist) and 
%            that surround 'peak'
%
%real_img = the original 3D-image
%
%ideal_particle = the 3D-image of the ideal particle


%Begin by assuming that all of the peaks returned by thresholding are no
%more than a distance D/6 away from the real particle center. 
%
%Therefore, the actual particle center is somewhere in a rectangular 
%box centered around 'peak' that is 2*test_dist pixels long, wide, and 
%tall and also centered 

%Find the difference between a cutout of the real image centered at any 
%point(i.e., pixel) within this box and the ideal particle image.  
%Finally, calculate and record the sum_of_squares for each 
%location in the corresponding element of the box called 'R_squared'
%
%Summary: Each element in 'R_squared' is the sum of the differences-squared
%between pixel values in a cutout section of the real image and 
%the ideal particle mask image (aka 'ideal_particle') squared. 
%An element with a low value in 'R_squared' means that the cutout 
%section of the real image centered at this location
%is very similar to the ideal particle mask image.  This in turn means,
%that a particle is centered at that location

ideal_img_squared = ideal_particle.^2;
diff_squared_max = sum(ideal_img_squared(:)); %this is the max sum-of-squares value for (ideal_particle - zeros(size(ideal_particle)))
background_value = diff_squared_max.^2; %use a value much larger than the 'diff_squared_max' so you can clearly see which pixels have a value due to (ideal_particle - zeros(size(ideal_particle)))
R_squared = background_value*ones(size(real_img));

dim = size(real_img);

spx_range = peak(1) + (-test_dist:test_dist);
spy_range = peak(2) + (-test_dist:test_dist);
spz_range = peak(3) + (-test_dist:test_dist);

%retain only those coordinates that are within the boudaries of the 
%3D-image
spx_range = spx_range((spx_range <= dim(1)) & (spx_range >= 1));
spy_range = spy_range((spy_range <= dim(2)) & (spy_range >= 1));
spz_range = spz_range((spz_range <= dim(3)) & (spz_range >= 1));

%jjj, kkk, and lll are the number of rows, columns, and slice
%(respectively) in the ideal particle mask image (i.e., 'ideal_particle).
%A value of zero denotes the center of the ideal particle mask image.
jjj = -fix(size(ideal_particle,1)/2):fix(size(ideal_particle,1)/2);
kkk = -fix(size(ideal_particle,2)/2):fix(size(ideal_particle,2)/2);
lll = -fix(size(ideal_particle,3)/2):fix(size(ideal_particle,3)/2);

%Technically, 'center_ideal_particle' is NOT the center of the spherical
%image. Instead, 'center_ideal_particle' is the center of 3D rectuanglar
%matrix that contains the spherical image of an ideal particle
center_ideal_particle = fix((size(ideal_particle)/2)+1); %use fix() to round to the nearest whole number which is the index for the row, column and slice

spx_ideal_particle = center_ideal_particle(1);
spy_ideal_particle = center_ideal_particle(2);
spz_ideal_particle = center_ideal_particle(3);

for iter1 = 1:length(spx_range)
    for iter2 = 1:length(spy_range)
        for iter3 = 1:length(spz_range)
            spx_test = spx_range(iter1);
            spy_test = spy_range(iter2);
            spz_test = spz_range(iter3);

            %ll, kk, and jj are the rows, columns, and slices respectively that are
            %centered on the pt. (spx_test,spy_test,spy_test)
            jj = round(spx_test) + jjj;
            kk = round(spy_test) + kkk;
            ll = round(spz_test) + lll;

            %the following edits jj,kk, and ll to only include elements that are
            %contained within the boundaries/dimensions of the 3D-image. 
            jj = jj((jj<=dim(1))&(jj>=1));
            kk = kk((kk<=dim(2))&(kk>=1));
            ll = ll((ll<=dim(3))&(ll>=1));
            
            %Use the center of the ideal particle as a reference point.  
            %And defining the rows, columns, and slices around this 
            %point using xx, yy, and zz

            front_rows = jj(jj<spx_test);
            back_rows = jj(jj>spx_test);
            xstart = spx_ideal_particle - numel(front_rows);
            xend = spx_ideal_particle + numel(back_rows);
            x_front = xstart:(spx_ideal_particle - 1);%finds all rows in front of the center row (i.e., 'spx_ideal_particle') for the ideal mask particle (i.e., 'ideal_particle')
            x_back = (spx_ideal_particle + 1):xend;%finds all rows behind the center row (i.e., 'spx_ideal_particle') for the ideal mask particle (i.e., 'ideal_particle')
              
            xx = horzcat(x_front, spx_ideal_particle, x_back); %'xx' is all the rows in front of, behind, and including the center row 

            front_cols = kk(kk<spy_test);
            back_cols = kk(kk>spy_test);
            ystart = spy_ideal_particle - numel(front_cols);
            yend = spy_ideal_particle + numel(back_cols);
            y_front = ystart:(spy_ideal_particle - 1);%finds all columns in front of the center column (i.e., 'spy_ideal_particle') for the ideal mask particle (i.e., 'ideal_particle')
            y_back = (spy_ideal_particle + 1):yend;%finds all columns behind the center column (i.e., 'spy_ideal_particle') for the ideal mask particle (i.e., 'ideal_particle')
            
            yy = horzcat(y_front, spy_ideal_particle, y_back); %'yy' is all the columns in front of, behind, and including the center column
            
            front_slices = ll(ll<spz_test);
            back_slices = ll(ll>spz_test);
            zstart = spz_ideal_particle - numel(front_slices);
            zend = spz_ideal_particle + numel(back_slices);
            z_front = zstart:(spz_ideal_particle - 1);%finds all columns in front of the center column (i.e., 'spy_ideal_particle') for the ideal mask particle (i.e., 'ideal_particle')
            z_back = (spz_ideal_particle + 1):zend;%finds all columns behind the center column (i.e., 'spy_ideal_particle') for the ideal mask particle (i.e., 'ideal_particle')
      
            zz = horzcat(z_front, spz_ideal_particle, z_back); %'zz' is all the slices in front of, behind, and including the center slice
    
            diff = real_img(jj,kk,ll) - ideal_particle(xx,yy,zz);
            diff_squared = diff.^2;
            R_squared(spx_test,spy_test,spz_test) = sum(diff_squared(:));
        end
    end
end

min_value = min(min(min(R_squared)));

idx_min = find(R_squared == min_value);

[spx_min, spy_min, spz_min] = ind2sub(size(R_squared),idx_min);
