function displacement = minimize_error(peaks,real_img,D,w)

% Primary objective: to locate particle centers at subpixel resolution

%objective: To minize error between the calculated image constructed using
%'peaks' and the real image and calculate how distant the peaks are from
%the centers of particles in the real image (note: this distance is denoted
%using 'dp')

%dp = particle centers in real image - particle centers in calculated image
%Note: the particle centers in the real image are unknown

%Author: Eru K.
%Date: 15-April-2013 (adopted from 'hist_1apr2013.m')
%19Oct2013, Modified for handling multiple peaks with different D and w values,

[over r x y z] = peak_placement(peaks,size(real_img),D,w);

r=abs(r)+eps; %'eps' is very nearly zero. Adding 'eps' prevents division by 0 in subsequent steps
mask = ipf_eru(r,D,w,over);

diff = real_img-mask;

%IMPORTANT CLARIFICATION: The pixel value in regions of the image not 
%occupied by particles is inconsequential because this algorithm
%only computes dp for regions of the image occupied by particles see last
%for-loop where 'idx = find(over = np)' is used

%Initialize:
mask_x = zeros(size(real_img)); 
mask_y = mask_x;
mask_z = mask_x;
mask_xx = mask_x;
mask_yy = mask_x;
mask_zz = mask_x;
mask_xy = mask_x;
mask_xz = mask_x;
mask_yz = mask_x;

for np = 1:size(peaks,1)
    %tanh1, sech2, etc. are for computation over the entire image
    %using the parameters for only one particle np...
    %Note: the values of all pixels changes with each iteration
    %   e.g, 
    %       when np = 1, tanh1, sech2, denom1, etc. == pixel values 
    %       for D(1) and w(1)
    % 
    %       when np = 2, tanh1, sech2, denom1, etc. == pixel values 
    %       for D(2) and w(2), the previous values for np = 1 are replaced
    %
    %       and so on...with pixel values being replaced each iteration
    
    tanh1 = tanh((r-(D(np)/2))/w(np));
    sech2 = (sech((r-(D(np)/2))/w(np))).^2;
    denom1 = (2*w(np)).*r;
    denom2 = ((r.^3)+eps).*(w(np)^2); 
    denom3 = ((r.^4)+eps).*(w(np)^2);
    
    %'idx' locates the region of tanh1, sech2, that corresponds to particle
    %np
    idx = find(over == np);

    %1st derivatives of the mask ideal particle image
    mask_x(idx) = (-x(idx).*sech2(idx))./denom1(idx);
    mask_y(idx) = (-y(idx).*sech2(idx))./denom1(idx);
    mask_z(idx) = (-z(idx).*sech2(idx))./denom1(idx);
    
    %2nd derivatives of the the mask ideal particle image
    mask_xx(idx) = (sech2(idx).*((-0.5*w(np).*r(idx).^2)+(0.5.*w(np).*x(idx).^2)+(r(idx).*(x(idx).^2).*tanh1(idx))))./denom2(idx);
    mask_yy(idx) = (sech2(idx).*((-0.5*w(np).*r(idx).^2)+(0.5.*w(np).*y(idx).^2)+(r(idx).*(y(idx).^2).*tanh1(idx))))./denom2(idx);
    mask_zz(idx) = (sech2(idx).*((-0.5*w(np).*r(idx).^2)+(0.5.*w(np).*z(idx).^2)+(r(idx).*(z(idx).^2).*tanh1(idx))))./denom2(idx);

    mask_xy(idx) = ((x(idx).*y(idx).*sech2(idx)).*(0.5.*r(idx).*w(np)+((r(idx).^2).*tanh1(idx))))./denom3(idx);
    mask_xz(idx) = ((x(idx).*z(idx).*sech2(idx)).*(0.5.*r(idx).*w(np)+((r(idx).^2).*tanh1(idx))))./denom3(idx);
    mask_yz(idx) = ((y(idx).*z(idx).*sech2(idx)).*(0.5.*r(idx).*w(np)+((r(idx).^2).*tanh1(idx))))./denom3(idx);
end

%chi = ?[real_img - mask(x,y,z)]^2 == 'sum of residual error squared' 
% note: really, mask = f(x,y,z,D,w) but you are assuming D and w are const.

%1st derivatives of chi 
chix = diff.*mask_x;
chiy = diff.*mask_y;
chiz = diff.*mask_z;

%2nd derivatives of chi
chixx = (mask_x.^2) - (diff.*mask_xx);
chixy = (mask_x.*mask_y) - (diff.*mask_xy);
chixz = (mask_x.*mask_z) - (diff.*mask_xz);

chiyy = (mask_y.^2) - (diff.*mask_yy);
chiyz = (mask_y.*mask_z) - (diff.*mask_yz);

chizz = (mask_z.^2) - (diff.*mask_zz);

%compute displacement,'dp', for each peak

%intialize:
dp = zeros(size(peaks));
for np = 1:size(peaks,1)
    %np = peak/particle number
    idx = find(over == np); %locate all pixels in the image that belong to a given peak using 'over'

    A = zeros(3,3);
    A(1,1) = sum(chixx(idx));
    A(1,2) = sum(chixy(idx));
    A(1,3) = sum(chixz(idx));

    A(2,1) = A(1,2);
    A(2,2) = sum(chiyy(idx));
    A(2,3) = sum(chiyz(idx));

    A(3,1) = A(1,3);
    A(3,2) = A(2,3);
    A(3,3) = sum(chizz(idx));

    b = [sum(chix(idx)), sum(chiy(idx)), sum(chiz(idx))];

    dp(np,:) = b*pinv(A);
end

displacement = dp;