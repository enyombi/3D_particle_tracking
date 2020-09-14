function ipi=ssf(param,cr)
% filename: ssf.m, 'sphere spread function'
% Adopted from Shattuck ipf().m 
%
% Objective: To calculates an ideal particle image ipi.  
% The particle has diameter D and
% width parameter w.  2w is the width of 76% of the fall off.

% Note: this function should be called via getCalcImg.m to ensure that the
% particle is drawn only in the spherical region that it occupies.

K = param(1); %prefactor to adjust the max height of the ssf
D = param(2); %particle diameter
w = param(3); %width of particle's hazy edges 
offset = param(4); %lowest value of particle's ssf

ipi=(K*((1-tanh((abs(cr)-D/2)/w))/2))+offset;
end

% %---------------- Successfully tested ----------------
% % In 1Dimension:
% x = -5:5;
% yExp = ipf([4,1.1],x); %D = 4, w =1
% 
% % now assume that you don't now D and w.  Use estimates D~5, w~1, and
% % nlinfit() to calculate them
% beta = nlinfit(x,yExp,'ipfA',[5 1]);
% % beta =
% % 
% %     4.0000    1.1000

% % In 2Dimensions:
% [x y] = ndgrid(-50:50,-50:50);
% r = sqrt(x.^2+y.^2);
% yExp = ipfA([20,4],r);
% 
% beta = nlinfit(r,yExp,'ipfA',[17 1]);
% Error: using nlinfit (line 101) Requires a vector second input argument.

% % nlinfit() only works for 1Dimension!
% % 
% % Recourse: 
% % 1) You either sample across a 3Dimensional particle to obtain linear
% % profiles to be used in nlinfit().  The drawback is that
% % nlinfit() finds the best value of D and w for the sampled linear 
% % profiles, and not necessarily the entire particle.  So this requires
% % strategic sampling.
% % 2)OR you can write your own nlinfit() function for mulitple dimensions.
% % This will be very hard but may be more accurate because it uses a 
% % larger sampling space. 
% % Note: you've already written minimizeDW().m which is essentially a
% % multidimenisional nlinfit().m but only for D and w!

% % DO OPTION 1. Just write a subfuntion for collecting linear profile
% % on a sphere, which is what you did yesterday!

% % Routine should be:
% % 1) identify the sphereical region occupied by a particle using
% %    hist18jun2014.m.  Keep this region as is and zero the background.
% %    This way neighboring particles do NOT retard linear profile
% %    measurements of D and w.
% % 2) Obtain axial and lateral linear profiles on the particle from step
% %    1. See office-work contained in e-mail from 17-June-2014
% % 3) Do nonlinear regression analysis on these linear profiles to
% %    estimate D and w.  Compare estimates from each linear profile...
% %    they should all be nearly the same so you can safe take an average
% %    and assign D_avg and w_avg to the particle 

% % next is sub-pixel resolution analysis using the entire 3Dimensional
% % particle image and minize_error.m. I am not sure if D and w should be
% % found before or after sub-pixel analysis.  Try all 3! 
% % ...estimate D and w before sub-pixel analysis, after sub-pixel 
% % analysis, and both before and after.
