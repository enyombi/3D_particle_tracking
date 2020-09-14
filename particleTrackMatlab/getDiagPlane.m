function diagPlane = getDiagPlane(cutout)

% Author: Eru K.
% Date: 4-Dec-2013

% Objective: To obtain the diagonal plane betweeen center1 and center2 of
% cutout
% 
%   cutout = 
%         ______  
%        /*    /|
%      1/__ *_/ |     
%       |.    |*|
%       |  .  | /2
%       |____.|/     
%             
%   diagPlane =
% 
%       1_____
%       .     *
%       .     * 
%       ._____* 
%              2 
%


% simplification of getDiagPlane_23Jun2014.m
% This is a simplification because you do NOT need to consider strange 
% cases, such as... (these cases have been accounted for in
% getLineProfile.m)
% 
%cutout is only 2-dimensional:
%(e.g., size(diagPlaneTemp) = [29 1 7] or [1 29 7] or [29 7 1])
%   Case1:      
%        .     .
%      1._____.      
%       |     | 
%       | .   | .
%       |.____|.     
%              2
%   Case2:
%         ......
%        /|    
%      1/.|....     
%       | |2....
%       | /   
%       |/.....
% 
%   Case3:
%         _______2               
%        /:     /:   
%      1/_:____/ :         
%       : :    :   
%       :      : 
% 
%cutout is only 1-dimensional:
%(e.g., size(diagPlaneTemp) = [29 1 1] or [1 29 1] or [1 1 29])
% 
% Case1:
% size(cutout) = [1 cols 1], where cols >= 1
%        .      .
%      1._____2.      
%       '     ' 
%       ' .   ' .
%       ._____'. 
% 
% Case2:
% size(cutout) = [rows 1 1], where rows >= 1
%         .    
%        .     .
%      1 -----.      
%       |     ' 
%       | .   ' .
%       2.----'.
% 
% Case3:
% size(cutout) = [1 1 slice], where slice >= 1
%         2    
%        /     .
%      1/-----.      
%       '     ' 
%       ' .   ' .
%       '.----'. 
% 

diagPlaneTemp = zeros([size(cutout,3), length(getDiag(cutout(:,:,1)))],'single');

for z = 1:size(cutout,3)
    diagPlaneTemp(z,:) = getDiag(cutout(:,:,z));
end
diagPlane = diagPlaneTemp;  
end