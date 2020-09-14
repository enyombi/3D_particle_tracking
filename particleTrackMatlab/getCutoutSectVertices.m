function [corners, peakNew] = getCutoutSectVertices(peak,dist,dimImage)
% Author: Eru K.
% Date: 2-July-2014

% corners = vertices of a rectangular region surrounding 'peak' in a
%            3D-image
% dist = distance from 'peak' that each face of the rectangular region
%        cutout from the 3D-image should be  
% dim = dimensions of the 3D-image

% objective: returns the corners/vertices of a rectangular section of a 
% larger 3D-image to the nearest whole pixel.  This rectangular region 
% is centered on 'peak' and has dimesions (dist)X(dist)X(dist) 

% 
%       /                                        /
%      /                                        /
%     /                                        /
%    /                                        /
%   /________________________________________/
%   |                                        |
%   |                                        |
%   |                                        |
%   |               2________3               |
%   |               /.       /|              |
%   |             1/_.______/ |              |
%   |              | .      4 |              |
%   |              | .  x   | |              |
%   |              | 6......|./7             |
%   |              |._______|/               |
%   |             5         8                |    /
%   |                                        |   /
%   |                                        |  /
%   |                                        | /
%   |________________________________________|/
% 
% x = peak's position at center of box

%round dist and peak to the nearest whole pixel value to get vertices that
%are whole pixel NOT fractions of pixels
% dist = round(dist); 
% peak = round(peak); 

% vertices of cutoutSect
v(1,:) = peak - dist; %topFrontCorner
v(7,:) = peak + dist; %bottomBackCorner

v(2,:) = [peak(1)-dist, peak(2)-dist, peak(3)+dist];
v(3,:) = [peak(1)-dist, peak(2)+dist, peak(3)+dist];
v(4,:) = [peak(1)-dist, peak(2)+dist, peak(3)-dist];

v(5,:) = [peak(1)+dist, peak(2)-dist, peak(3)-dist];
v(6,:) = [peak(1)+dist, peak(2)-dist, peak(3)+dist];
v(8,:) = [peak(1)+dist, peak(2)+dist, peak(3)-dist];

v = round(v);
% now make certain that the rectangular box around the peak does not
% exceed the bounds of the 3D image...
for np = 1:length(v) %particle #
    for dim = 1:3 %each particle has 3 dimensions/coords
        if(v(np,dim) <= 0)
            v(np,dim) = 1;
        end

        if(v(np,dim) > dimImage(dim))
            v(np,dim) = dimImage(dim);
        end
    end
end

% Now find the position of 'peak' inside the rectangular section
%       ________  
%      /.       /|
%[1,1,1].______/ |
%     | .      | |
%     | .  x   | |
%     | .......|./
%     |._______|/
%             
displacement = peak-v(1,:);
peakNew = [1,1,1] + displacement; %the position of 'peak' inside the rectangular cutout section

% % % % This correction should NOT be necessary!
% % % for dim = 1:3 %each particle has 3 dimensions/coords
% % %     if(peakNew(dim) <= 0)
% % %         peakNew(dim) = 1;
% % %     end
% % % 
% % %     if(peakNew(dim) > dimImage(dim))
% % %         peakNew = dimImage(dim);
% % %     end
% % % end

corners = v;
end