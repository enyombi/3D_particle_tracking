function [profile, diagPlaneNew] = getDiag(diagPlane)
% Author: Eru K.
% Objective: return the diagonal elements of 'diagPlane'.  If 'diagPlane'
% is not a square matrix than use linear interpolation to transform it into
% a square matrix
% 
% Note: strange cases like finding the diagonal of a 1D array...
% 
%    _____
%   /    /|
% 1/____/ |
%  |    | /
%  |____|/
% 2 
% 
%   2_____
%   /    /|
% 1/____/ |
%  |    | /
%  |____|/
%
%    _____
%   /    /|
% 1/____/2|
%  |    | /
%  |____|/
%  
%    _____
%   /    /|
% 1/2___/ |
%  |    | /
%  |____|/
% 
% 
% ...are dealt with in getLineProfile.m and therefore should NEVER reach
% getDiag.m


rows = size(diagPlane,1);
columns = size(diagPlane,2);

[yy xx] = ndgrid(1:rows,1:columns);

if(rows < columns)
    %add more rows...
    %                       1___________
    %  1___________         |           |
    %  |           |  ==>   |           |
    %  |___________|        |           |
    %              2        |           |
    %                       |___________|
    %                                   2
    %
    [yyi xxi] = ndgrid(linspace(1,rows,columns),1:columns);
    diagPlaneNew = interp2(xx,yy,diagPlane,xxi,yyi);
    profile = diag(diagPlaneNew); 
end

if(rows > columns) 
    %add more columns...
    %
    %    ___         __________
    %   |   |       |          |
    %   |   | ==>   |          |
    %   |   |       |          | 
    %   |___|       |__________|
    %
    [yyi xxi] = ndgrid(1:rows,linspace(1,columns,rows));
    diagPlaneNew = interp2(xx,yy,diagPlane,xxi,yyi);
    profile = diag(diagPlaneNew); 
end

if(rows == columns) %regular square matrix
    diagPlaneNew = diagPlane;
    profile = diag(diagPlaneNew);
end
end