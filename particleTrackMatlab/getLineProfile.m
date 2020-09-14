function [profile, diagPlane] = getLineProfile(center1, center2, imageData)
%Author: Eru K.
% adopted from get_linear_profile.m on 19-June-2014

% objective: To get the line profile between center1 and center2. Note: the
% linear profile, i.e. 'profile', is returned with pixel values in order
% from center1 to center2

% diagPlane = diagonal plane of REAL pixel values between center1 and
% center2 
% 
% note: getDiag().m returns the linear profile between centers1&2
% regardless if 'diagPlane' is a square matrix or not (if 'diagPlane' is
% not square then it uses linear interpolation) BUT 'diagPlane' returne by
% getLineProfile.m is NOT the linearly interpolated plane. It is the real
% pixel values

%         ______  
%        /*    /|
%      1/__ *_/ |     
%       |.    |*|
%       |  .  | /2
%       |____.|/     
%             
%   diagPlane =
% 
%       1______
%       .a1    *
%       .  a2  * 
%       .____a3* 
%              2 
% 
% profile = [a1 a2 a3];

r1 = uint16(center1(1));
r2 = uint16(center2(1));

c1 = uint16(center1(2));
c2 = uint16(center2(2));

s1 = uint16(center1(3));
s2 = uint16(center2(3));

% --------------------------(correction)---------------------------
if(r1 <= 0)
    r1 = 1;
end

if(r1 > size(imageData,1))
    r1 = size(imageData,1);
end

if(r2 <= 0)
    r2 = 1;
end

if(r2 > size(imageData,1))
    r2 = size(imageData,1);
end

if(c1 <= 0)
    c1 = 1;
end

if(c1 > size(imageData,2))
    c1 = size(imageData,2);
end

if(c2 <= 0)
    c2 = 1;
end

if(c2 > size(imageData,2))
    c2 = size(imageData,2);
end

if(s1 <= 0)
    s1 = 1;
end

if(s1 > size(imageData,3))
    s1 = size(imageData,3);
end

if(s2 <= 0)
    s2 = 1;
end

if(s2 > size(imageData,3))
    s2 = size(imageData,3);
end

% --------------------------(correction)---------------------------

if(s2 == s1)
    if(r1 == r2)
        if(c1 == c2)
            %A config (checked 23-June-2014, 24Jun14)
            profile = imageData(r1,c1,s1); %a point 
            diagPlane = profile;
        end
        
        if(c1 < c2) 
            cStart = c1;
            cStop = c2;
            % B config (checked 23-June-2014, 24Jun14)
            % 1-----2             
            profile = imageData(r1,cStart:cStop,s1); %a column vector 
            diagPlane = profile;
        end
        
        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            % C config(checked 23-June-2014, 24Jun14)
            % 2----1 ==> flipdim( ,2) = 1---2            
            profile = flipdim(imageData(r1,cStart:cStop,s1),2); %a column vector
            diagPlane = profile;
        end
    end
    
    if(r1 < r2) 
        rStart = r1;
        rStop = r2;
        if(c1 == c2)
            % D config (checked 23-June-2014, 24Jun14)
            %
            % 1
            % |
            % | 
            % |
            % 2
            diagPlane = imageData(rStart:rStop,c1,s1); % a row vector
            profile = permute(diagPlane,[2 1]); %aka profile = diagPlane';
        end
        
        if(c1 < c2)
            cStart = c1;
            cStop = c2;
            % E config (checked 23-June-2014, 24Jun14)
            %
            % 1_____
            % |     |
            % |     |
            % |_____|2
            % 
            diagPlane = imageData(rStart:rStop,cStart:cStop,s1);
            profile = getDiag(diagPlane); %getDiag() == diag_eru()
        end
        
        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            % F config (checked 23-June-2014, 24Jun2014)
            %
            %  _____1                      1_____  
            % |     |                      |     |
            % |     |  ===> rot90( ,1) ==> |     | (checked, 20Jun2014)
            % |_____|                      |_____|    
            % 2                                  2
            %
            % note: rot90(,+) rotates counter-clockwise
            %       rot90(,-) rotates clockwise
            diagPlane = rot90(imageData(rStart:rStop,cStart:cStop,s1),-1);
            profile = getDiag(diagPlane);
        end
    end
    
    if(r1 > r2)
        rStart = r2;
        rStop = r1;
        if(c1 == c2)
            % G config (checked 23-June-2014, 24-Jun-2014)
            %
            % 2                       1
            % |                       |
            % | ==> flipdim( , 1) ==> |
            % |                       |
            % 1                       2
            %
            profile = flipdim(imageData(rStart:rStop,c1,s1),1);
            diagPlane = profile;     
        end
        
        if(c1 < c2) 
            cStart = c1;
            cStop = c2;
            %H config (checked 23-June-2014, 24-June-2014)
            %
            %  ____2                     1_____    
            % |    |                     |     |
            % |    | ==> rot90( ,-1) ==> |     | (checked,20June2014)
            % |____|                     |_____|
            % 1                                2
            %
            diagPlane = rot90(imageData(rStart:rStop,cStart:cStop,s1),1);
            profile = getDiag(diagPlane);
        end
        
        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            % I config (checked 23-June-2014, 24-Jun-2014)
            %
            % 2____                     1____
            % |    |                    |    |
            % |    | ==> rot90( ,2) ==> |    |
            % |____|                    |____|
            %       1                         2
            %
            diagPlane = rot90(imageData(rStart:rStop,cStart:cStop,s1),2);
            profile = getDiag(diagPlane);
        end
    end
end

if(s1 < s2)
    sStart = s1;
    sStop = s2;
    if(r1 == r2)
        if(c1 == c2)
            %J config (checked 23-June-2014, 24Jun2014)
            %
            %    2
            %   /
            %  /   ==> permute(,[2 1]) ==> 1-----2 (checked, 20Jun2014)
            % /
            % 1
            diagPlane = permute(imageData(r1,c1,sStart:sStop),[1 3 2]);
            profile = diagPlane;            
        end

        if(c1 < c2) 
            cStart = c1;
            cStop = c2;
            %K config (checked 23-June-2014, 24-Jun2014 on both square and
            %non-square matrices)
            %   y______2
            %   /     /|                      
            % 1/_____/x|                     1 _____ y
            %  |     | |  ==> squeeze( ) ==>  |     |
            %  |     | /                      |     | (checked, 20Jun14)
            %  |_____|/                       |_____|
            %                                x       2
            %
            diagPlane = squeeze(imageData(r1,cStart:cStop,sStart:sStop));
            profile = getDiag(diagPlane);
        end

        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            %L config (checked 23-June-2014, 24Jun2014)
            %   2______y                              _____1
            %  x/____ /|                             |     |                   
            %  |     1 |  ==> permute( ,[3 2 1]) ==> |     |  
            %  |     | |                             |_____|                    
            %  |_____|/                             2       
            %                 
            %note: permute(,[3 2 1]) takes the 3rd dimension and places
            %it in the 1st. Does nothing to the 2nd dimension. And
            %moves the 1st dimension to the 3rd
            %
            %                    1_____ y
            %                    |     |
            % ==> rot90(, 1) ==> |     |  (checked 20jun2014)
            %                    |_____|
            %                   x       2
            diagPlane = rot90(permute(imageData(r1,cStart:cStop,sStart:sStop),[3 2 1]),1);
            profile = getDiag(diagPlane);
        end
    end

    if(r1 < r2)
        rStart = r1;
        rStop = r2;
        if(c1 == c2)
            %M config (checked 23-June-2014, 24Jun2014)
            %    ______                              1_____
            %   /:    /|                             |     |
            % 1/_:___/ |  ==> permute( ,[1 3 2]) ==> |     | (checked 20Jun14)
            %  | :   | |     (switch 3rd and 2nd     |_____|
            %  | /2  | /        dimensions)                 2
            %  |/____|/
            %
            diagPlane = permute(imageData(rStart:rStop,c1,sStart:sStop),[1 3 2]);
            profile = getDiag(diagPlane);
        end

        if(c1 < c2)
            cStart = c1;
            cStop = c2;
            %N config (checked 23-June-2014, 24-Jun-2014)
            %    ______
            %   /   . /|                          1_____
            % 1/._'__/ |  ===> getDiagPlane() ==> |     |
            %  |     |.|2                         |     |
            %  |   . | /                          |_____|
            %  |.____|/                                  2
            %
            diagPlane = getDiagPlane(imageData(rStart:rStop,cStart:cStop,sStart:sStop));
            profile = getDiag(diagPlane);
        end

        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            %O config (checked 24-June-2014, square & non-square matrix)
            %     _______                          ______
            %    /:     /|                        /      /|
            %   /_:_"_./1|  ==> flipdim(,2) ==> 1/______/ |    
            %   |2:    | |                       |      | |2
            %   | / .  | /                       |      | /
            %   |/____.|/                        |______|/
            %                        1 _____ 
            %                         |     |
            %  ==> getDiagPlane() ==> |     |
            %                         |_____|
            %                               2
            %
            diagPlane = getDiagPlane(flipdim(imageData(rStart:rStop,cStart:cStop,sStart:sStop),2));
            profile = getDiag(diagPlane);
        end
    end

    if(r1 > r2)
        rStart = r2;
        rStop = r1;
        if(c1 == c2)
            %P config (checked 24-June-2014, square & non-square matrix)
            %   2______
            %   /:    /|                             _____2
            %  /_:___/ | ==> permute(,[1 3 2]) ==>  |     |
            %  | :   | |                            |     |
            %  | :   | /                            |_____|
            %  |/____|/                             1
            %  1
            %                     1______
            %                     |      |
            % ==> flipdim(,1) ==> |      |  (checked, 20Jun2014)
            %                     |______|
            %                            2
            %
            diagPlane = flipdim(permute(imageData(rStart:rStop,c1,sStart:sStop),[1 3 2]),1);
            profile = getDiag(diagPlane);
        end

        if(c1 < c2)
            cStart = c1;
            cStop = c2;
            %Q config (checked 24-June-2014, square & non-square matrices)
            %                                   ______
            %     ______2                      /     /| 
            %    /     /|                    1/_____/ |
            %   /_____/ | ==> flipdim(,1) ==> |     | |
            %   |     | |                     |     | /2
            %   |     | /                     |_____|/
            %   |_____|/
            %   1
            %                          1_____
            %                          |     |
            %  ==> getDiagPlane() ==>  |     |
            %                          |_____|2
            %
            diagPlane = getDiagPlane(flipdim(imageData(rStart:rStop,cStart:cStop,sStart:sStop),1));
            profile = getDiag(diagPlane);
        end

        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            %R config (checked 24-June-2014, square & non-square matrices)
            %   2______                              1 _____
            %   /     /|                             /     /|
            %  /_____/ |  ==> permute(,[2 1 3]) ==> /_____/ |
            %  |     | |      (switch columns       |     | |
            %  |     | /        and rows)           |     | /
            %  |_____|/                             |_____|/
            %         1                                     2
            %                        _____                      
            %                       /     /|
            %                     1/_____/ |
            %                      |     | |2 
            %  ==> flipdim(,3) ==> |     | /
            %                      |_____|/
            %
            %
            %                       1 _____
            %                        |     |
            %  ==> getDiagPlane ==>  |     |
            %                        |_____|
            %                               2
            %
            diagPlane = getDiagPlane(flipdim(permute(imageData(rStart:rStop,cStart:cStop,sStart:sStop),[2 1 3]),3));
            profile = getDiag(diagPlane);
        end
    end
end

if(s1 > s2)
    sStart = s2;
    sStop = s1;
    if(r1 == r2)
        if(c1 == c2)
            %S config (checked 24-June-2014)
            %    1_____                        2_____
            %   /     /|                       /     /|
            %  /_____/ | ==> flipdim(,3) ==> 1/_____/ |
            % 2|     | |                      |     | |
            %  |     | /                      |     | /
            %  |_____|/                       |_____|/
            %
            %
            %
            %
            %  ==> squeeze() ==> 1----2  (checked 20June2014)
            %
            %
            diagPlane = squeeze(flipdim(imageData(r1,c1,sStart:sStop),3));
            profile = diagPlane;
        end
        
        if(c1 < c2)
            cStart = c1;
            cStop = c2;
            %T config (checked 24-June-2014)
            %   1______                         x______2   
            %   /     /|                        /     /|
            %  /_____/ |                      1/_____/ |
            %  |     2 |  ==> flipdim(,3) ==>  |     y |
            %  |     | /                       |     | /
            %  |_____|/                        |_____|/
            %
            %
            %                    1_____x   
            %  ==> squeeze() ==> |     |
            %                    |     |    (checked 20Jun2014)
            %                    y_____|
            %                           2
            %
            diagPlane = squeeze(flipdim(imageData(r1,cStart:cStop,sStart:sStop),3));
            profile = getDiag(diagPlane);
        end
        
        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            %U config (checked 24-June-2014)
            %  x______1                      
            %  /_____/|                     2 _____x 
            % 2|     y|  ==> squeeze( ) ==>  |     |
            %  |     ||                      |     | 
            %  |_____|/                      |_____|
            %                               y       1
            %
            %                    1 ____ y
            %                    |     |
            %  ==> rot90(,2) ==> |     |  (checked 20Jun2014)
            %                    |_____|
            %                    x     2     
            %
            diagPlane = rot90(squeeze(imageData(r1,cStart:cStop,sStart:sStop)),2);
            profile = getDiag(diagPlane);
        end
    end
    
    if(r1 < r2)
        rStart = r1;
        rStop = r2;
        if(c1 == c2)
            %V config (opposite of P config) (checked 24-Jun-2014)
            %   1______
            %   /:    /|                             _____1
            %  /_:___/ | ==> permute(,[1 3 2]) ==>  |     |
            %  | :   | |                            |     |
            %  | :   | /                            |_____|
            %  |/____|/                             2
            %  2
            %
            %                     1______
            %                     |      |
            % ==> flipdim(,2) ==> |      |  
            %                     |______|
            %                            2
            %
            diagPlane = flipdim(permute(imageData(rStart:rStop,c1,sStart:sStop),[1 3 2]),2);
            profile = getDiag(diagPlane);
        end
        
        if(c1 < c2)
            cStart = c1;
            cStop = c2;
            %X config (checked 24-Jun-2014, square and non-square matrices)
            %   1______                           _____
            %   /     /|                         /     /|
            %  /_____/ |  ==> flipdim(,3)==>   1/_____/ |
            %  |     | |                        |     | |2 
            %  |     | /                        |     | /
            %  |_____|/                         |_____|/
            %         2                                   
            %                       1 _____   
            %                        |     |
            % ==> getDiagPlane() ==> |     |
            %                        |_____|
            %                               2
            %
            diagPlane = getDiagPlane(flipdim(imageData(rStart:rStop,cStart:cStop,sStart:sStop),3));
            profile = getDiag(diagPlane);
        end
        
        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            %Y config (opposite of Q) (checked 24Jun14, square &
            %non-square matrices)
            %
            %                                  1______
            %     ______1                      /     /| 
            %    /     /|                     /_____/ |
            %   /_____/ | ==> flipdim(,2) ==> |     | |
            %   |     | |                     |     | /
            %   |     | /                     |_____|/
            %   |_____|/                            2
            %   2
            %                           ______
            %                          /     /|
            %                        1/_____/ |  (checked 20Jun2014)
            %                         |     | |
            %   ===> flipdim(,3) ==>  |     | /2
            %                         |_____|/
            %
            %
            %                         1 _____  
            %  ==> getDiagPlane() ==>  |     |
            %                          |     |
            %                          |_____|
            %                                2
            %
            diagPlane = getDiagPlane(flipdim(flipdim(imageData(rStart:rStop,cStart:cStop,sStart:sStop),2),3));
            profile = getDiag(diagPlane);
        end
    end
    
    if(r1 > r2) 
        rStart = r2;
        rStop = r1;
        if(c1 == c2)
            %Z config (checked 24-Jun-14, square and non-square matrices)
            %     ______ y                             ______
            %    /     /|                             /     /|
            %   /_____/ | ==> permute(,[1 3 2]) ==> 2/_____ y|
            %   |    2| |     (keep rows the same    |     | |
            %   |     | /1     switch columns with   |     | /
            %   |_____|/       slices)               |_____|/
            %         x                             x       1 
            %                           
            %                    1 _____ x
            %  ==> rot90(,2) ==>  |     |    (Checked 20Jun2014)
            %                     |     | 
            %                     |_____|
            %                    y       2
            %           
            diagPlane = rot90(permute(imageData(rStart:rStop,c1,sStart:sStop),[1 3 2]),2);
            profile = getDiag(diagPlane);
        end
        
        if(c1 < c2)
            cStart = c1;
            cStop = c2;
            %AA config (checked 23Jun2014) (checked 24Jun14, square and
            %non-square matrices)
            %    ______                           ______ 2
            %   /:    /|                         /     /|
            %  /_:___2 |  ==> flipdim(,3)==>    /_____/ |
            %  | :   | |                        |     | | 
            %  | /1  | /                        |     | /
            %  |/____|/                         |_____|/
            %                                  1
            %
            %                          _____
            %                         /     /|
            %  ==> flipdim(,1) ==>  1/_____/ |
            %                        |     | |
            %                        |     | /2
            %                        |_____|/
            %       
            %
            %
            %                         1_____
            %                         |     |
            % ==> getDiagPlane() ==>  |     |
            %                         |_____|
            %                               2
            %
            diagPlane = getDiagPlane(flipdim(flipdim(imageData(rStart:rStop,cStart:cStop,sStart:sStop),3),1));
            profile = getDiag(diagPlane);
        end
        
        if(c1 > c2)
            cStart = c2;
            cStop = c1;
            %BB config (checked 24Jun2014, square & non-square matrices)
            %     _____                          _____
            %    /     /|                       /:    /|
            %  2/_____/ |  ==> flipdim(,2) ==> /_:___/2|
            %   |     | | 1                    | :   | |
            %   |     | /                      | /1  | /
            %   |_____|/                       |/____|/
            %
            %                       1 _____
            %                       /     /|
            %  ==> flipdim(,1) ==> /_____/ |
            %                      |     | |
            %                      |     | /
            %                      |_____|/
            %                             2
            %                                                       
            %                       x_____                             
            %                      /     /|                       1 _____ y 
            % ==>flipdim(,3) ==> 1/_____/ |  ==> getDiagPlane()==> |     | 
            %                     |     | |                        |     | 
            %                     |     | /2                       |_____|   
            %                     |_____|/                        x      2
            %                           y
            %
            diagPlane = getDiagPlane(flipdim(flipdim(flipdim(imageData(rStart:rStop,cStart:cStop,sStart:sStop),2),1),3));
            profile = getDiag(diagPlane);
        end
    end
end
end