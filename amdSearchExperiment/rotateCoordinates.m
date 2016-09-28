function rotCoords = rotateCoordinates(coords,degAngle)
% ROTATECOORDINATES rotates coordinates by an angle given in degrees.
% The pivot point is at [0,0]. If you want to rotate around another point
% translate the coordinates by the appropriate amount first, rotate them,
% and translate them back.

rotMatrix = [cosd(degAngle), -sind(degAngle);...
             sind(degAngle),  cosd(degAngle)];
         
rotCoords = rotMatrix * coords;