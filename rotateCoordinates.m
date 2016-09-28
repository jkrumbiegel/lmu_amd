function rotCoords = rotateCoordinates(coords,degAngle)

rotMatrix = [cosd(degAngle), -sind(degAngle);...
             sind(degAngle),  cosd(degAngle)];
         
rotCoords = rotMatrix * coords;