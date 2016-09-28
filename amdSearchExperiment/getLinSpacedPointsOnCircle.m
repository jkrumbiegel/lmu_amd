function coords = getLinSpacedPointsOnCircle(n, radius, degOffset)
% Angles at which our stimuli's centers will be. We start at degOffset
% and then equally space points around the edge of a circle
anglesDeg = linspace(0+degOffset, 360+degOffset, n+1);
%Delete last (duplicate) point
anglesDeg(end)=[];
%Convert to Radians
anglesRad = anglesDeg * (pi / 180);

% X and Y coordinates of the points on the circle
yPosVector = sin(anglesRad) .* radius;
xPosVector = cos(anglesRad) .* radius;

coords = [xPosVector; yPosVector];