function coords = getMMTfixPoly(n, rotation, widthDeg, innerCircle, outerCircle)
% coords = GETMMTFIXPOLY(n, rotation, widthDeg, innerCircle, outerCircle)
%
% Get coordinates of polygons forming a macular mapping test star pattern
% with 'n' spikes, 'rotation' rotation in degrees, 'widthDeg' width of each
% spike in degrees, 'innerCircle' inner radius on which the inner spike
% tips lie and 'outerCircle' outer radius on which the outer spike edges
% lie.
%
% The coordinates are meant to be drawn with Psychtoolbox's
% Screen('FillPoly') function.

innerPoints = getLinSpacedPointsOnCircle(n,innerCircle,rotation);

outerPoints = getLinSpacedPointsOnCircle(n, outerCircle, rotation);
outerPointsCW = rotateCoordinates(outerPoints,widthDeg/2);
outerPointsCCW = rotateCoordinates(outerPoints,-widthDeg/2);

coords = [reshape(innerPoints,    [1,n*2]);...
          reshape(outerPointsCCW, [1,n*2]);...
          reshape(outerPointsCW,  [1,n*2])];