function coords = getMMTfixPoly(n, rotation, widthDeg, innerCircle, outerCircle)

innerPoints = getLinSpacedPointsOnCircle(n,innerCircle,rotation);

outerPoints = getLinSpacedPointsOnCircle(n, outerCircle, rotation);
outerPointsCW = rotateCoordinates(outerPoints,widthDeg/2);
outerPointsCCW = rotateCoordinates(outerPoints,-widthDeg/2);

coords = [reshape(innerPoints,    [1,n*2]);...
          reshape(outerPointsCCW, [1,n*2]);...
          reshape(outerPointsCW,  [1,n*2])];