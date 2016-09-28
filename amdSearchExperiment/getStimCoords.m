function coords = getStimCoords(stimulus, Scale, rotation, correctionPx)

coords = [];

if nargin == 1
    Scale = 1;
end

switch stimulus
    case 'E'
        %a capital 'E', standard size is width of 1 (unitless)
        sq2 = sqrt(2);
        xCoords = Scale .* [-0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 -0.5];
        yCoords = Scale .* [sq2/2 sq2/2 0 0 -sq2/2 -sq2/2 sq2/2 -sq2/2];
        coords = [xCoords; yCoords];
        if nargin == 4
            cb = correctionPx/2;
            coords = coords + [cb 0 cb 0 cb 0 0   0;
                               0  0 0  0 0  0 cb -cb];
        end
        if nargin == 3
            coords = rotateCoordinates(coords,rotation);
        end
    case 'F'
        %a capital 'F', standard size is width of 1 (unitless)
        sq2 = sqrt(2);
        xCoords = Scale .* [-0.5 0.5 -0.5 0.5 -0.5 -0.5];
        yCoords = Scale .* [-sq2/2 -sq2/2 0 0 sq2/2 -sq2/2];
        coords = [xCoords; yCoords];
        if nargin == 4
            cb = correctionPx/2;
            coords = coords + [cb 0 cb 0 0   0;
                               0  0 0  0 cb -cb];
        end
        if nargin == 3
            coords = rotateCoordinates(coords,rotation);
        end
    case 'X'
        sq2 = sqrt(2);
        xCoords = Scale .* [-sq2/2 sq2/2 -sq2/2 sq2/2];
        yCoords = Scale .* [-sq2/2 sq2/2 sq2/2 -sq2/2];
        coords = [xCoords; yCoords];
        if nargin == 3
            coords = rotateCoordinates(coords,rotation);
        end
    case 'O'
        %a circle with a standard diameter of 1
        left = Scale * (-1);
        top = Scale * (-1);
        right = Scale * 1;
        bottom = Scale * 1;
        coords = [left; top; right; bottom];
    case '/'
        %a vertical line with a standard length of 1
        xCoords = Scale .* [0 0];
        yCoords = Scale .* [-0.5 0.5];
        coords = [xCoords; yCoords];
        if nargin == 3
            coords = rotateCoordinates(coords,rotation);
        end
end
