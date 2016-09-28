function coords = getLinSpacedPointsOnRectGrid(rect,xdist,ydist)

x = rect(1);
y = rect(2);
n = x*y;

coords = NaN(2,n);

counter = 1;
for i = 0 : x-1
    for j = 0 : y-1
        coords(1, counter) = i;
        coords(2, counter) = j;
        counter = counter+1;
    end
end

coords = coords + repmat([ -(x-1)/2; -(y-1)/2],[1,n]);
coords(1,:) = coords(1,:) .* xdist;
coords(2,:) = coords(2,:) .* ydist;
