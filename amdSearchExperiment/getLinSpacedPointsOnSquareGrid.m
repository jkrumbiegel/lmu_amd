function coords = getLinSpacedPointsOnSquareGrid(n, distance)

if mod(sqrt(n),1) == 0
    x = sqrt(n);
    y = x;
else
    error([num2str(n),' elements not possible. Only perfect squares.']);
end

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
coords = coords .* distance;

