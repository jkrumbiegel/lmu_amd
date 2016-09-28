function coords = getRandomPointsInRect(n, x, y, safety)

safety = ceil(safety);

leftSafe = safety;
topSafe = safety;
rightSafe = x-safety;
bottomSafe = y-safety;

iterations = 1000000; %could be a more intelligent solution than brute forcing but I don't feel like it

%create first random point
points(1,:)=[randi([leftSafe,rightSafe]),randi([topSafe,bottomSafe])];

for i=2:n
    success = false;
    for tries = 1:iterations
        %create preliminary random point
        curPoint = [randi([leftSafe,rightSafe]),randi([topSafe,bottomSafe])];
        isSafe = true;
        %loop through all previously generated points
        for c=1:size(points,1)
            
            dist = getDistance(points(c,:),curPoint);
            
            %check new point for safety, if it fails try to generate another
            %point, otherwise continue checking
            if dist >= safety
                continue
            else
                isSafe = false;
                break
            end
        end
        %if the point passes all safety tests, add it to the list,
        %otherwise do another round
        if isSafe
            points(i,:) = curPoint;
            success = true;
            break
        else
            continue
        end
    end
    
    if success
        continue
    else
        warning('No safe random point could be found in ',num2str(iterations),' iterations. Point list is empty.');
        break
    end
end

if success
    coords = points;
else
    coords = [];
end

    function d = getDistance(a,b)
        x1 = a(1,1);
        x2 = b(1,1);
        y1 = a(1,2);
        y2 = b(1,2);
        d = sqrt((x2-x1)^2+(y2-y1)^2);
    end
end


