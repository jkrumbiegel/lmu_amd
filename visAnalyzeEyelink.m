function out = visAnalyzeEyelink(filename)

if nargin == 0
    [filename,pathname,~] = uigetfile();
end
filename = strsplit(filename,'.');
filename = filename{1};
fnparts = strsplit(filename,'-');
edfname = strcat(fnparts{1},fnparts{5},fnparts{6},fnparts{7});

expData = load([pathname,filename,'.mat']);
expData = expData.exp;
elData = load([pathname,edfname,'.mat']);
elData = elData.elData;

tMax = length(elData);
cur = 1;

showGaze = true;
dataSelect = 1; % 0 all, 1 only search, 2 only fixation
showFixation = false;
showTouch = false;
showSaccades = false;
playingAnimation = false;

curGaze = [];
curFix = [];
curSac = [];

gazePlotSize = 25;

textCounter = 0;

f = figure;
warning('off', 'Images:initSize:adjustingMag');
set(f, 'KeyPressFcn', @(h_obj, evt) KeyPress(h_obj, evt));
drawGraph();
    %%
    function KeyPress(h_obj, evt)
        changed = false;
        if(strcmp(evt.Key, 'rightarrow'))
            cur = cur+1;
            if cur > tMax
                cur = 1;
            end
            changed = true;
        end
        
        if(strcmp(evt.Key, 'leftarrow'))
            cur = cur-1;
            if cur < 1
                cur = tMax;
            end
            changed = true;
        end
        
        if(strcmp(evt.Key, 'g'))
            showGaze = ~showGaze;
            changed = true;
        end
        
        if(strcmp(evt.Key, 'f'))
            showFixation = ~showFixation;
            changed = true;
        end
        
        if(strcmp(evt.Key, 't'))
            showTouch = ~showTouch;
            changed = true;
        end
        
        if(strcmp(evt.Key, 'd'))
            dataSelect = mod(dataSelect+1,3);
            changed = true;
        end
        
        if(strcmp(evt.Key, 's'))
            showSaccades = ~showSaccades;
            changed = true;
        end
        
        if(strcmp(evt.Key, 'space'))
            playAnimation();
        end
        
        if changed
            playingAnimation = false;
            drawGraph();
        end
    end
    %%
    function showScreenshot(i)
        imshow(expData.trial(i).Screenshot);
    end
    %%
    function drawGraph()        
        showScreenshot(cur);
        hold on
        if showGaze
            curGaze = getGaze();            
            plotGaze();
        end
        if showSaccades
            curSac = getSaccades();
            plotSaccades();
        end
        if showFixation
            curFix = getFixation();
            plotFixation();
        end
        if showTouch
            plotTouch();
        end
        
        drawInfoText();
        
        hold off
    end
    %%
    function g = getGaze()
        g = elData(cur).Data(:,2:3);
        if dataSelect > 0
            searchTime = getSyncTime();
            searchSample = find(elData(cur).Data(:,1) == searchTime);
            switch dataSelect
                case 1
                    g(1:searchSample-1,:)=[];
                case 2
                    g(searchSample:end,:)=[];
            end
        end
    end
    %%
    function f = getFixation()
        f = horzcat(elData(cur).fix.H,elData(cur).fix.V,elData(cur).fix.dur,elData(cur).fix.Tstart,elData(cur).fix.Tend);
        if dataSelect > 0
            searchTime = getSyncTime();
            switch dataSelect
                case 1
                    f(f(:,5)<searchTime,:)=[]; % show only fixations in search period
                case 2
                    f(f(:,4)>=searchTime,:)=[]; % show only fixations in fixation period
            end
        end
    end
    %%
    function drawInfoText()
        leftmargin = 10;
        vertdistance = 25;
        %height = size(expData.trial(cur).Screenshot,1);
        text(leftmargin,next(1)*vertdistance,['Subject: ',fnparts{1}]);
        switch elData(cur).fix.eye(1)
            case 0
                eyeText = 'Left';
            case 1
                eyeText = 'Right';
        end
        text(leftmargin,next()*vertdistance,['Eye: ',eyeText]);
        text(leftmargin,next()*vertdistance,['Trial: ',num2str(cur), ' of ',num2str(tMax)]);
        text(leftmargin,next()*vertdistance,['Search Time: ',num2str(elData(cur).Data(end,1)-getSyncTime()),' ms']);
        switch dataSelect
            case 0
                dataText = 'Data: Fixation and Search';
            case 1
                dataText = 'Data: Search only';
            case 2
                dataText = 'Data: Fixation only';
        end
        text(leftmargin,next()*vertdistance,dataText);
        switch showGaze
            case false
                gazeText = 'Gaze: hidden';
            case true
                gazeText = 'Gaze: visible';
        end
        text(leftmargin,next()*vertdistance,gazeText);
        switch showFixation
            case false
                fixText = 'Fixation: hidden';
            case true
                fixText = 'Fixation: visible';
        end
        text(leftmargin,next()*vertdistance,fixText);
        switch showTouch
            case false
                touchText = 'Touch: hidden';
            case true
                touchText = 'Touch: visible';
        end
        text(leftmargin,next()*vertdistance,touchText);
        switch showSaccades
            case false
                sacText = 'Saccades: hidden';
            case true
                sacText = 'Saccades: visible';
        end
        text(leftmargin,next()*vertdistance,sacText);
    end
    %%
    function s = getSaccades()
        s = horzcat(elData(cur).sac.Hstart,elData(cur).sac.Hend,...
            elData(cur).sac.Vstart,elData(cur).sac.Vend,...
            elData(cur).sac.dur,elData(cur).sac.amp,elData(cur).sac.pv,...
            elData(cur).sac.Tstart,elData(cur).sac.Tend);
        if dataSelect > 0
            searchTime = getSyncTime();
            switch dataSelect
                case 1
                    s(s(:,9)<searchTime,:)=[]; % show only saccades in search period
                case 2
                    s(s(:,8)>=searchTime,:)=[]; % show only fixations in fixation period
            end
        end
    end
    %%
    function t = getSyncTime()
        t = elData(cur).msg.time(find(not(cellfun('isempty', strfind(elData(cur).msg.text,'SYNCTIME'))))); % Very pretty
    end
    %%
    function plotGaze()        
        colormap parula;
        c = cumsum(getDistanceVector(curGaze));
        scatter(curGaze(:,1),curGaze(:,2),gazePlotSize,c,'filled');
    end

    %%
    function playAnimation()
        
        showScreenshot(cur);
        hold on
        
        drawInfoText();
        
        n = size(curGaze,1);
        lengthSecs = n/1000;
        c = cumsum(getDistanceVector(curGaze));
        colors = getColorVector(n,parula,c);        
        speed = 20; %fps
        waittime = 1/speed;
        %frames = ceil(n/speed);
        samplesPerFrame = (n/lengthSecs)/speed;
        frames = ceil(n/samplesPerFrame);
        
        playingAnimation = true;
        for frame = 1:frames
            if playingAnimation
                before = GetSecs();
                curBegin = 1+((frame-1)*samplesPerFrame);
                curEnd = frame * samplesPerFrame;
                if curEnd > n
                    curEnd = n;
                end
                curCol = colors(curBegin:curEnd,:);
                scatter(curGaze(curBegin:curEnd,1),curGaze(curBegin:curEnd,2),gazePlotSize,curCol,'filled');
                drawnow;
                after = GetSecs();
                %disp(after-before);
                toWait = waittime - (after - before);
                if toWait > 0
                    WaitSecs(toWait);
                end
            else
                break;
            end
        end
        playingAnimation = false;
        hold off
        
    end

    %%
    function plotFixation()
        n = size(curFix,1);
        colors = getColorVector(n,parula);
        scatter(curFix(:,1),curFix(:,2),curFix(:,3),colors,'linewidth',3);
    end

    %%
    function cv = getColorVector(n, cmap, samplepoints)
        csize = size(cmap,1);
        if nargin == 2
            points = linspace(1,csize,n);            
        end
        if nargin == 3
            high = samplepoints(end);
            low = samplepoints(1);            
            points = (((samplepoints-low)./high).*(csize-1))+1;
        end
        cv = interp1(1:csize,cmap,points);
            
    end

    %%
    function plotSaccades()
        n = size(curSac,1);
        map = parula;
        colors = interp1(1:size(map,1),map,linspace(1,size(map,1),n));
        for i=1:n
            xs = curSac(i,1);
            xe = curSac(i,2);
            ys = curSac(i,3);
            ye = curSac(i,4);
            quiver( xs,ys,xe-xs,ye-ys,0,'linewidth',3,'MaxHeadSize',8,'Color',colors(i,:));
        end
    end

    %%
    function plotTouch()
        xy = expData.trial(cur).touchPos;
        scatter(xy(1),xy(2),1000,'black','x','LineWidth',3);
    end
    %%
    function d = getDistanceVector(g)
        for i = length(g):-1:2
            d(i) = sqrt((g(i,1)-g(i-1,1))^2 + (g(i,2)-g(i-1,2))^2);
        end
        d(1) = 0;
    end
    %
    function n = next(init)
        if nargin == 1
            textCounter = init;
        else
            textCounter = textCounter+1;
        end
        n = textCounter;
    end
end