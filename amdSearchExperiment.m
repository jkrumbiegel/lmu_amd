% Clear the workspace and the screen
sca;
close all;
clearvars;

exp.motiontracking = false;
exp.eyetracking = false;

cd='C:\Users\ru35pec\Desktop\Julius\Macular Degeneration\AMD Experimente';

exp.subjectinitials = input('Enter participant code (2 letters): ','s');
if isempty(exp.subjectinitials)
    exp.subjectinitials = 'temp';
end
exp.timestamp = datestr(datetime('now','TimeZone','local'),'yyyy-mm-dd-HH-MM-ss');
filename = strcat(exp.subjectinitials,'-',exp.timestamp);
eyelinkfilename = strcat(exp.subjectinitials,datestr(datetime('now','TimeZone','local'),'HHMMss'))  ;

% put the current time into the random number generator as a seed value
rng(now);

% Access the arduino that sends 5v pulses upon completion to the
% MotionMonitor software
if exp.motiontracking
    arduino = serial('COM5'); % Check the correct COM port in device manager beforehand
    fopen(arduino);
    pause(3); % Allow the arduino a short moment to establish the connection
              % Don't start the motion tracking before, some voltage spikes
              % that occur when accessing the arduino the first time will
              % already trigger the recording
end


% Some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Select the external screen if there are two (or more)
screenNumber = max(screens);
HideCursor(screenNumber);

% Define some colors
white = [1;1;1];
black = [0;0;0];
grey = white / 2;
red = [1;0;0];
green = [0;1;0];
blue = [0;0;1];

% Define draw type codes. These are used to mark stimuli with the type of
% drawing routine they require, e.g., if they consist of lines they need
% line drawing, if they consist of center points and a radius they need a
% circle drawing method etc. This way all stimuli of the same type can be
% grouped and drawn in one command later.
drawTypeLine = 1;
drawTypeFrameOval = 2;

% Open an on screen window using PsychImaging and color it grey.
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey,[],[],[],[],4);
[xCenter, yCenter] = RectCenter(windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Measure the vertical refresh rate of the monitor
ifi = Screen('GetFlipInterval', window);

% Retreive the maximum priority number and set max priority
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

% Get screen resolution
resolution = Screen('Resolution', screenNumber);
% Enter the measured screen width
exp.screenWidthCM = 52.13; %iiyama prolite t2452mts
% Enter the subject's distance to the monitor in cm
exp.screenDistanceCM = 57;
% Calculate how many pixels are 1 degree of visual angle
pixPerCM = resolution.width/exp.screenWidthCM;
va = 2 * atand(0.5/exp.screenDistanceCM); % Visual Angle for 1 cm
pva = pixPerCM / va;
exp.pixPerVA = pva;

% Set trial parameters
exp.searchMode = 2; % 1: pointing, 2: yes/no
exp.conjunctionSearch = true;
exp.noTargetPossible = true;
exp.nTrialsMissingTargets = 5;
exp.numStimuli = 35;
exp.numTrials = 10;
exp.fixCrossType = 2; % 1: normal, 2: MMT star pattern
exp.fixCrossWidth = 5*pva;
exp.fixCrossLineWidth = 2; % in px
exp.fixCrossDuration = 3;
exp.fixCrossMMTn = 8;
exp.fixCrossMMTrot = 360/exp.fixCrossMMTn/2;
exp.fixCrossMMTvarRot = true;
exp.fixCrossMMTwidth = 0.2*360/exp.fixCrossMMTn;
exp.fixCrossMMTrad = 0;
exp.fixCrossMMTColor = black;
exp.jitterDistribution = true;
exp.jitterStrength = 1*pva;
exp.distributionMode = 4; %1: circle, 2: grid pattern, 3: random placement, 4: rectangular grid
exp.rectDistribution = [7,5];
exp.rectGridXDistance = 6*pva;
exp.rectGridYDistance = 4*pva;
exp.stimGridDistance = 4*pva;
exp.stimCircleRadius = 10*pva;
exp.randomStimuliSafetyDistance = 2.5*pva;

lineWidthPix = 5;

if exp.noTargetPossible
        exp.existingTargetIndex = Shuffle([zeros(1,exp.nTrialsMissingTargets),ones(1,exp.numTrials-exp.nTrialsMissingTargets)]);
    else
        exp.existingTargetIndex = ones(1,exp.numStimuli);
end

if exp.rectDistribution(1)*exp.rectDistribution(2)~=exp.numStimuli
                error(['Rectangle distribution doesn''t make sense. You tried ',num2str(exp.numStimuli),' stimuli in a ',num2str(exp.rectDistribution(1)),' by ',num2str(exp.rectDistribution(2)),' pattern']);
end

%% Search stimulus definitions
% Normal Search Distractor
exp.s1DrawType = drawTypeLine;
exp.s1Stimulus = 'E';
exp.s1Size = 1*pva;
exp.s1Rot = 0;
exp.s1LineWidth = lineWidthPix;
exp.s1Color = white;
exp.s1StimCoords = getStimCoords(exp.s1Stimulus,exp.s1Size,exp.s1Rot,exp.s1LineWidth);

% Normal Search Target
exp.s2DrawType = drawTypeLine;
exp.s2Stimulus = 'F';
exp.s2Size = 1*pva;
exp.s2Rot = 0;
exp.s2LineWidth = lineWidthPix;
exp.s2Color = white;
exp.s2StimCoords = getStimCoords(exp.s2Stimulus,exp.s2Size,exp.s2Rot,exp.s2LineWidth);

% Conjunction Search Distractor 1
exp.s3DrawType = drawTypeLine;
exp.s3Stimulus = 'E';
exp.s3Size = 1*pva;
exp.s3Rot = 0;
exp.s3LineWidth = lineWidthPix;
exp.s3Color = blue;
exp.s3StimCoords = getStimCoords(exp.s3Stimulus,exp.s3Size,exp.s3Rot,exp.s3LineWidth);

% Conjunction Search Target
exp.s4DrawType = drawTypeLine;
exp.s4Stimulus = 'F';
exp.s4Size = 1*pva;
exp.s4Rot = 0;
exp.s4LineWidth = lineWidthPix;
exp.s4Color = blue;
exp.s4StimCoords = getStimCoords(exp.s4Stimulus,exp.s4Size,exp.s4Rot,exp.s4LineWidth);

%Conjunction Search Distractor 2
exp.s5DrawType = drawTypeLine;
exp.s5Stimulus = 'F';
exp.s5Size = 1*pva;
exp.s5Rot = 0;
exp.s5LineWidth = lineWidthPix;
exp.s5Color = red;
exp.s5StimCoords = getStimCoords(exp.s5Stimulus,exp.s5Size,exp.s5Rot,exp.s5LineWidth);

%% Eyelink configuration
if exp.eyetracking
    dummymode=0;
    
    el=EyelinkInitDefaults(window);
    el.MMTfixation = true;
    el.MMTouterRadius = 5*pva;
    el.MMTinnerRadius = 0;
    el.MMTnspikes = 8;
    el.MMTspikewidth = (360/el.MMTnspikes)/3;
    EyelinkUpdateDefaults(el); % Without this no changes to el will be picked up by the callback function

    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        return;
    end
    
    [v, vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % Open file for recording data
    edfFile=[eyelinkfilename,'.edf'];
    Eyelink('Openfile', edfFile);

    % Do setup and calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    % Do a final check of calibration using driftcorrection
    % You have to hit esc before return.
    EyelinkDoDriftCorrection(el);
end

%% Show intro text
DrawFormattedText(window, 'Press any key to start the experiment', 'center', 'center', black);
Screen('Flip',window);
KbStrokeWait();

%% Start of trial loop
for t=1:exp.numTrials
    % Define a vector with number of stimuli with desired type numbers
    
    if ~exp.conjunctionSearch
        if ~exp.existingTargetIndex(t)
            exp.trial(t).stimulusTypes = ones(1,exp.numStimuli); % No Target
        else
            exp.trial(t).stimulusTypes = [1 * ones(1,exp.numStimuli-1),...
                                          2]; % Target
        end
    else
        if ~exp.existingTargetIndex(t)
            exp.trial(t).stimulusTypes = [3 * ones(1,ceil(exp.numStimuli/2)),...
                                          5 * ones(1,floor(exp.numStimuli/2))];
        else
            exp.trial(t).stimulusTypes = [3 * ones(1,ceil(exp.numStimuli/2)-1),...
                                          4,... % Target
                                          5 * ones(1,floor(exp.numStimuli/2))];
        end
        
    end
    exp.trial(t).stimulusTypes = Shuffle(exp.trial(t).stimulusTypes);
    
    switch exp.distributionMode
        case 1
            % Circle
            exp.trial(t).stimCenters = getLinSpacedPointsOnCircle(exp.numStimuli, exp.stimCircleRadius, 0);
        case 2
            % Square Grid
            exp.trial(t).stimCenters = getLinSpacedPointsOnSquareGrid(exp.numStimuli, exp.stimGridDistance);
        case 3
            % Random distribution with safety distance
            exp.trial(t).stimCenters = getRandomPointsInRect(exp.numStimuli, resolution.width, resolution.height, exp.randomStimuliSafetyDistance)'...
                - repmat([xCenter;yCenter],[1,exp.numStimuli]);
        case 4
            % Rectangular Distribution            
            exp.trial(t).stimCenters = getLinSpacedPointsOnRectGrid(exp.rectDistribution,exp.rectGridXDistance,exp.rectGridYDistance);
    end
    
    % Jitter the stimulus positions radially
    
    if exp.jitterDistribution && exp.distributionMode ~= 3
        jitterValues = NaN(2,exp.numStimuli);
        for i = 1:exp.numStimuli
            curDeg = randi(360);
            jitterValues(:,i) = [cos(curDeg); sin(curDeg)] .* rand * exp.jitterStrength;
        end
        exp.trial(t).stimCenters = exp.trial(t).stimCenters + jitterValues;
    end
    
    % Define a cell array to hold the stimuli
    stimuli = cell(4, exp.numStimuli);
    
    % Loop through all the stimuli and add their properties to the stimuli
    % cell array
    for i=1:exp.numStimuli
        switch exp.trial(t).stimulusTypes(i)
            case 1 % Normal Search Distractor
                % Define drawing type
                stimuli{1, i} = exp.s1DrawType;
                % Get stimulus coordinates
                stimuli{2, i} = exp.s1StimCoords;
                % Define color
                stimuli{3, i} = repmat(exp.s1Color,[1,size(stimuli{2, i}, 2)]);
                % Define centers of stimuli
                stimuli{4, i } = repmat(exp.trial(t).stimCenters(:,i),[1,size(stimuli{2, i}, 2)]);
                
            case 2 % Normal Search Target
                % Save target index
                exp.trial(t).targetIndex = i;
                % Define drawing type
                stimuli{1, i} = exp.s2DrawType;
                % Get stimulus coordinates
                stimuli{2, i} = exp.s2StimCoords;
                % Define color
                stimuli{3, i} = repmat(exp.s2Color,[1,size(stimuli{2, i}, 2)]);
                % Define centers of stimuli
                stimuli{4, i } = repmat(exp.trial(t).stimCenters(:,i),[1,size(stimuli{2, i}, 2)]);
                % Save target coordinates
                exp.trial(t).targetPos = exp.trial(t).stimCenters(:,i) + [xCenter yCenter]';
                
            case 3 % Conjunction Search Distractor 1
                % Define drawing type
                stimuli{1, i} = exp.s3DrawType;
                % Get stimulus coordinates
                stimuli{2, i} = exp.s3StimCoords;
                % Define colors
                stimuli{3, i} = repmat(exp.s3Color,[1,size(stimuli{2, i}, 2)]);
                % Define centers of stimuli
                stimuli{4, i } = repmat(exp.trial(t).stimCenters(:,i),[1,size(stimuli{2, i}, 2)]);
                
            case 4 % Conjunction Search Target
                % Save target index
                exp.trial(t).targetIndex = i;
                % Define drawing type
                stimuli{1, i} = exp.s4DrawType;
                % Get stimulus coordinates
                stimuli{2, i} = exp.s4StimCoords;
                % Define colors
                stimuli{3, i} = repmat(exp.s4Color,[1,size(stimuli{2, i}, 2)]);
                % Define centers of stimuli
                stimuli{4, i } = repmat(exp.trial(t).stimCenters(:,i),[1,size(stimuli{2, i}, 2)]);
                % Save target coordinates
                exp.trial(t).targetPos = exp.trial(t).stimCenters(:,i) + [xCenter yCenter]';
                
                
            case 5 %Conjunction Search Distractor 2
                % Define drawing type
                stimuli{1, i} = exp.s5DrawType;
                % Get stimulus coordinates
                stimuli{2, i} = exp.s5StimCoords;
                % Define colors
                stimuli{3, i} = repmat(exp.s5Color,[1,size(stimuli{2, i}, 2)]);
                % Define centers of stimuli
                stimuli{4, i } = repmat(exp.trial(t).stimCenters(:,i),[1,size(stimuli{2, i}, 2)]);
        end
    end
        
    % Draw fixation cross
    Screen('FillRect', window, grey);
    switch exp.fixCrossType
        case 1
            Screen('DrawLines', window,...
                     [-exp.fixCrossWidth/2, exp.fixCrossWidth/2, 0, 0;...
                      0, 0, -exp.fixCrossWidth/2, exp.fixCrossWidth/2],...
                      exp.fixCrossLineWidth, white, [xCenter yCenter], 2);
        case 2
            if exp.fixCrossMMTvarRot
                rot = randi(360);
            else
                rot = exp.fixCrossMMTrot;
            end
            fixpolycoords = getMMTfixPoly(exp.fixCrossMMTn,...
                rot,exp.fixCrossMMTwidth,exp.fixCrossMMTrad,sqrt(xCenter^2 + yCenter^2));
            for poly = 1:size(fixpolycoords,2)/2
                Screen('FillPoly', window, exp.fixCrossMMTColor, fixpolycoords(1:3,(2*poly-1):(2*poly))+repmat([xCenter yCenter],[3 1]));
            end
    end
    % Show fixation cross
    Screen('Flip',window);
    if exp.eyetracking
        Eyelink('StartRecording');
    end
    WaitSecs(exp.fixCrossDuration);

 
    %% Loop through all used drawing types and draw every groups members
    % with one command
    for drawType = unique([stimuli{1,:}])
        typeIndex = find([stimuli{1,:}] == drawType);
        % Get all stimuli with specified draw type
        curStimuli = stimuli(2:4,typeIndex);
        curStimCoords = cat(2,curStimuli{1,:});
        curStimColors = cat(2,curStimuli{2,:});
        curexp.trial(t).stimCenters = cat(2,curStimuli{3,:});
        
        switch drawType
            case drawTypeLine
                Screen('DrawLines', window, curStimCoords + curexp.trial(t).stimCenters,...
                    lineWidthPix, curStimColors, [xCenter yCenter], 2);
            case drawTypeFrameOval
                Screen('FrameOval', window , curStimColors,...
                    curStimCoords...
                    + repmat(curexp.trial(t).stimCenters...
                    + repmat([xCenter; yCenter],[1, size(curexp.trial(t).stimCenters,2)])...
                    ,[2 1]), lineWidthPix);
        end
    end
    
    % Get mouse position before showing stimulus
    % This is a workaround for click detection on the iiyama touchscreen
    % It will not work if the participant touches the top left pixel (which
    % is not really possible on most displays and no stimuli should appear
    % there)
    if exp.searchMode == 1
        SetMouse(0, 0, screenNumber);
        [xTouchscreen,yTouchscreen,buttons] = GetMouse(screenNumber);
        xTouchscreenPrev = xTouchscreen;
        yTouchscreenPrev = yTouchscreen;
    end
    
    %% Show the stimuli
    Screen('Flip', window);
    tFlip = GetSecs();
    
    if exp.eyetracking
        Eyelink('Message', 'SYNCTIME');
    end
    
    %% Check for a response
    switch exp.searchMode
        case 1
            while xTouchscreen == xTouchscreenPrev && yTouchscreen == yTouchscreenPrev
                xTouchscreenPrev = xTouchscreen;
                yTouchscreenPrev = yTouchscreen;
                [xTouchscreen,yTouchscreen,buttons] = GetMouse(screenNumber);
            end
        case 2
            RestrictKeysForKbCheck([KbName('y'),KbName('n')]);
            [~, keyCode, ~] = KbPressWait();
            if keyCode(KbName('y'))
            	exp.trial(t).yesnoResponse = 'yes';
            elseif keyCode(KbName('n'))
                exp.trial(t).yesnoResponse = 'no';
            end
            RestrictKeysForKbCheck([]);
    end
    
    tTouch = GetSecs();
    
    %% Give the arduino the signal to send a 5v pulse to the AD board that
    % finishes the recording of each activity    
    if exp.motiontracking
        fprintf(arduino,'%c','p');
    end
    
    if exp.eyetracking
        Eyelink('Message', 'TRIAL_END');
        Eyelink('StopRecording');
    end
    
    %% Save trial related data
    exp.trial(t).stimuli = stimuli;
    if exp.searchMode == 1
        exp.trial(t).touchPos = [xTouchscreen yTouchscreen]';
        if exp.existingTargetIndex(t)
            exp.trial(t).touchError = exp.trial(t).touchPos - exp.trial(t).targetPos;
        else
            exp.trial(t).touchError = NaN;
        end
    elseif exp.searchMode == 2
        if exp.existingTargetIndex(t)
            exp.trial(t).targetExistent = 'yes';
        else
            exp.trial(t).targetExistent = 'no';
        end
    end
    exp.trial(t).tTouch = tTouch;
    exp.trial(t).tFlip = tFlip;
    exp.trial(t).RT = tTouch - tFlip;
    
    % Save a screenshot of the stimuli for easier visual analysis later
    exp.trial(t).Screenshot = Screen('GetImage', window);
    
    % Save mat file
    save([cd,'\results\',filename,'.mat'],'exp');
    
    [keyIsDown,secs,keyCode] = KbCheck; 
    if keyCode(KbName('shift')) % Quit if shift is pressed
        break
    end

    if keyCode(KbName('control')) % Pause if control is pressed
        DrawFormattedText(window, 'Experiment paused\nPress any key to continue', 'center', 'center', black);
        Screen('Flip',window);
        WaitSecs(3)
        KbStrokeWait();
    end
end

%% Finish the experiment
DrawFormattedText(window, 'Experiment finished.\nClosing...', 'center', 'center', black);
Screen('Flip',window);
WaitSecs(1);

% Show the cursor again
ShowCursor(screenNumber);

% Release the arduino (when stopping the script before this point, you
% might need to delete(instrfindall) )
if exp.motiontracking
    fclose(arduino);
end

if exp.eyetracking
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    WaitSecs(0.5);
    try
        status=Eyelink('ReceiveFile',edfFile,[cd '\results\'],1);
        if status > 0
            disp(['ReceiveFile status ', status]);
        end
        if 2==exist(edfFile, 'file')
            disp(['Data file can be found in ',pwd]);
        end
    catch %#ok<*CTCH>
        disp('Problem receiving data file');
    end
    
    Eyelink('Shutdown');
end

% Clear the screen.
sca;