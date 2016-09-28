% Experiment to explore visual search strategies, mainly in patients with
% age-related macular degeneration (AMD). Multiple search paradigms can be
% selected or changed quickly to be able to react to the specific visual
% abilities patients have left.
%
% The experiment includes options to enable motion and gaze tracking. The
% gaze tracking assumes a connected EyeLink 1000 device. The motion
% tracking assumes an Arduino connected via USB, set up to wait for input
% of a 'p' to send a 5v pulse to an A/D board connected to a PC running the
% MotionMonitor software, to trigger recording of motion data in sync with
% the experiment.
%
% A star pattern to infer a fixation point from centrally converging spikes
% is implemented to allow patients with losses of visual ability in parts
% of their visual field to fixate centrally without being able to see the
% center of the fixation cross.
%
% Programmed by Julius Krumbiegel at Ludwig-Maximilians-Universitaet
% Munich starting in July 2016


% Clear the workspace and the screen
sca;
close all;
clearvars;

% Turn motion tracking on or off. If turned on, this will send the trigger
% signal to the arduino which gives the MotionMonitor its recording trigger
exp.motiontracking = false;
% Turn eye tracking on or off. If this is turned on but no working eye
% tracker is connected, you are asked if you want to start the experiment
% with eye tracking dummy mode.
exp.eyetracking = false;

% Define a folder in which the experiment and eyetracking data will be
% stored. Create the folder if it doesn't exist.
resultsFolder = 'C:\Users\ru35pec\Desktop\Julius\Macular Degeneration\AMD Experimente\results';
[status,message,messageid] = mkdir(resultsFolder);
if ~status
    error('The specified results folder doesn''t exist and couldn''t be created');
end

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

% Define some colors for convenience
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
pva = pixPerCM / va; % This value is used throughout the script to set distances in visual angles in functions that require pixel units
exp.pixPerVA = pva;

%% Set trial parameters
exp.searchMode = 2; % 1: pointing, 2: yes/no
exp.conjunctionSearch = true; % Enable or disable conjunction search (uses stimuli 3,4,5 instead of 1,2)
exp.noTargetPossible = true; % Enable or disable trials without target (don't use in pointing mode)
exp.nTrialsMissingTargets = 5; % Number of trials without target
exp.numStimuli = 35; % Number of stimuli per trial (target and distractors)
exp.numTrials = 10; % Number of trials per experiment

exp.fixCrossType = 2; % Type of fixation cross. 1: normal, 2: MMT star pattern (Macular Mapping Test, to infer fixation point from converging lines)
exp.fixCrossWidth = 5*pva; % Width of standard fixation cross in pixels
exp.fixCrossLineWidth = 2; % Line width of standard fixation cross in pixels
exp.fixCrossDuration = 3; % Fixation cross duration in seconds
exp.fixCrossMMTn = 8; % Number of star spikes in the MMT fixation cross
exp.fixCrossMMTrot = 360/exp.fixCrossMMTn/2; % Rotation of the MMT fixation cross in degrees
exp.fixCrossMMTvarRot = true; % If true, the MMT fixation cross is randomly rotated in each trial. This takes precedence over the fixed rotation setting
exp.fixCrossMMTwidth = 0.2*360/exp.fixCrossMMTn; % Width of each MMT fixation cross star spikes in degrees
exp.fixCrossMMTinnerRad = 0; % Inner radius of the MMT fixation cross pattern in pixels. If 0, the spikes meet in the center point.
exp.fixCrossMMTColor = black; % Color of the MMT fixation cross

exp.jitterDistribution = true; % Enable or disable radial area jitter of the stimulus positions. Doesn't have an effect on randomly placed stimuli
exp.jitterStrength = 1*pva; % Maximum radial displacement by the stimulus position jitter in pixels
exp.distributionMode = 4; % Distribution pattern for stimulus positions. 1: circle, 2: square grid, 3: random placement, 4: rectangular grid
exp.rectDistribution = [7,5]; % Dimensions of the rectangular grid in X and Y. The multiple of both must equal the number of stimuli
exp.rectGridXDistance = 6*pva; % Horizontal distance between stimuli in the rectangular grid pattern in pixels
exp.rectGridYDistance = 4*pva; % Vertical distance between stimuli in the rectangular grid pattern in pixels
exp.stimGridDistance = 4*pva; % Horizontal and vertical distance between stimuli in the square grid pattern in pixels
exp.stimCircleRadius = 10*pva; % Radius of the circle pattern in pixels
exp.randomStimuliSafetyDistance = 2.5*pva; % Minimum distance between the points generated in the random grid pattern in pixels. This distance will be kept to the boundary rectangle borders as well

lineWidthPix = 5; % A standard line width in pixels that may be used in stimulus definitions below

% Create a vector whose elements indicate which trials in the experiment
% have targets (1) and which don't (0)
if exp.noTargetPossible
        exp.existingTargetIndex = Shuffle([zeros(1,exp.nTrialsMissingTargets),ones(1,exp.numTrials-exp.nTrialsMissingTargets)]);
    else
        exp.existingTargetIndex = ones(1,exp.numStimuli);
end

% Check if rectangle dimensions and number of stimuli agree
if exp.rectDistribution(1)*exp.rectDistribution(2)~=exp.numStimuli
                error(['Rectangle distribution doesn''t make sense. You tried ',num2str(exp.numStimuli),' stimuli in a ',num2str(exp.rectDistribution(1)),' by ',num2str(exp.rectDistribution(2)),' pattern']);
end

%% Search stimulus definitions

% Normal Search Distractor
exp.s1DrawType = drawTypeLine; % The drawing type used for the selected stimulus. Currently, only stimuli drawn with a single drawing routine are supported
exp.s1Stimulus = 'E'; % The chosen stimulus. Must be one of the predefined options in getStimCoords()
exp.s1Size = 1*pva; % The size of the stimulus which is multiplied with the stimulus coordinates to scale it
exp.s1Rot = 0; % The rotation of the stimulus in degrees
exp.s1LineWidth = lineWidthPix; % The line width of the stimulus, only relevant for drawing routines relying on lines instead of filled shapes. Is also relevant for some stimuli whose coordinates have to be adjusted to their line thickness
exp.s1Color = white; % The color of the stimulus. Currently, no multicolored stimuli are supported
exp.s1StimCoords = getStimCoords(exp.s1Stimulus,exp.s1Size,exp.s1Rot,exp.s1LineWidth); % Get the desired stimulus coordinates by feeding the specified settings into the getStimCoords() function

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
    % Set eyetracking dummymode. Normal eyetracking if set to 0, dummymode
    % if set to 1
    dummymode=0;
    
    % Initialize the eyelink data structure with standard settings
    el=EyelinkInitDefaults(window);
    
    % Write custom settings to the eyelink data structure. These will only
    % be used if the customized PsychEyelinkDispatchCallback.m file in the
    % folder C:\toolbox\Psychtoolbox\PsychHardware\EyelinkToolbox\EyelinkBasic
    % is present. Otherwise none of these settings will have any effect but
    % they shouldn't disrupt the process
    el.MMTfixation = true; % Enable or disable the MMT star calibration pattern
    el.MMTouterRadius = 5*pva; % The outer radius of the MMT star pattern in pixels
    el.MMTinnerRadius = 0; % The inner radius of the MMT star pattern in pixels
    el.MMTnspikes = 8; % The number of star spikes in the MMT star pattern
    el.MMTspikewidth = (360/el.MMTnspikes)/3; % The width of each MMT star spike in degrees
    
    % Without this no changes to the eyelink data structure will be picked
    % up by the callback function
    EyelinkUpdateDefaults(el); 
    
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

%% Show intro text and continue after any key is pressed
DrawFormattedText(window, 'Press any key to start the experiment', 'center', 'center', black);
Screen('Flip',window);
KbStrokeWait();

%% Loop through all trials
for t=1:exp.numTrials
    % Define a vector with as many entries as there are stimuli, where
    % the type of stimulus is specified by its code number. These are:
    % 1: normal search distractor
    % 2: normal search target
    % 
    % 3: conjunction search distractor 1
    % 4: conjunction search target
    % 5: conjunction search distractor 2
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
    
    % Jitter the stimulus positions radially, but not if random mode is selected    
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
    % cell array, which will be the input for the psychtoolbox drawing
    % routines further down
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
        % This draws a simple fixation cross made of two intersecting lines
        case 1
            Screen('DrawLines', window,...
                     [-exp.fixCrossWidth/2, exp.fixCrossWidth/2, 0, 0;...
                      0, 0, -exp.fixCrossWidth/2, exp.fixCrossWidth/2],...
                      exp.fixCrossLineWidth, white, [xCenter yCenter], 2);
                  
        % This draws a star pattern which allows patients with imperfect
        % vision to infer the fixation point from the converging lines. It
        % will randomly rotate between trials if the respective option is
        % turned on to reduce afterimages from constant exposure to the
        % same stimulus.
        case 2
            if exp.fixCrossMMTvarRot
                rot = randi(360);
            else
                rot = exp.fixCrossMMTrot;
            end
            
            % Get the star pattern polygon coordinates for a star pattern
            % with the parameters specified in the settings section
            fixpolycoords = getMMTfixPoly(exp.fixCrossMMTn,...
                rot,exp.fixCrossMMTwidth,exp.fixCrossMMTinnerRad,sqrt(xCenter^2 + yCenter^2));
            for poly = 1:size(fixpolycoords,2)/2
                Screen('FillPoly', window, exp.fixCrossMMTColor, fixpolycoords(1:3,(2*poly-1):(2*poly))+repmat([xCenter yCenter],[3 1]));
            end
    end
    % Show fixation cross
    Screen('Flip',window);
    
    % Send the eyelink computer the signal to start recording gaze data
    if exp.eyetracking
        Eyelink('StartRecording');
    end
    
    % Show the fixation cross for the amount of time specified in the
    % settings section
    WaitSecs(exp.fixCrossDuration);

 
    %% Loop through all used drawing types and draw all groups' members
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
        % Set the cursor position (invisibly) to the top left pixel
        SetMouse(0, 0, screenNumber);
        % Read cursor position to ensure correct behavior
        [xTouchscreen,yTouchscreen,buttons] = GetMouse(screenNumber);
        xTouchscreenPrev = xTouchscreen;
        yTouchscreenPrev = yTouchscreen;
    end
    
    %% Show the stimuli
    Screen('Flip', window);
    
    % Get the time when the search stimuli were shown
    tFlip = GetSecs();
    
    % Send the eyelink computer a marker that the search array is now
    % visible
    if exp.eyetracking
        Eyelink('Message', 'SYNCTIME');
    end
    
    %% Check for a response
    switch exp.searchMode
        case 1
        % Pointing task
        % In this case a touchscreen cursor position change will mark the
        % reaction of the subject 
            while xTouchscreen == xTouchscreenPrev && yTouchscreen == yTouchscreenPrev
                xTouchscreenPrev = xTouchscreen;
                yTouchscreenPrev = yTouchscreen;
                [xTouchscreen,yTouchscreen,buttons] = GetMouse(screenNumber);
            end
        case 2
        % Yes / No task
        % In this case pressing the y or n key will mark the reaction of
        % the subject
            RestrictKeysForKbCheck([KbName('y'),KbName('n')]);
            [~, keyCode, ~] = KbPressWait();
            if keyCode(KbName('y'))
            	exp.trial(t).yesnoResponse = 'yes';
            elseif keyCode(KbName('n'))
                exp.trial(t).yesnoResponse = 'no';
            end
            RestrictKeysForKbCheck([]);
    end
    % Record reaction time
    tTouch = GetSecs();
    
    %% Give the arduino the signal to send a 5v pulse to the AD board that
    % finishes the recording of each activity    
    if exp.motiontracking
        fprintf(arduino,'%c','p');
    end
    
    % Send the eyelink computer the signal to stop the current recording
    if exp.eyetracking
        Eyelink('Message', 'TRIAL_END');
        Eyelink('StopRecording');
    end
    
    %% Save trial related data
    % Save all stimuli parameters and coordinates in exp structure
    exp.trial(t).stimuli = stimuli;
    
    % Pointing task
    % Save touch positions, but not if no target was present in a trial
    if exp.searchMode == 1
        exp.trial(t).touchPos = [xTouchscreen yTouchscreen]';
        if exp.existingTargetIndex(t)
            exp.trial(t).touchError = exp.trial(t).touchPos - exp.trial(t).targetPos;
        else
            exp.trial(t).touchError = NaN;
        end
        
    % Yes / No task
    % Save yes no response, but no touch position
    elseif exp.searchMode == 2
        if exp.existingTargetIndex(t)
            exp.trial(t).targetExistent = 'yes';
        else
            exp.trial(t).targetExistent = 'no';
        end
    end
    
    % Save touch times and flip times into exp structure, calculate
    % reaction time as well
    exp.trial(t).tTouch = tTouch;
    exp.trial(t).tFlip = tFlip;
    exp.trial(t).RT = tTouch - tFlip;
    
    % Save a screenshot of the stimuli for easier visual analysis later
    exp.trial(t).Screenshot = Screen('GetImage', window);
    
    % Save exp struct to a mat file in the results folder
    save([resultsFolder,filename,'.mat'],'exp');
    
    % Quit if shift is pressed at trial end
    [keyIsDown,secs,keyCode] = KbCheck; 
    if keyCode(KbName('shift')) 
        break
    end

    % Pause if control is pressed at trial end
    if keyCode(KbName('control')) 
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

% Receive the recorded eyelink data from eyelink computer. This will fail
% if the chosen folder does not exist yet. The eyelink computer will flash
% red text and show the file explorer in that case.
if exp.eyetracking
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    WaitSecs(0.5);
    try
        status=Eyelink('ReceiveFile',edfFile,resultsFolder,1);
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