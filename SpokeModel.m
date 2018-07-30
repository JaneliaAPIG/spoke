classdef SpokeModel < most.Model
    %SPIKEGRID Application class for Spoke; all its functionality is implemplemented here
    % See project wiki for both user & developer overview documentation 
    % See end-of-file for nitty-gritty developer notes
    
    
    
    %% PUBLIC PROPERTIES
    properties (SetObservable)
        refreshRate = 5; %Refresh rate, in Hz, at which data is reviewed and newly detected spikes are plotted
        
        thresholdType = 'rmsMultiple'; %One of {'volts' 'rmsMultiple'}. If rmsMultiple is specified, it is assumed that signals' DC is already filtered.
        thresholdVal = 5; %Threshold, in volts or SD units, to use for spike detection
        thresholdAbsolute = false; % Logical indicating whether threshold should be considered an absolute value
        
        baselineStatsRefreshPeriod = 2; %Period, in seconds, at which signal baseline stats (mean & RMS) values are recomputed, when needed (e.g. for 'rmsMultiple' threshold detection)
        %baselineStatsRefreshOnRetrigger = true; %Logical indicating whether RMS value should be recomputed on SpikeGL 'retriggers' (applies for 'rmsMultiple' threshold detection)
        
        globalMeanSubtraction = false; %Logical indicating whether to compute/apply global mean subtraction
        
        dataReadMode = 'file'; %One of {'file' 'spikeGL'}. Specifies whether to read data from file (faster, but assumes only one file) or via spikeGL function.
        
        filterWindow = [0 inf]; %Frequency range in Hz allowed to pass, using a 1-pole Butterworth bandpass filter
        
        tabDisplayed = 1; %Number of 'tab' -- a batch of (PLOTS_PER_TAB) channels -- currently diplayed in the spike waveform grid.
        
        displayMode = 'waveform'; %One of {'waveform' 'raster'}. Specifies type of information to display on plot grid.
        
        %Spike waveform display properties
        horizontalRange = [-1 2] * 1e-3; %2 element vector ([pre post]) indicating times, in seconds, to display before and after threshold crossing
        waveformAmpUnits = 'volts'; %One of {'volts' 'rmsMultiple'} indicating units to use for waveform plot display. Value of 'rmsMultiple' only available if thresholdType='rmsMultiple'
        
        waveformsPerPlot = 100; %Number of waveforms to display (overlay) in each channel subplot before clearing
        waveformsPerPlotClearMode = 'all'; %One of {'all' 'oldest'}. Specify waveform-clearing behavior when the waveformsPerPlot limit is reached.
        waveformPlotClearPeriod = inf; % (TODO) Time, in seconds, after which to clear all or oldest waveform if no waveforms have been received
        
        %Raster/PSTH display properties
        stimStartChannel = []; %Channel number (zero-indexed) to use for signaling start of stim
        stimStartThreshold = 0.5; %Value, in volts, to use for stim start signal
        horizontalRangeRaster = [-1 1]; %2 element vector ([pre post]) indicating times, in seconds, to display before and after stim start threshold crossing in raster/psth plots
        verticalRangeRaster = [1 inf]; %2 element vector indicating which stim numbers to display, counted from last start/reset. Set second element to Inf to specify all stimuli should be displayed.
        verticalRangeRasterInfIncrement = 20; %When verticalRangeRaster(2)=Inf, specify the increment in which the number of stimuli displayed are incremented by.
        stimEventTypesDisplayed = {}; %String or string cell array indicating which event type(s) to display raster/PSTH data for.
        
        spikeRefractoryPeriod = 2e-3; %Time, in seconds, to prevent detection of spike following previously detected spike. When displayMode='waveform', this value is given by second element of horizontalRange.
        
        psthTimeBin = 10e-3; %Time, in seconds, over which to bin spikes for PSTH summary plot
        psthAmpRange = [0 120]; %Amplitude range to display, in units of spikes/second
        
        psthTimerPeriod = inf; %Period, in seconds, at which plotPSTH() method is called automatically when displayMode='raster'
        
        %channelSubset = inf; %DEPRECATED BY HIDDEN DISPLAYCHANNELS PROP Subset of channels to acquire from
        
        % The following are properties that were SetAccess=protected, but
        % have been moved out of protected to allow config file saves with
        % most.
        refreshPeriodMaxWaveformRate = inf; %Maximum waveform rate (Hz) to detect/plot for each refresh period; spikes/stimuli above this rate are discarded (in spike-/stim-triggered waveform modes, respectively)
        
        % The following are properties that back up dependent properties.
        % This is for properties that need to be saved to disk.
        verticalRange_; % backs up property verticalRange
    end
    
    properties (SetObservable, Transient)
        stimEventClassifyFcn = []; %Function handle specifying function used to classify stimuli into event types. See class notes for further information.
    end
    
    properties (SetObservable, Dependent)
        verticalRange; %2 element vector ([min max]) indicating voltage bounds or RMSMultiple (depending on thresholdType) for each spike plot
        refreshPeriodMaxNumWaveforms = inf; %Maximum number of spikes to detect/plot during a given refresh period
        numAuxChans; %Number of auxiliary
    end
    
    properties (Dependent)
        numTabs;
        stimEventCount; %Count of stimuli that have been detected since start/restart/rollover at current stimEventTypesDisplayed (when displayMode='raster')
    end
    
    %Following are determined at start of acquisition, via the
    %stimEventClassifyFcn (no arg case). This function must return a structure
    %with fields stimEventTypes & stimEventClassifyTime.
    properties (SetAccess=protected)
        stimEventTypes; %String cell array indicating names of distinct event types which can be signaled by stimulus start pulse. If empty, event Classify is not used
        stimEventClassifyNumScans = 0; %Number of scans required to classify event type using stimulus start signal.
        stimEventClassifyChannel = []; %SpikeGL (demuxed) channel number containing signal used for stimulus event classification (immediately following stimulus detection). Must be an auxiliary channel. If empty, the stimStartChannel is used.
        
        stimTotalCount = 0; %Count of stimuli that have been detected since start/restart/rollover, for /all/ stimEventTypes (when displayMode='raster')
        
        configFileName = ''; %Full filename (incl. path) of last-saved/loaded config file, if any
    end
    
    properties (SetAccess=private,SetObservable)
        running = false;
    end
    
    %% HIDDEN PROPERTIES
    
    properties (Hidden)
        maxBufSizeSeconds = 1;
    end
    
    
    properties (Hidden,SetAccess=protected)
        hSGL; %Handle to SpikeGLX object
        sglParamCache; %cached SpikeGLX parameters
        
        hTimer; %Timer to use for periodically checking Spoke's input data stream
        hPSTHTimer; %Timer to periodically call plotPSTH() method
        
        %Following structs have 3 keys {'waveform' 'raster' 'psth'}
        hFigs; %Struct of handles to grid figure
        hPanels; %Struct of handles to subpanels in grid figure, one for each channel per tab
        hChanLabels; %Struct of arrays of text labels, identifying each channel
        hButtons; %Struct of handles to buttons on the waveform and raster grids.
        
        %Handle graphics specific to waveform display
        hPlots; %Array of axes handles, one for each axes in grid
        hThresholdLines; %Cell array of line handles, marking threshold level for each plot in grid
        hWaveforms; %Array of animated line objects for each plotted waveform
        
        %Handle graphics specific to raster display
        hRasters; %Array of axes handles, one for each axes in grid
        hPSTHs; %Array of axes handles, one for each axes in grid
        hRasterLines; %Cell array of arrays, one per stimEventType, containing animated line objects for each raster plot in grid
        
        numActiveTabs;

        %neuralChans; %Acquisition channel numbers corresponding to neural inputs. Ordered as in SpikeGLX connection.
        %auxChans; %Acquisition channel numbers corresponding to auxiliary inputs, i.e. suited for gating/stimulus. These are /not/ displayed. Ordered as in SpikeGLX connection.
        sglChanSubset; %subset of SpikeGLX channels  (includes neural and non-neural channels)
        
        neuralChanDispOrder; %Order in which neuralChans should be displayed. Includes all acquired neural chans, whether included in subset or not.
        
        neuralChanAcqList; % List of neural chans included in SpikeGLX chan subset
        neuralChanDispList; %Ordered list of neural chans displayed. Display ordering applied. Channel subset, if any, is applied.
        auxChanProcList; %List of auxiliary chans processed. Channel subset, if any, is applied.
        
        baselineRMS; %Array of RMS values, one per channel, computed from all non-spike-window scans from within last baselineRMSTime
        baselineMean; %Array of mean values, one per channel, computed from all non-spike-window scans from within last baselineRMSTime
        baselineRMSLastScan = 0; %Last scan number at which threshold RMS value was updated
        
        reducedData = {}; %Cell array, one cell per available neural channel, of structures storing reduced data from the full data (timestamps, waveform snippets, stimulus-tagged timestamps), based on detected stimuli and/or spike(s)
        lastPlottedWaveformCount = []; %Scalar array specifying the spike count, for each channel that's been last plotted
        lastPlottedWaveformCountSinceClear = []; %Scalar array specificying the spike count, for each channel that's been last plotted (cleared every time a channel is cleared) - only used for 'all' plot clear mode.
        %globalMean; %Global mean value across acquisition channels, computed from all non-spike-window scan from within last baselineRMSTime, if globalMeanSubtraction=true is nonempty
        
        %baselineRMSExcludedScans; %Array of numbers, one per channel, indicating number of scans included into RMS calculation, from within last baselineRMSTime
        %baselineRMSNumSamples;
        
        bufScanNumEnd; %Scan number of the last element in the fullDataBuffer
        fullDataBuffer; %Array to cache full channel data to process during each timer cycle. Grows & contracts each cycle.
        
        voltsPerBitAux; %Scaling factor between A/D values and voltage for Auxiliary Channels
        voltsPerBitNeural; %Scaling factor between A/D values and voltage for Neural Channels
        
        filterCoefficients = {}; %Cell array of a,b filter coefficients for current filterWindow
        filterCondition; %Maintain 'initial condition' of filter between calls to filter()
        
        tabChanNumbers; %Array of pad channel numbers (one-based) for currently specified tabDisplayed value
        
        blockTimer = false; %Flag used to block timer callback actions during certain operations
        
        maxReadableScanNum = 0; %Maximum scan number in currently logged-to file counted since start of acquisition (count uninterrupted across start & retrigger gaps)
        lastMaxReadableScanNum = 0;  %Value of maximum scan number in currently logged-to file, measured from start of file, at /last/ timer period
        priorfileMaxReadableScanNum = 0; %Value of maximum scan number, measured from start of acquisition (uninterrupted count), reached for /prior/ file
        
        %Following arrays are all the same length (could be a structure, but kept separate for ease of access/indexing)
        stimScanNums; %Array of scan numbers corresponding to detected stimulus edges
        stimWindowStartScanNums; %Array of scan numbers marking starting edge of stimulus display events
        stimWindowEndScanNums; %Array of scan numbers marking ending edge of stimulus display events
        stimEventTypeNames = {}; %Cell string array of event type names
        
        stimNumsPlotted; %Structure array, of size neuralChansAvailable and with fields given by stimEventTypeNames, indicating number of stims that have been plotted so far for each event
        stimLastEventScanNumWindow; %1x2 array indicating start/end scan numbers for last stimulus trial
        
        stimEventCount_; %Struct var containing stimEventCount value for each of the stimEventTypes_
        
        bmarkReadTimeStats = zeros(3,1); %Array of [mean std n]
        bmarkPreProcessTimeStats = zeros(3,1);; %Array of [mean std n]
        bmarkPlotTimeStats = zeros(3,1);; %Array of [mean std n]
        bmarkPostProcessTimeStats = zeros(3,1);; %Array of [mean std n]
        
        %waveformWrap/partialWaveformBuffer props handle edge-cases in stim-triggered waveform mode. they could potentially be used for spike-triggered waveform mode.
        waveformWrap = []; %waveform detected towards end of timer processing period; specifies number of samples in the next processing period needed to complete the waveform
        partialWaveformBuffer = {}; % Buffer that holds part of a waveform from prior processing period. Used when waveformWrap is true.
        
        verticalRangeCache; % A struct used to persist the last value used for verticalRange in each mode.

        % Debug related properties.
        diagramm; %Temporary figure for debugging.
        diagrammm; %Temporary figure 2 for debugging.
        debug = false; % Set to true to enable debugging.
    end
    
    properties (Hidden, SetAccess=immutable)
        sglIPAddress; %IP Address for remote control connection to SpikeGLX
        
        %Number of logical chans of each of the sub-types available, as configured via the SpikeGLX NI Configuration
        neuralChansAvailable;
        analogMuxChansAvailable;
        analogSoloChansAvailable;        
    end
    
    properties (SetAccess=protected,Hidden,SetObservable,AbortSet)
        maxNumWaveformsApplied = false; %Logical indicating if the refreshPeriodMaxNumWaveforms clamp was applied for any channel on the last refresh period
        
    end
    
    properties (Hidden, Dependent)
        horizontalRangeScans; %horizontalRange in 'scan' units, which correspond to A/D 'scans' (a sample for each channel). Note that the window includes this # of pre & post scan, PLUS one additional scan for the spike itself
        
        baselineStatsRefreshPeriodScans; %baselineStatsRefreshPeriod in scan units
        
        maxBufSizeScans; %Maximum number of scans to process during a refresh period
        refreshPeriodAvgScans; %Average number of scans to process at each refresh period
        
        displayModeAxes; %Returns either hPlots/hRaster, depending on displayMode value
        
        stimEventTypes_; %Returns {'allstim'} in case that stimEventTypes is empty
        
        gridFigPosition; %Figure position of raster/waveform grid figures (same position for both..only one shown at a time)
        psthFigPosition; %Figure position of PSTH grid figure
        
        maxPointsPerAnimatedLine; %Used to set the number of max points per animated line.
                        
        %neuralChanAcqList; %Used to get the number of neural channels (used instead of sglChanSubset, which returns all channels, not just MN chans)
    end
    
    %Constants
    
    
    % ********** NOTES *******
    % 1. Is it possible to change plot order in tab to 2x16 (column)
    % 2. Can we do 2x64?
    % 3. 
    
    properties (Hidden,Constant)
        PLOTS_PER_TAB = 32;
        MAX_NUM_TABS = 8;
        
        INIT_RMS_THRESHOLD = 10; %Initial rmsMultiple to use when detecting spikes to exclude from initial RMS value determination
        RASTER_DISP_STIM_RANGE_INCREMENT_FRACTION = 0.85; %fraction of verticalRangeRasterInfIncrement which must be reached before display range is auto-incremented.
        
        SGL_BITS_PER_SAMPLE = 16; %A constant from SpikeGLX currently; should perhaps update SpikeGL to pull this from the DAQmx API
    end
    
    %% CONSTRUCTOR/DESTRUCTOR
    methods
        
        function obj = SpokeModel(sglIPAddress)
            
            %Process inputs
            obj.sglIPAddress = sglIPAddress;
            obj.hSGL = SpikeGL(sglIPAddress);
            
            obj.sglParamCache = GetParams(obj.hSGL); %initialize SpikeGL param cache
            
            %Create class-data file
            s.lastConfigFileName = '';
            obj.ensureClassDataFile(s,mfilename('class'));
            
            %Initialize resources
            obj.ziniCreateGrids(); %Create spike waveform grid figure
            
            obj.hTimer = timer('Name','Spoke Waveform Grid Timer','ExecutionMode','fixedRate','TimerFcn',@obj.zcbkTimerFcn,'BusyMode','queue','StartDelay',0.1);
            obj.hPSTHTimer = timer('Name','Spoke Plot PSTH Timer','ExecutionMode','fixedRate','TimerFcn',@(src,evnt)obj.plotPSTH,'BusyMode','drop','StartDelay',0.1);
            
            %Immutable prop initializations
            [obj.neuralChansAvailable, obj.analogMuxChansAvailable, obj.analogSoloChansAvailable] = obj.zprvGetAvailAcqChans();
            
            %Programmatic prop intializations
            aiRangeMax = obj.sglParamCache.niAiRangeMax;
            niMNGain = obj.sglParamCache.niMNGain;
            niMAGain = obj.sglParamCache.niMAGain;
            
            %obj.voltsPerBitNeural = aiRangeMax / 2^(obj.SGL_BITS_PER_SAMPLE - 1) / niMNGain;
            obj.voltsPerBitNeural = (aiRangeMax / niMNGain) / ( 2^(obj.SGL_BITS_PER_SAMPLE - 1));
            obj.voltsPerBitAux = (aiRangeMax / niMAGain) / ( 2^(obj.SGL_BITS_PER_SAMPLE - 1));
            obj.refreshRate = obj.refreshRate; %apply default value
            
            obj.zprvResetReducedData();
            
            %Initialize a default display for appearances (some aspects gets overridden by processing on start())
            %             obj.sglChanSubset = GetChannelSubset(obj.hSGL); %channel subset as specified in SpikeGLX. Wierd - this /has/ to be done here, outside of zprvZpplyChanOrderAndSubset() to avoid a hang.
            % - MOVED TO zprvApplyChanOrderAndSubset. obj.sglChanSubset = GetSaveChansNi(obj.hSGL); %channel subset as specified in SpikeGLX. Wierd - this /has/ to be done here, outside of zprvZpplyChanOrderAndSubset() to avoid a hang.
            obj.zprvApplyChanOrderAndSubset();
            
            numNeuralChans = numel(obj.neuralChansAvailable);
            obj.hThresholdLines = repmat({ones(numNeuralChans,1) * -1},2,1);
            
            %Allocate spike waveforms & raster plots
            obj.hWaveforms = gobjects(numNeuralChans,1);
            for i=1:obj.PLOTS_PER_TAB
                obj.hWaveforms(i) = animatedline('Parent',obj.hPlots(i),'Color','k','MaximumNumPoints',Inf,'Marker','.','MarkerSize',3,'LineStyle','-');
            end
            obj.zprvInitializeRasterGridLines();
            
            obj.verticalRange = [-aiRangeMax aiRangeMax] / 1000; % Use millivolts for spikes.
            obj.tabDisplayed = 1;
            
            %Clean-up
            Close(obj.hSGL);
            obj.hSGL = [];
            
        end
        
        function initialize(obj)
            initialize@most.Model(obj);
            
            if ~isempty(obj.getClassDataVar('lastConfigFileName'))
                obj.loadConfig();
            end
        end
        
        function delete(obj)
            
            %delete(obj.hFigs); %Why doesn't this happen automatically?
            % Can't do the above because Matlab won't let you delete handle struct.
            
            obj.quit();
        end
        
    end
    
    methods (Hidden)
        %Helper methods
        
        function ziniCreateGrids(obj)
            
            numNeuralChans = numel(obj.neuralChansAvailable);
            
            obj.numActiveTabs = ceil(numNeuralChans/obj.PLOTS_PER_TAB);
            assert(obj.numActiveTabs <= obj.MAX_NUM_TABS,'Exceeded maximum number of tabs (%d)',obj.MAX_NUM_TABS); %TODO: Deal more gracefully
            
            %This is a good idea...but let's keep it simple for now
            %       if obj.numActiveTabs > 1
            %         numPlots = obj.PLOTS_PER_TAB;
            %       else
            %         numPlots = numDispChans;
            %       end
            numPlots = obj.PLOTS_PER_TAB;
            
            gridDimension = ceil(sqrt(numPlots));
            gridPanelSize = 1/gridDimension;
            
            gridDimensionX = 2;
            gridDimensionY = 16;
            
            gridPanelSizeX = 1/gridDimensionX;
            gridPanelSizeY = 1/gridDimensionY;
            
            %Create waveform & raster grids
            obj.hFigs.waveform = most.idioms.figureScaled(1.6,'Name','Spoke Waveform Grid','CloseRequestFcn',@(src,evnt)set(src,'Visible','off'));
            obj.hFigs.raster = most.idioms.figureScaled(1.6,'Name','Spoke Raster Grid','CloseRequestFcn',@(src,evnt)set(src,'Visible','off'));
            obj.hFigs.psth = most.idioms.figureScaled(1.6,'Name','Spoke PSTH Grid','Visible','off','CloseRequestFcn',@(src,evnt)set(src,'Visible','off'));
            structfun(@(hFig)set(hFig,'NumberTitle','off','Menubar','none','Toolbar','none'),obj.hFigs);
            
            set([obj.hFigs.waveform obj.hFigs.raster obj.hFigs.psth],'Units','normalized');
            
            %Add pushbuttons to the waveform and raster grids.
            obj.hButtons.toolbar = uitoolbar(obj.hFigs.waveform);
            
            % Read an image
            [img,map] = imread(fullfile(matlabroot, 'toolbox','matlab','icons','matlabicon.gif'));

            % Convert image from indexed to truecolor
            icon = ind2rgb(img,map);

            obj.hButtons.screenshot = uipushtool(obj.hButtons.toolbar,'TooltipString','Click to Save a Screenshot',...
                 'ClickedCallback',...
                 'hGrid.zprvSnapshot');
            obj.hButtons.pause =  uipushtool(obj.hButtons.toolbar,'TooltipString','Click to Pause or Unpause Operation',...
                 'ClickedCallback',...
                 'hGrid.zprvPause');

             % Set the button icon
            obj.hButtons.screenshot.CData = icon;
            obj.hButtons.pause.CData = icon;
            
            fprintf('gridDimension: %d, gridPanelSize: %d, rowCount: %d, colCount: %d\n',gridDimension,gridPanelSize,ceil(1/gridDimension),mod(0,gridDimension) + 1);            
            
            %TODO: Use gobjects
            for i=1:numPlots
                
                %Place panel from top-left, counting right and then down
                %rowCount = ceil(i/gridDimension);
                %colCount = mod(i-1,gridDimension) + 1;                
                %panelPosn = [(colCount-1)*gridPanelSize  (gridDimension-rowCount)*gridPanelSize gridPanelSize gridPanelSize];
                
                rowCount = ceil(i/gridDimensionX);
                colCount = mod(i-1,gridDimensionX) + 1;
                
                panelPosn = [(colCount-1)*gridPanelSizeX (gridDimensionY-rowCount)*gridPanelSizeY gridPanelSizeX gridPanelSizeY];
                
                obj.hPanels.waveform(i) = uipanel(obj.hFigs.waveform,'Position',panelPosn);
                obj.hPanels.raster(i) = uipanel(obj.hFigs.raster,'Position',panelPosn);
                obj.hPanels.psth(i) = uipanel(obj.hFigs.psth,'Position',panelPosn);
                
                %Places axes in panel and configure
                obj.hPlots(i) = axes('Parent',obj.hPanels.waveform(i),'Position',[0 0 1 1],'XLim',obj.horizontalRange); %,'YLim',obj.verticalRange);
                obj.hRasters(i) = axes('Parent',obj.hPanels.raster(i),'Position',[0 0 1 1],'XLim',obj.horizontalRangeRaster);
                obj.hPSTHs(i) = axes('Parent',obj.hPanels.psth(i),'Position',[0 0 1 1]);
                
                %TODO: Remove this function & just set here
                obj.zprvSetAxesProps(obj.hPlots(i));
                obj.zprvSetAxesProps(obj.hRasters(i));
                obj.zprvSetAxesProps(obj.hPSTHs(i));
                
            end
            
            obj.hChanLabels = struct('waveform',[],'raster',[],'psth',[]);
            
        end
    end
    
    %% PROPERTY ACCESS
    methods
        
        function set.dataReadMode(obj,val)
            obj.zprpAssertNotRunning('dataReadMode');
            obj.validatePropArg('dataReadMode',val);
            obj.dataReadMode = val;
        end
        
        function set.displayMode(obj,val)
            obj.zprpAssertNotRunning('displayMode');
            obj.validatePropArg('displayMode',val);
            
            %Impose requirements before allowing switch to 'raster' mode
            disallowRaster = false;
            if strcmpi(val,'raster')
                if isempty(obj.stimStartChannel)
                    fprintf(2,'WARNING: Must specify stimStartChannel in order to use spike raster display mode\n');
                    disallowRaster = true;
                end
                
                if ~isempty(obj.stimEventTypes) && isempty(obj.stimEventTypesDisplayed)
                    fprintf(2,'WARNING: A valid stimEventTypesDisplayed value must be specified in order to use spike raster display mode\n');
                    disallowRaster = true;
                end
                
                if disallowRaster
                    val = 'waveform';
                end
                
            end
            
            obj.displayMode = val;

            obj.zprvShowDisplayFig();
                       
            %Side-effects
            obj.zprvResetReducedData();
            obj.tabDisplayed = obj.tabDisplayed;
            
        end
        
        function set.filterWindow(obj,val)
            obj.validatePropArg('filterWindow',val);
            assert(val(2) > val(1),'Specified filter window must contain 2 elements in ascending order, specifying lower and upper frequency bound');
            
            obj.filterWindow = val;
            
            %Side-effects
            if isequal(val(:),[0;inf])
                obj.filterCoefficients = {};
            else
                if val(1) == 0
                    [b,a] = butter(1,obj.filterWindow(2) * 2 / obj.sglParamCache.niSampRate,'low');
                elseif isinf(val(2)) %high-pass filter
                    [b,a] = butter(1,obj.filterWindow(1) * 2 / obj.sglParamCache.niSampRate,'high');
                else
                    [b,a] = butter(1,obj.filterWindow * 2 / obj.sglParamCache.niSampRate,'bandpass');
                end
                
                obj.filterCoefficients = {a b};
            end
            
        end
        
        
        function set.globalMeanSubtraction(obj,val)
            obj.validatePropArg('globalMeanSubtraction',val);
            obj.globalMeanSubtraction = val;
        end
        
        function val = get.displayModeAxes(obj)
            switch obj.displayMode
                case 'raster'
                    val = obj.hRasters;
                case 'waveform'
                    val = obj.hPlots;
            end
        end
        
        function val = get.gridFigPosition(obj)
            val = get(obj.hFigs.(obj.displayMode),'Position');
        end
        
        function set.gridFigPosition(obj,val)
            set(obj.hFigs.(obj.displayMode),'Position',val);
        end
        
        function val = get.maxBufSizeScans(obj)
            val = round(obj.maxBufSizeSeconds * obj.sglParamCache.niSampRate);
        end
        
        
        function val = get.psthFigPosition(obj)
            val = get(obj.hFigs.psth,'Position');
        end
        
        function set.psthFigPosition(obj,val)
            set(obj.hFigs.psth,'Position',val);
        end
        
        function set.psthTimeBin(obj,val)
            obj.validatePropArg('psthTimeBin',val);
            obj.psthTimeBin = val;
        end
        
        function set.psthTimerPeriod(obj,val)
            obj.validatePropArg('psthTimerPeriod',val);
            
            changeWhileRunning = strcmpi(get(obj.hPSTHTimer,'Running'),'on');
            
            if changeWhileRunning
                stop(obj.hPSTHTimer);
            end
            
            if ~isinf(val)
                set(obj.hPSTHTimer,'Period',val,'StartDelay',val);
                
                if changeWhileRunning
                    start(obj.hPSTHTimer);
                end
            end
            
            obj.psthTimerPeriod = val;
        end
        
        function set.psthAmpRange(obj,val)
            obj.validatePropArg('psthAmpRange',val);
            
            set(obj.hPSTHs,'YLim',val);
            obj.psthAmpRange = get(obj.hPSTHs(1),'YLim');
        end
        
        function val = get.numTabs(obj)
            numNeuralChans = numel(obj.neuralChansAvailable);
            val = ceil(numNeuralChans/obj.PLOTS_PER_TAB);
        end
        
        function val = get.refreshPeriodMaxNumWaveforms(obj)
            val = obj.refreshPeriodMaxWaveformRate / obj.refreshRate;
        end
        
        function set.refreshPeriodMaxWaveformRate(obj,val)
            obj.validatePropArg('refreshPeriodMaxWaveformRate',val);
            
            lclVar = ceil(val / obj.refreshRate);
            obj.refreshPeriodMaxWaveformRate = lclVar * obj.refreshRate;
        end
        
        function val = get.refreshPeriodAvgScans(obj)
            val = round(get(obj.hTimer,'Period') * obj.sglParamCache.niSampRate);
        end
        %
        %         function val = get.neuralChanAcqList(obj)
        %             % obj.sglChanSubset does not discriminate between neural, aux, and other types of channels.
        %             % return obj.sglChanSubset's neural chans only. hGrid.sglParamCache.niMNChans1, hGrid.sglParamCache.niMNChans2
        %             lclMNChans = str2num(num2str(obj.sglParamCache.niMNChans1));
        %             lclChanLim = (max(lclMNChans) + 1) * obj.sglParamCache.niMuxFactor; % DO NOT HARDCODE THIS TO NUMBER OF WAVEFORMS PER TAB.
        %             %val = obj.sglChanSubset;
        %             val = obj.sglChanSubset(obj.sglChanSubset < lclChanLim);
        %         end
        
        function set.refreshRate(obj,val)
            obj.zprpAssertNotRunning('refreshRate');
            obj.validatePropArg('refreshRate',val);
            
            %Ensure value does not exceed the processing refresh period
            hrng = diff(obj.horizontalRange);
            f_samp = obj.sglParamCache.niSampRate;
            assert(ceil(hrng * f_samp) < floor(f_samp/val),'horizontalRange must be shorter than the processing refresh period');
            
            refreshPeriodRounded = round(1e3 * 1/val) * 1e-3; %Make an integer number of milliseconds
            set(obj.hTimer,'Period',refreshPeriodRounded);
            
            currMaxSpikeRate = obj.refreshPeriodMaxWaveformRate;
            
            obj.refreshRate = 1/get(obj.hTimer,'Period');
            
            %Side-effects
            obj.refreshPeriodMaxWaveformRate = currMaxSpikeRate;
        end
        
        function val = get.horizontalRangeScans(obj)
            val = round(obj.horizontalRange * obj.sglParamCache.niSampRate);
        end
        
        function val = get.stimEventTypes_(obj)
            if isempty(obj.stimEventTypes)
                val = {'allstim'};
            else
                val = obj.stimEventTypes;
            end
        end
        
        function val = get.stimEventClassifyChannel(obj)
            if isempty(obj.stimEventClassifyChannel)
                val = obj.stimStartChannel;
            else
                val = obj.stimEventClassifyChannel;
            end
        end
        
        function val = get.stimEventCount(obj)
            if isempty(obj.stimEventCount_) %Acq hasn't yet been started
                for i=1:length(obj.stimEventTypes_)
                    obj.stimEventCount_.(obj.stimEventTypes_{i}) = 0;
                end
            end
            
            if isempty(obj.stimEventTypes)
                val = obj.stimEventCount_.allstim;
            else
                val = sum(cellfun(@(eventType)obj.stimEventCount_.(eventType),obj.stimEventTypesDisplayed));
            end
        end
        
        function set.waveformAmpUnits(obj,val)
            obj.validatePropArg('waveformAmpUnits',val);

            % Save old value of vertical range to cache.
            oldVal = obj.waveformAmpUnits;          
            obj.verticalRangeCache.(obj.waveformAmpUnits) = obj.verticalRange;

            % Update the property value.
            obj.waveformAmpUnits = val;

            % Check to see if cached value exists in struct.
            if isfield(obj.verticalRangeCache,(val))
                %Restore cached value of vertical range if it exists.
                obj.verticalRange = obj.verticalRangeCache.(val);
            end            
        end
               
        function val = get.verticalRange(obj)
            val = get(obj.hPlots(1),'YLim');
        end
        
        function set.verticalRange(obj,val)
            obj.validatePropArg('verticalRange',val);
            
            if strcmpi(obj.waveformAmpUnits,'volts');
               aiRangeMax = obj.sglParamCache.niAiRangeMax;
               if any(abs(val) > 1.1 * aiRangeMax)
                   warning('Specified range exceeded input channel voltage range by greater than 10% -- spike amplitude axis limits clamped');
                   val = min(val,1.1 * aiRangeMax);
                   val = max(val,-1.1 * aiRangeMax);
               end
            end
            
            set(obj.hPlots,'YLim',val);
            
            %Set real property
            obj.verticalRange_ = val;
        end
        
        function set.verticalRange_(obj,val)
            %force recalc of dependent property
            obj.verticalRange_ = val;
        end
        
        function set.horizontalRange(obj,val)
            obj.zprpAssertNotRunning('horizontalRange');
            obj.validatePropArg('horizontalRange',val);
            
            %Ensure value does not exceed processing refresh period
            dval = diff(val);
            f_samp = obj.sglParamCache.niSampRate;
            
            fprintf('Horizontal range = ceil(%f * %f) = %f\n',dval,f_samp,ceil(dval * f_samp) );
            fprintf('Processing refresh period = floor(%f / %f) = %f\n',f_samp,obj.refreshRate,floor(f_samp/obj.refreshRate) );
            
            assert(dval > 0,'Horizontal range must be specified from minimum to maximum');
            assert(ceil(dval * f_samp) < floor(f_samp/obj.refreshRate),'horizontalRange must be shorter than the processing refresh period');
                        
            obj.horizontalRange = val;
            
            %Side-effects
            set(obj.hPlots,'XLim',val);
            
        end
        
        function val = get.maxPointsPerAnimatedLine(obj)
            %Calculate MaximumNumPoints
            
            if strcmp(obj.waveformsPerPlotClearMode,'oldest')
                spikeSampleRate = obj.sglParamCache.niSampRate;
                numPointsPerWindow = spikeSampleRate * (obj.horizontalRange(2)-obj.horizontalRange(1));
                val = ceil(obj.waveformsPerPlot * numPointsPerWindow);
            else
                val = Inf;
            end
        end
        
        %         function val = get.sglChanSubset(obj)
        %              val = GetSaveChansNi(obj.hSGL); %channel subset as specified in SpikeGLX. Wierd - this /has/ to be done here, outside of zprvZpplyChanOrderAndSubset() to avoid a hang.
        %         end
        
        
        % Why does setting refreshPeriodMaxWaveformRate > refreshRate cause
        % waveformsPerPlot to be bypassed? Even if I set the
        % waveformsPerPlot to 1, if refreshPeriodMaxWaveformRate = 16 and
        % refreshRate = 4, I get 4 lines per plot - and also my waveform
        % processing gets screwed up.
        
        % well, apparently, it has to do with some kind of max waveform
        % processing in zlcldetectspikes.
        
        
        function set.waveformsPerPlot(obj,val)
            obj.validatePropArg('waveformsPerPlot',val);
            obj.waveformsPerPlot = val;
        end
        
        function set.waveformsPerPlotClearMode(obj,val)
            obj.zprpAssertNotRunning('waveformsPerPlotClearMode');
            obj.validatePropArg('waveformsPerPlotClearMode',val);
            
            obj.waveformsPerPlotClearMode = val;
            
            %side-effects
            %TODO: add check here to reset max points in animatedline when changed to anything but 'oldest'.
            %TODO: add code to recompute max points on animatedlines when changed to 'oldest'.
            obj.zprvClearPlots('waveform');
        end
        
        function val = get.spikeRefractoryPeriod(obj)
            switch obj.displayMode
                case 'waveform'
                    val = obj.horizontalRange(2);
                case 'raster'
                    val = obj.spikeRefractoryPeriod;
            end
        end
        
        function set.spikeRefractoryPeriod(obj,val)
            obj.validatePropArg('spikeRefractoryPeriod',val);
            
            oldVal = obj.spikeRefractoryPeriod;
            switch obj.displayMode
                case 'waveform'
                    if obj.mdlInitialized && ~isequal(oldVal,val)
                        fprintf(2,'WARNING: When displayMode = ''waveform'', the ''spikeRefractoryPeriod'' cannot be directly set. Value specified ignored.\n')
                    end
                case 'raster'
                    obj.spikeRefractoryPeriod = val;
            end
        end
        
        %     function val = get.stimEventTypes(obj)
        %       if isempty(obj.stimEventClassifyFcn)
        %         val = '';
        %       else
        %         val = feval(obj.stimEventClassifyFcn);
        %       end
        %     end
        
        function set.stimEventTypesDisplayed(obj,val)
            assert(ischar(val) || iscellstr(val),'Value of ''stimEventTypesDisplayed'' must be either a string or string cell array');
            
            if isempty(obj.stimEventTypes)
                assert(isempty(val),'Value of ''stimEventTypesDisplayed'' must be empty when ''stimEventTypes'' is empty');
                return;
            else
                assert(~isempty(val),'Valid value of ''stimEventTypesDisplayed'' must be supplied.');
            end
            
            if ~iscell(val)
                val = {val};
            end
            
            assert(all(ismember(val,obj.stimEventTypes)),'One or more of the specified stim event types not recognized');
            
            oldVal = obj.stimEventTypesDisplayed;
            obj.stimEventTypesDisplayed = val;
            
            %Side-effects
            if ~isequal(oldVal,val) && strcmpi(obj.displayMode,'raster')
                obj.zprvClearPlots('raster');
                if obj.running
                    obj.zprvUpdateRasterPlot();
                end
            end
            
        end
        
        function set.stimEventClassifyFcn(obj,val)
            obj.zprpAssertNotRunning('stimEventClassifyFcn');
            
            %Custom validation
            isFcn = isscalar(val) && isa(val,'function_handle');
            assert(isempty(val) || isFcn, 'Property ''stimEventClassifyFcn'' must be a function handle (or empty)');
            
            errorString = '';
            if isFcn
                try
                    s = feval(val);
                    obj.stimEventTypes = s.stimEventTypes;
                    obj.stimEventClassifyNumScans = s.stimEventClassifyNumScans;
                    
                    if isfield(s,'stimEventClassifyChannel')
                        obj.stimEventClassifyChannel = s.stimEventClassifyChannel;
                    else
                        obj.stimEventClassifyChannel = [];
                    end
                    
                    if isempty(obj.stimEventTypes)
                        obj.stimEventTypesDisplayed = '';
                    elseif isempty(obj.stimEventTypesDisplayed) || ~all(ismember(obj.stimEventTypesDisplayed,obj.stimEventTypes))
                        obj.stimEventTypesDisplayed = obj.stimEventTypes{1};
                    end
                    
                    obj.stimEventCount_ = struct();
                    for i=1:length(obj.stimEventTypes)
                        obj.stimEventCount_.(obj.stimEventTypes{i}) = 0;
                    end
                    
                catch ME
                    errorString = 'Specified ''stimEventClassifyFcn'' must return structure with fields ''stimEventTypes'' and ''stimEventClassifyNumScans'' when called with no arguments.';
                end
            end
            
            if ~isFcn || ~isempty(errorString)
                val = [];
                obj.stimEventTypes = '';
                obj.stimEventClassifyNumScans = 0;
                
                obj.stimEventCount_ = struct();
                obj.stimEventCount_.allstim = 0;
            end
            
            if errorString
                error(errorString);
            else
                obj.stimEventClassifyFcn = val;
            end
            
            %Side-effects
            obj.zprvInitializeRasterGridLines(); %Initialize raster grid animated line objects
            
            
        end
        
        function set.verticalRangeRaster(obj,val)
            obj.validatePropArg('verticalRangeRaster',val);
            
            ylim = obj.zprvverticalRangeRaster2YLim(val);
            
            set(obj.hRasters,'YLim',ylim);
            obj.verticalRangeRaster = val;
        end
        
        function set.verticalRangeRasterInfIncrement(obj,val)
            obj.validatePropArg('verticalRangeRasterInfIncrement',val);
            obj.verticalRangeRasterInfIncrement = val;
        end
        
        function set.stimStartChannel(obj,val)
            obj.zprpAssertNotRunning('stimStartChannel');
            obj.validatePropArg('stimStartChannel',val);
            obj.zprpAssertAuxChan(val,'stimStartChannel'); %Assert any stimulus channel specified is a valid auxiliary channel
            obj.stimStartChannel = val;
        end
        
        function set.stimStartThreshold(obj,val)
            obj.validatePropArg('stimStartThreshold',val);
            obj.stimStartThreshold = val;
        end
        
        function set.horizontalRangeRaster(obj,val)
            obj.zprpAssertNotRunning('horizontalRangeRaster');
            obj.validatePropArg('horizontalRangeRaster',val);
            
            obj.zprvClearPlots('raster');
            set(obj.hRasters,'XLim',val);
            
            obj.zprvClearPlots('psth');
            
            obj.horizontalRangeRaster = val;
        end
        
        
        function set.tabDisplayed(obj,val)
            obj.validatePropArg('tabDisplayed',val);
            assert(val <= obj.numTabs,'Value specified (%d) exceeds the number of available tabs (%d)',val,obj.numTabs);
            
            obj.zprvAssertAvailChansConstant();
            
            obj.blockTimer = true;
            
            try
                obj.tabDisplayed = val;
                dispType = obj.displayMode;
                
                %Update tabChanNumbers
                numNeuralChans = numel(obj.neuralChansAvailable); %#ok<*MCSUP>
                
                tcn = (1:obj.PLOTS_PER_TAB) + (val-1)*obj.PLOTS_PER_TAB;
                tcn(tcn > numNeuralChans) = [];
                obj.tabChanNumbers = tcn; %One-based channel index
                
                %Clear grid waveform/raster plots
                obj.zprvClearPlots({obj.displayMode 'psth'},false); %Don't reuse existing threshold lines
                
                %Update plot channel labels
                most.idioms.deleteHandle(obj.hChanLabels.(dispType));
                most.idioms.deleteHandle(obj.hChanLabels.psth);
                
                hAxes = obj.displayModeAxes;

                for i=1:length(tcn)
                    actualChannelNumber = tcn(i) - 1;
                    if ~isequal(obj.neuralChanDispOrder, obj.neuralChansAvailable)
                        %Use channel mapping numbers (if they exist)
                        tempChannelNumber = obj.neuralChanDispOrder( tcn(i) );
                        if tempChannelNumber ~= tcn(i) - 1
                            channelNumberString = strcat(num2str(tempChannelNumber),' (',num2str(tcn(i) - 1),')');
                            textPosition = [.13 .92];
                        else
                            channelNumberString = num2str(tcn(i) - 1);
                            textPosition = [.08 .92];
                        end
                    else
                        %Display with 0-based channel index
                        tempChannelNumber = tcn(i) - 1;
                        channelNumberString = num2str(tcn(i) - 1);
                        textPosition = [.08 .92];
                    end
                    
                    obj.hChanLabels.(dispType)(i) = text('Parent',hAxes(i),'String',channelNumberString,'HandleVisibility','off','FontWeight','bold','Color','k','Units','normalized','Position',textPosition,'HorizontalAlignment','center');
                    obj.hChanLabels.psth(i) = text('Parent',obj.hPSTHs(i),'String',channelNumberString,'HandleVisibility','off','FontWeight','bold','Color','k','Units','normalized','Position',textPosition,'HorizontalAlignment','center');
                end
                
                %Update threshold lines/raster plots, as appropriate
                if strcmpi(obj.displayMode,'raster')
                    obj.zprvUpdateRasterPlot();
                end
                
                drawnow expose update;
                
                obj.blockTimer = false;
                
            catch ME
                obj.blockTimer = false;
                ME.rethrow();
            end
        end
        
        function set.thresholdAbsolute(obj,val)
            obj.validatePropArg('thresholdAbsolute',val);
            obj.thresholdAbsolute = val;
            
            %Side-effects
            obj.zprvDrawThresholdLines();
        end
        
        function val = get.baselineStatsRefreshPeriodScans(obj)
            val = round(obj.baselineStatsRefreshPeriod * obj.sglParamCache.niSampRate);
        end
        
        function set.thresholdType(obj,val)
            obj.zprpAssertNotRunning('thresholdType');
            obj.validatePropArg('thresholdType',val);
            
            oldVal = obj.thresholdType;
            obj.thresholdType = val;
            
            %Side Effects
            if ~strcmpi(oldVal,val)
                
                %Adjust thresholdVal & verticalRange
                
                %TODO(?): A smarter adjustment based on the last-cached RMS values, somehow handlign the variety across channels
                switch val
                    case 'volts'
                        aiRangeMax = obj.sglParamCache.niAiRangeMax;
                        obj.thresholdVal = .1 * aiRangeMax;
                    case 'rmsMultiple'
                        obj.thresholdVal = 5;
                end
                
                %Redraw threshold lines
                obj.zprvDrawThresholdLines();
            end
        end
        
        function set.thresholdVal(obj,val)
            %obj.zprpAssertNotRunning('thresholdVal');
            obj.validatePropArg('thresholdVal',val);
            
            obj.thresholdVal = val;
            
            %Side-effects
            obj.zprvDrawThresholdLines();
        end
        
        function set.baselineStatsRefreshPeriod(obj,val)
            obj.zprpAssertNotRunning('baselineStatsRefreshPeriod');
            obj.validatePropArg('baselineStatsRefreshPeriod',val);
            
            obj.baselineStatsRefreshPeriod = val;
        end
        
    end
    
    
    %Property-access helper functions
    methods (Access=protected)
        
        function zprpAssertNotRunning(obj,propName)
            assert(~obj.running,'Property ''%s'' cannot be set while grid waveform display is running',propName);
        end
        
        function zprpAssertAuxChan(obj,chanNum,propName)
            assert(isempty(chanNum) || ismember(chanNum,union(obj.analogMuxChansAvailable, obj.analogSoloChansAvailable)), 'The property ''%s'' must specify one of the available auxiliary channels',propName);
            %assert(isempty(chanNum) || (chanNum <= (obj.numChansTotal-1) && chanNum >= (obj.numChansTotal - obj.numAuxChans)),'The property ''%s'' must specify one of the %d auxiliary channels',propName,obj.numAuxChans);
        end
        
        
    end
    
    %% PUBLIC METHODS
    methods
        function start(obj)
            
            if obj.running
                return;
            end
            
            if obj.debug
                obj.diagramm = figure('name','raw spike plot (set to debug to false to get rid of this.');
                %obj.diagrammm = figure('name','temporary plot 2');
            end
            
            %Open SpikeGL connection & updateparameter cache
            obj.hSGL = SpikeGL(obj.sglIPAddress);
            obj.sglParamCache = GetParams(obj.hSGL);
            
            %Initializations
            obj.maxReadableScanNum = 0;
            obj.lastMaxReadableScanNum = 0;
            obj.priorfileMaxReadableScanNum = 0;
            obj.waveformWrap = [];
            obj.partialWaveformBuffer = {};
            
            hTimers = obj.hTimer;
            
            %Apply channel ordering & subsetting, if specified in SpikeGLX
            %TODO: Consider to allow further subsetting by Spoke user, to give faster Spoke processing
            obj.zprvAssertAvailChansConstant();
            obj.zprvApplyChanOrderAndSubset();
            
            %Reset various state vars -- RMS/mean, filterCondition, spike data, etc
            obj.zprvResetAcquisition();
            
            %Clear previous lines
            handlesToClear = [];
            switch obj.displayMode
                case 'waveform'
                    obj.zprvClearPlots('waveform');
                    %                     for i=1:numDispChans
                    %                         handlesToClear = [handlesToClear; obj.hWaveforms{i}(isgraphics(obj.hWaveforms{i}))'];
                    %                     end
                    %                     %VI051012: Seems like we should probably clear out obj.hWaveforms here too
                    %                     %set(handlesToClear,'EraseMode','normal');
                    %                     delete(handlesToClear);
                    
                case 'raster'
                    obj.zprvClearPlots({'raster' 'psth'});
            end
            
            %Show display fig
            obj.zprvShowDisplayFig();
            
            %Display-type specific initialization
            switch obj.displayMode
                case 'waveform'
                    
                    %Draw new threshold lines, if applicable
                    obj.zprvDrawThresholdLines();
                    
                    [obj.stimScanNums, obj.stimWindowStartScanNums, obj.stimWindowEndScanNums] = deal([]);
                    obj.stimEventTypeNames = {};
                    
                case 'raster'
                    
                    obj.stimTotalCount = 0;
                    for i=1:length(obj.stimEventTypes_)
                        obj.stimEventCount_.(obj.stimEventTypes_{i}) = 0;
                    end
                    
                    [obj.stimScanNums, obj.stimWindowStartScanNums, obj.stimWindowEndScanNums] = deal([]);
                    obj.stimEventTypeNames = {};
                    
                    %Initialize reducedData stimEventTypeStruct
                    
                    if isempty(obj.stimEventTypes)
                        fnames = {'allstim'};
                    else
                        fnames = obj.stimEventTypes;
                    end
                    
                    for i=1:numel(obj.neuralChansAvailable)
                        for j=1:length(fnames)
                            obj.reducedData{i}.stimEventTypeStruct.(fnames{j}) = [];
                        end
                    end
                    
                    obj.zprvResetStimNumsPlotted();
                    
                    if ~isinf(obj.psthTimerPeriod)
                        hTimers(end+1) = obj.hPSTHTimer;
                    end
                    
            end
            
            %Apply graphics updates
            drawnow expose update;
            
            %Reset stats
            obj.bmarkReadTimeStats = zeros(3,1);
            obj.bmarkPreProcessTimeStats = zeros(3,1);
            obj.bmarkPlotTimeStats = zeros(3,1);
            obj.bmarkPostProcessTimeStats = zeros(3,1);
            
            %Start timer(s)
            start(hTimers);
            
            obj.running = true;
            
            disp('started');
            
        end
        
        function stop(obj)
            if ~obj.running
                return;
            end
            
            stop([obj.hTimer obj.hPSTHTimer]);
            if ~isempty(obj.hSGL)
                ME = [];
                try
                    Close(obj.hSGL);
                catch MEtemp
                    ME = MEtemp;
                end
                
                obj.hSGL = [];
                
                if ~isempty(ME)
                    ME.rethrow();
                end
                
            end
            
            obj.running = false;
            
            %Reset to initial conditions
            obj.maxNumWaveformsApplied = false;
        end
        
        function loadConfig(obj,filename)
            
            if nargin < 2
                if isempty(obj.configFileName)
                    startName = obj.getClassDataVar('lastConfigFileName');
                else
                    startName = obj.configFileName;
                end
                [f,p] = uigetfile('*.cfg','Load Spoke configuration file',startName);
                if isnumeric(f)
                    return;
                end
                filename = fullfile(p,f);
            else
                errorMessage = 'Argument ''filename'' supplied was not a valid string or did not specify an existing configuration (CFG) file';
                assert(most.idioms.isstring(filename),errorMessage);
                assert(exist(filename,'file'),errorMessage)
                [~,~,e] = fileparts(filename);
                assert(strcmpi(e,'.cfg'),errorMessage);
            end
            
            obj.mdlLoadConfig(filename);
            
            if ~isempty(obj.hController)
                assert(isscalar(obj.hController));
                obj.hController{1}.ctlrLoadGUILayout(filename);
            end
            
            %Cache name of last-saved/loaded config filename
            obj.configFileName = filename;
            obj.setClassDataVar('lastConfigFileName',obj.configFileName);
            
        end
        
        function saveConfigAs(obj,filename)
            
            if nargin < 2
                if isempty(obj.configFileName)
                    startName = obj.getClassDataVar('lastConfigFileName');
                else
                    startName = obj.configFileName;
                end
                [f,p] = uiputfile('*.cfg','Save Spoke configuration file',startName);
                if isnumeric(f)
                    return;
                end
                filename = fullfile(p,f);
            else
                errorMessage = 'Argument ''filename'' supplied was not a valid string or did not specify a valid configuration (CFG) file';
                assert(most.idioms.isstring(filename),errorMessage);
                [p,~,e] = fileparts(filename);
                assert(exist(p,'dir'),errorMessage);
                assert(strcmpi(e,'.cfg'),errorMessage);
            end
            
            %Save model properties
            obj.mdlSaveConfig(filename,'include',{'verticalRange' 'gridFigPosition' 'psthFigPosition'});
            
            %Save controller fig layout
            if ~isempty(obj.hController)
                assert(isscalar(obj.hController));
                obj.hController{1}.ctlrSaveGUILayout(filename);
            end
            
            %Cache name of last-saved/loaded config filename
            obj.configFileName = filename;
            obj.setClassDataVar('lastConfigFileName',obj.configFileName);
            
        end
        
        function plotPSTH(obj)
            if (~strcmpi(obj.displayMode,'raster'))
                %assert(strcmpi(obj.displayMode,'raster'),'Plot PSTH operation only available when displayMode=''raster''');
                fprintf(2,'WARNING: Plot PSTH operation only available when displayMode=''raster''\n');
                return;
            end
            
            obj.blockTimer = true;
            
            ME = [];
            
            try
                obj.zprvClearPlots('psth');
                
                set(obj.hFigs.psth,'Visible','on');
                colorOrder = get(0,'DefaultAxesColorOrder');
                
                %VVV: use a hidden prop instead of this logic, here and elsewhere
                if isempty(obj.stimEventTypes)
                    eventTypes = {'allstim'};
                else
                    eventTypes = obj.stimEventTypesDisplayed;
                end
                
                if length(eventTypes) < length(obj.stimEventTypes)
                    [tf,colorIdxs] = ismember(eventTypes,obj.stimEventTypes);
                    assert(all(tf));
                else
                    colorIdxs = 1:length(eventTypes);
                end
                
                scanPeriod = 1/ obj.sglParamCache.niSampRate;
                scansToBin = max(1,round(obj.psthTimeBin/scanPeriod));
                scanPeriodBinned = scanPeriod * scansToBin;
                
                histogramBins = floor(obj.horizontalRangeRaster(1)/scanPeriod):scansToBin:ceil(obj.horizontalRangeRaster(2)/scanPeriod);
                
                plotCount = 0;
                for c=obj.tabChanNumbers
                    plotCount = plotCount + 1;
                    
                    if isempty(obj.reducedData{c}.scanNums)
                        continue;
                    end
                    
                    for e=1:length(eventTypes)
                        eventType = eventTypes{e};
                        
                        eventIdxs = obj.reducedData{c}.stimEventTypeStruct.(eventType);
                        
                        histData = hist(obj.reducedData{c}.stimRelScanNums(eventIdxs),histogramBins);
                        histIdxs = (histData > 0);
                        
                        axes(obj.hPSTHs(plotCount));
                        if ~isempty(histIdxs)
                            line('XData',histogramBins(histIdxs),'YData',histData(histIdxs)/(obj.stimEventCount_.(eventType)*scanPeriodBinned),'LineStyle','-','Color',colorOrder(colorIdxs(e),:));
                        end
                        
                    end
                end
                
            catch ME_
                ME = ME_;
            end
            
            %Clean-up
            obj.blockTimer = false;
            
            if ~isempty(ME)
                ME.rethrow();
            end
        end
        
        
        function sortRastergram(obj)
            %Method to sort stimuli displayed in rastergram by stimulus event type
            
            assert(iscell(obj.stimEventTypesDisplayed) && length(obj.stimEventTypesDisplayed) > 1,'Only possible (or useful) to reorganize rastergram when stimEventTypesDisplayed specifies more than one stimulus event type');
            
            obj.blockTimer = true;
            
            ME = [];
            
            
            try
                obj.zprvClearPlots('raster');
                obj.zprvUpdateRasterPlot();
            catch ME_
                ME = ME_;
            end
            
            %Clean-up
            obj.blockTimer = false;
            
            if ~isempty(ME)
                ME.rethrow();
            end
            
            
        end

        function outputChanArray = parseChanMapFile(obj,filename)
            % obj.sglParamCache.snsNiChanMapFile
            fileId = fopen(filename,'r');
            tempString = fread(fileId,'*char');
            
            % Assume that the file orders channels 0 to 255. So to create a
            % mapping, just read in the second channel to the mapping
            % array.
           
            tempInputArray = [];
            inputChanArray = [];

            tempOutputArray = [];
            outputChanArray = [];
            
            outputChan = false;
            inputChan = false;
            processChan=false;
            
            for iter = 1:size(tempString,1)
                % Gate whether we are looking at input channels or output
                % channels.
                
                if tempString(iter) == ' '
                   outputChan = true;
                   inputChan = false;
                   processChan=false;
                elseif tempString(iter) == ';'
                   outputChan = false;
                   inputChan = true;
                   processChan=false;
                elseif tempString(iter) == 'M' || ismember(tempString(iter), char([10 13]))
                   outputChan = false;
                   inputChan = false;
                   processChan = true;
                end
                                
                % Set up the arrays as necessary.
                if outputChan                    
                    if isstrprop(tempString(iter),'digit') % tempString(iter) ~= ' ' && tempString(iter) ~= ';' && ~isletter(tempString(iter))
                        tempOutputArray = [tempOutputArray tempString(iter)];
                    end
                    
                    %Convert the temp arrays to integers.
                    if ~isempty(tempInputArray) %&& ismember(str2num(tempInputArray), obj.neuralChansAvailable)
                        if ismember(str2num(tempInputArray), [obj.neuralChansAvailable, obj.analogMuxChansAvailable ])
                            inputChanArray = [inputChanArray str2num(tempInputArray)];
                            tempInputArray = [];
                        end
                    end                    
                elseif inputChan
                    if isstrprop(tempString(iter),'digit') % tempString(iter) ~= ' ' && tempString(iter) ~= ';' && ~isletter(tempString(iter))
                        tempInputArray = [tempInputArray tempString(iter)];
                    end
                                    
                    if ~isempty(tempOutputArray) %&& ismember(str2num(tempOutputArray), obj.neuralChansAvailable)
                        if ismember(str2num(tempOutputArray), [obj.neuralChansAvailable, obj.analogMuxChansAvailable ])
                            outputChanArray = [outputChanArray str2num(tempOutputArray)];
                            tempOutputArray = [];
                        end
                    end
                elseif processChan
                end
            end
        end
        
        function quit(obj)
            % Delete timer objects.
            if isvalid(obj.hTimer)
                stop(obj.hTimer)
                delete(obj.hTimer)
            end
            
            if isvalid(obj.hPSTHTimer)
                stop(obj.hPSTHTimer)
                delete(obj.hPSTHTimer)
            end
            
            % Delete figures.
            delete(obj.hFigs.waveform)
            delete(obj.hFigs.raster)
            delete(obj.hFigs.psth)
            delete(obj.hPanels.waveform)
            delete(obj.hPanels.raster)
            delete(obj.hPanels.psth)
            %            delete(obj.hPlots)
            %     delete(obj.hRasters)
            %  delete(obj.hPSTHs)
        end
    end
    
    %% HIDDEN METHODS
    methods (Hidden)
        
        function zcbkTimerFcn(obj,src,evnt)
            if obj.blockTimer
                return;
            end
            
            try
                
                newWaveformWrapVal = [];
                                
                t0 = tic;
                
                % cnt = GetScanCount(obj.hSGL);
                
                cnt = GetScanCountNi(obj.hSGL);
                
                %Use current scan number as reference scan number on first timer entry following start/restart
                if ismember(obj.bufScanNumEnd,[0 -1])
                    %fprintf('bufScanNumEnd: %d\n',obj.bufScanNumEnd);
                    obj.bufScanNumEnd = cnt;
                    obj.lastMaxReadableScanNum = cnt;
                    return;
                else
                    obj.maxReadableScanNum = cnt;
                end
                
                numNeuralChans = numel(obj.neuralChansAvailable);
                
                rmsMultipleThresh = strcmpi(obj.thresholdType,'rmsMultiple');
                rmsMultipleInitializing = rmsMultipleThresh && obj.baselineRMSLastScan == 0;
                
                waveformDisplay = strcmpi(obj.displayMode,'waveform'); %either in spike-triggered or stimulus-triggered waveform display mode
                rasterDisplayMode = strcmpi(obj.displayMode,'raster');
                stimulusMode = ~isempty(obj.stimStartChannel) && ~isempty(obj.stimStartThreshold);
                stimulusTriggeredWaveformMode = strcmpi(obj.displayMode,'waveform') && stimulusMode;
                
                sampRate = obj.sglParamCache.niSampRate;
                sampPeriod = 1 / sampRate;
                
                %Handle case of no new data
                if obj.maxReadableScanNum == obj.lastMaxReadableScanNum %no new data
                    
                    %Check if SpikeGL has stopped
                    if ~IsRunning(obj.hSGL)
                        obj.stop();
                        return;
                    else % Logging has been stopped...may resume anytime
                        return;
                    end
                end
                
                %Previous approach to detecting SpikeGL stop/restart - now
                %defunct. Superceded by changedFilename test.
                if obj.maxReadableScanNum < obj.lastMaxReadableScanNum
                    assert(false);
                end
                
                %dfprintf('obj.maxReadableScanNum: %d     obj.lastMaxReadableScanNum: %d\n',obj.maxReadableScanNum,obj.lastMaxReadableScanNum);
                scansToRead_ = obj.maxReadableScanNum - obj.lastMaxReadableScanNum;
                
                meanScansPerTimerTick = sampRate / obj.refreshRate;
                scansToRead = min(scansToRead_, round(2  * meanScansPerTimerTick));
                
                if scansToRead < scansToRead_
                    fprintf(2,'WARNING. A large number of queued-up samples to read detected: %d. If intermittent, this should not cause a problem.\n',scansToRead_ - scansToRead);
                end
                
                %STAGE 1: Read newly available data
                %[scansToRead, newData] = znstReadAvailableData(obj.maxReadableScanNum-obj.priorfileMaxReadableScanNum-scansToRead,scansToRead); %obj.bufScanNumEnd will be 0 in case of file-rollover or SpikeGL stop/restart
                %                 newData = GetDAQData(obj.hSGL,obj.lastMaxReadableScanNum,scansToRead,obj.sglChanSubset);
                newData = FetchNi(obj.hSGL,obj.lastMaxReadableScanNum,scansToRead,obj.sglChanSubset);
                obj.lastMaxReadableScanNum = obj.lastMaxReadableScanNum + scansToRead;
                t1 = toc(t0);
                
                
                %STAGE 2: Apply global mean subtraction, if applicable. Applies only to neural channels. TODO: Restrict to displayed channels
                if obj.globalMeanSubtraction
                    %                     newData(:,1:numNeuralChans) = newData(:,1:numNeuralChans) - mean(mean(newData(:,1:numNeuralChans)));
                    newData(:,1:obj.sglChanSubset) = newData(:,1:obj.sglChanSubset) - mean(mean(newData(:,1:obj.sglChanSubset)));
                end
                t2 = toc(t0);
                
                %STAGE 3: Filter data if needed. Also: housekeeping to form fullDataBuffer & partialWaveformBuffer for subsequent processing stages.
                if ~isempty(obj.filterCoefficients)
                    [newData,obj.filterCondition] = filter(obj.filterCoefficients{2},obj.filterCoefficients{1},double(newData),obj.filterCondition); %Convert to double..but still in A/D count values, not voltages
                end
                
                %Housekeeping: Form partialWaveformBuffer, appending new data
                if ~isempty(obj.waveformWrap)
                    for hiter=1:numel(obj.neuralChanAcqList)
                        iter = obj.sglChanSubset(hiter)+1;
                        obj.partialWaveformBuffer{end, iter} = [obj.partialWaveformBuffer{end, iter};newData(1:obj.waveformWrap(end),iter)];
                    end
                end
                
                %Housekeeping: Form fullDataBuffer, appending new data
                
                %fprintf('size of full data buffer before: %s\n',sprintf('%d ', size(obj.fullDataBuffer)));
                bufStartScanNum = znstAugmentFullDataBuffer(scansToRead, newData);
                %fprintf('size of full data buffer after: %s\n',sprintf('%d ', size(obj.fullDataBuffer)));
                
                %Not sure what actual edge case this checks for. Let's warn
                %instead of returning preemptively (which causes a
                %guaranteed error). TODO: review what cases this might pertain to.
                if size(obj.fullDataBuffer,1) < (obj.horizontalRangeScans(2) + 2)
                    %return; %warn instead of returning preemptively
                    if ~rasterDisplayMode
                        warning('Size fullDataBuffer is small relative to horizontalRange. This is suspected to potentially cause an issue/error.\n');
                    end
                end
                
                t3 = toc(t0);
                
                %STAGE 4: Detect stimulus and/or spike(s) within data buffer
                newSpikeScanNums = cell(numNeuralChans,1);
                
                %STAGE 4a: Detect spike(s) within data buffer
                %  (Modes: Spike-triggered Waveform & Raster)
                %  Excludes spikes during:
                %   * 'refractory period'(discarded)
                %   * final 'post-trigger' time, as spec'd by horizontalRange  (processed during next timer period)
                
                if ~stimulusTriggeredWaveformMode
                    if rmsMultipleInitializing %Handle case where no RMS data has been computed yet
                        %obj.fullDataBuffer((obj.refreshPeriodAvgScans+1):end,:) = []; %VVV062812: Is this needed/wanted?
                        
                        znstUpdateBaselineStats([],[]); %Compute rms/mean without spikes
                        newSpikeScanNums = zprvDetectNewSpikes(obj,bufStartScanNum); %Detect spikes using rmsMultiple=obj.INIT_RMS_THRESHOLD
                        znstUpdateBaselineStats(newSpikeScanNums,bufStartScanNum); %Recompute a rms value with, if anything, excess spike exclusion
                        
                        %newSpikeScanNums = cell(numDispChans,1); %Don't plot/store these spikes, though
                        newSpikeScanNums = zprvDetectNewSpikes(obj,bufStartScanNum);
                    else
                        newSpikeScanNums = zprvDetectNewSpikes(obj,bufStartScanNum);
                    end
                end
                
                %STAGE 4b: Detect stimulus within data buffer
                % (modes: Stimulus-triggered Waveform & Raster)
                % Currently only detects up to one stimulus per data period
                if stimulusTriggeredWaveformMode
                    obj.stimScanNums=[];
                    znstDetectStimulus(bufStartScanNum);
                elseif rasterDisplayMode
                    znstDetectStimulus(bufStartScanNum);
                end                
                
                t4 = toc(t0);
                
                %STAGE 5: Data reduction of curent full data batch, storing only needed features (waveform snippets, timestamps, spike-stimulus associations) for subsequent stages & timer periods
                                
                %STAGE 5a: (Waveform modes only) Determine and store timestamps & waveform snippets                
                % Raster & Spike-triggered Waveform modes: store detected spike scan numbers
                % Stimulus-triggered Waveform mode: store detected stimulus scan numbers
                % Waveform modes: store waveform snippets
                % also: initialize stimulus spike-tagging storage for Raster mode (TODO:  Review if STAGE 6 could be combined with STAGE 4, updating this step accordingly)
                if ~stimulusTriggeredWaveformMode
                    timestampOffsets = newSpikeScanNums;
                else
                    timestampOffsets = obj.stimScanNums;
                end
                
                znstStoreReducedData(timestampOffsets,bufStartScanNum);
                                
                %STAGE 5b: (Raster mode only) Determine & store spike-stimulus associations, and stimulus classifications if applicable 
                if rasterDisplayMode
                    if ~isempty(obj.stimScanNums)
                        znstClassifyStimulus(bufStartScanNum); %Classify stimulus event type, if possible
                    end
                    chanNewSpikes = znstAssociateSpikesToStimuli(); %Tag spike data with stimulus/event info, as needed/possible
                end
                t5 = toc(t0);
                                                
                %STAGE 6: Plot newly detected spike(s) that were stored for display - will always be enough post data, and generally enough pre-data except for spikes at very beginning
                if rasterDisplayMode
                    obj.zprvUpdateRasterPlot(chanNewSpikes);
                elseif stimulusTriggeredWaveformMode
                    obj.zprvUpdateWaveformPlot();
                else
                    obj.zprvUpdateWaveformPlot();
                end
                                
                %Housekeeping: Update waveformWrap val, for subsequent timer period
                if ~isempty(newWaveformWrapVal)
                    obj.waveformWrap(end+1) = newWaveformWrapVal;
                end
                
                t6 = toc(t0);               
                
                %STAGE 7: Update current baseline stats values (mean & RMS), if needed
                if (rmsMultipleThresh || obj.filterWindow(1) == 0) && ...
                        (obj.bufScanNumEnd - obj.baselineRMSLastScan) > obj.baselineStatsRefreshPeriodScans % enough time has elapsed since last RMS sampling
                    znstUpdateBaselineStats(newSpikeScanNums,bufStartScanNum);
                end
                t7 = toc(t0);
                
                %Housekeeping: reset the fullDataBuffer
                %TODO: Review moving this to fullDataBuffer formation step above in STAGE 3
                try
                    if rasterDisplayMode
                        %Raster mode: leave all but the number of scans required to classify
                        obj.fullDataBuffer(1:end-max(1,obj.stimEventClassifyNumScans)+1,:) = [];
                    else
                        %Waveform mode: leave only one full horizontalRange (pre+post+1 sample) at the end
                       % fprintf('obj.horizontalRangeScans: %s\n',sprintf('%d ', obj.horizontalRangeScans));
                        obj.fullDataBuffer(1:end-(diff(obj.horizontalRangeScans)+1),:) = [];
                    end
                catch ME
                    fprintf(2,'Error contracting the data buffer, which has size: %s\n',mat2str(size(obj.fullDataBuffer)));
                    ME.rethrow();
                end
                t8 = toc(t0);
                
                %fprintf('Total Time (%d scans of %d channels): %g\tRead: %g\tMeanSubtract: %g\tFilter: %g\tDetect: %g\tStore: %g\tStimTag: %g\tPlot: %g \t Mean Compute: %g\tContraction: %g\n',scansToRead,size(newData,2),t8*1000,t1*1000,(t2-t1)*1000,(t3-t2)*1000,(t4-t3)*1000,(t5-t4)*1000,(t6-t5)*1000,(t7-t6)*1000,(t8-t7)*1000,(t9-t8)*1000);
                readTime = 1000 * t1;
                procTimePre = 1000 * (t5-t1);
                plotTime = 1000 * (t6-t5);
                procTimePost = 1000 * (t8-t6);
                
                n = obj.bmarkReadTimeStats(3);
                mu = (obj.bmarkReadTimeStats(1) * n + readTime) / (n+1);
                std = sqrt((obj.bmarkReadTimeStats(2)^2 * n + (readTime - mu)^2)/ (n+1));
                obj.bmarkReadTimeStats(1) = mu;
                obj.bmarkReadTimeStats(2) = std;
                obj.bmarkReadTimeStats(3) = obj.bmarkReadTimeStats(3) + 1;
                
                n = obj.bmarkPreProcessTimeStats(3);
                mu = (obj.bmarkPreProcessTimeStats(1) * n + procTimePre) / (n+1);
                std = sqrt((obj.bmarkPreProcessTimeStats(2)^2 * n + (procTimePre - mu)^2)/ (n+1));
                obj.bmarkPreProcessTimeStats(1) = mu;
                obj.bmarkPreProcessTimeStats(2) = std;
                obj.bmarkPreProcessTimeStats(3) = obj.bmarkPreProcessTimeStats(3) + 1;
                
                n = obj.bmarkPlotTimeStats(3);
                mu = (obj.bmarkPlotTimeStats(1) * n + plotTime) / (n+1);
                std = sqrt((obj.bmarkPlotTimeStats(2)^2 * n + (plotTime - mu)^2)/ (n+1));
                obj.bmarkPlotTimeStats(1) = mu;
                obj.bmarkPlotTimeStats(2) = std;
                obj.bmarkPlotTimeStats(3) = obj.bmarkPlotTimeStats(3) + 1;
                
                n = obj.bmarkPostProcessTimeStats(3);
                mu = (obj.bmarkPostProcessTimeStats(1) * n + procTimePost) / (n+1);
                std = sqrt((obj.bmarkPostProcessTimeStats(2)^2 * n + (procTimePost - mu)^2)/ (n+1));
                obj.bmarkPostProcessTimeStats(1) = mu;
                obj.bmarkPostProcessTimeStats(2) = std;
                obj.bmarkPostProcessTimeStats(3) = obj.bmarkPostProcessTimeStats(3) + 1;
                %                fprintf('Total Time (%d scans/%d chans): %g -- Read: <%g, %g, %g> Plot: <%g, %g, %g> PreProc: <%g, %g, %g> PostProc: <%g, %g, %g> (Format: <last, mean, std>)\n', ...
                %                    scansToRead,size(newData,2),1000*t9,...
                %                    readTime,obj.bmarkReadTimeStats(1),obj.bmarkReadTimeStats(2),...
                %                    plotTime,obj.bmarkPlotTimeStats(1),obj.bmarkPlotTimeStats(2),...
                %                    procTimePre,obj.bmarkPreProcessTimeStats(1),obj.bmarkPreProcessTimeStats(2),...
                %                    procTimePost,obj.bmarkPostProcessTimeStats(1),obj.bmarkPostProcessTimeStats(2));
                
            catch ME %Handle Timer CB errors
                most.idioms.reportError(ME);
                ME.rethrow();
            end
            
            
            return;
            
            function bufStartScanNum = znstAugmentFullDataBuffer(scansToRead, newData)
                %dfprintf('augmentFullDataBuffer! scansToRead: %d\n', scansToRead);
                %dfprintf('newdata size: %d\n', size(newData,1));
                assert(ismember(size(obj.fullDataBuffer,1),[0 diff(obj.horizontalRangeScans)+1 obj.stimEventClassifyNumScans - 1]),'Expected fullDataBuffer to be empty or exactly equal to size of spike window');
                
                %         if obj.bufScanNumEnd == 0
                %         else
                %           obj.bufScanNumEnd = obj.bufScanNumEnd + scansToRead; %End index of augmented fullDataBuffer
                %         end
                
                obj.bufScanNumEnd = obj.maxReadableScanNum;
                bufStartScanNum = obj.bufScanNumEnd - scansToRead - size(obj.fullDataBuffer,1); %Start index of fullDataBuffer (including previously read samples carried over from last timer batch, the last post-window worth not yet processed)
                
                %dfprintf('before augment: fulldatabuffer size: %d\n', size(obj.fullDataBuffer,1));
                fprintf('scanstoRead: %d, size of full data buffer: %d, size of array: %d\n',scansToRead, obj.fullDataBuffer, size(newData));

                obj.fullDataBuffer = [obj.fullDataBuffer; newData];
                
                fprintf('scansToRead: %d, size of full data buffer: %d, size of array: %d\n',scansToRead, obj.fullDataBuffer, size(newData));
                
                %dfprintf('after augment: fulldatabuffer size: %d\n', size(obj.fullDataBuffer,1));
            end
            
            function znstStoreReducedData(timestampOffsets,bufStartScanNum)
                
                %timestampOffsets is a scalar in
                %stimulusTriggeredWaveformMode;and a (multi-channel) cell
                %array of (multi-spike) variable-length arrays in
                %spike-triggered waveform mode
                
                if obj.debug && ~isempty(timestampOffsets{1})
                    % ******************* DEBUG ********************
                    fprintf('bufStartScanNum: %d, timestampOffsets: %s, delta: %d\n',bufStartScanNum, sprintf('%d ',timestampOffsets{1}),timestampOffsets{1}-bufStartScanNum);
                    % ******************* DEBUG ********************
                end
                
                if stimulusTriggeredWaveformMode
                    assert(length(timestampOffsets) <= 1, 'Unexpectedly detected more than one stimulus in stimulus triggered waveform mode.');
                else
                    if isempty(timestampOffsets)
                        return;
                    end
                end
                
                try
                   
                    scanWindowRelative = obj.horizontalRangeScans(1):obj.horizontalRangeScans(2);
                    
                    for h=1:numel(obj.neuralChanAcqList)
                        i = obj.sglChanSubset(h)+1;
                        %TODO: Where possible, short-circuit storage for
                        %channels not being displayed, to reduce processing
                        %time
                        
                        if stimulusTriggeredWaveformMode
                            timestampOffsets_ = timestampOffsets;
                        else
                            timestampOffsets_ = timestampOffsets{i};
                        end
                        numNewTimestamps = length(timestampOffsets_);                                                
                        
                        %In waveform display modes - clear all previous channel data
                        if waveformDisplay
                            obj.reducedData{i}.scanNums = [];
                            obj.reducedData{i}.waveforms = cell(numNewTimestamps,1);
                        end
                        
                        if numNewTimestamps == 0
                            continue;
                        end
                        
                        %Store timestamp data
                        if stimulusTriggeredWaveformMode
                            obj.reducedData{i}.scanNums = [obj.reducedData{i}.scanNums timestampOffsets];
                        elseif waveformDisplay %spike-triggered waveform mode
                            obj.reducedData{i}.scanNums = timestampOffsets{i};
                        else %raster mode
                            obj.reducedData{i}.scanNums = [obj.reducedData{i}.scanNums timestampOffsets{i}];
                            obj.reducedData{i}.stimRelScanNums = [obj.reducedData{i}.stimRelScanNums zeros(1,numNewTimestamps)];
                            obj.reducedData{i}.stimNums = [obj.reducedData{i}.stimNums zeros(1,numNewTimestamps)];
                            obj.reducedData{i}.stimEventTypes = [obj.reducedData{i}.stimEventTypes repmat({''},1,numNewTimestamps)];
                        end
                        
                        %Store waveform data
                        if waveformDisplay
                            znstStoreReducedData_Waveforms();
                        end
                        
                    end
                    
                catch ME
                    fprintf(1,'Error: %s\n',ME.message);
                    ME.rethrow();
                end
                
                
                function znstStoreReducedData_Waveforms()
                    idxWindowMin = scanWindowRelative + timestampOffsets_(1) - bufStartScanNum;
                    idxWindowMax = scanWindowRelative + timestampOffsets_(numNewTimestamps) - bufStartScanNum;
                    
                    %size(obj.fullDataBuffer)
                    
                    if obj.debug && h==1 % Limit debug to channel 1 to keep things simple.
                        % ****************************** DEBUG *****************************
                        fprintf('\nscanNums: %d, numNewTimestamps: %d\n', obj.reducedData{i}.scanNums, numNewTimestamps);
                        fprintf('timestampOffsets: %s\n',sprintf('%d ',timestampOffsets_));
                        fprintf('scanWindowRelative: %s\n',sprintf('%d ',scanWindowRelative));
                        fprintf('idxWindowMin: %s\n',sprintf('%d ',idxWindowMin));
                        fprintf('idxWindowMax: %s\n',sprintf('%d ',idxWindowMax));
                        fprintf('\n************************************************************************\n');
                        % ****************************** DEBUG *****************************
                    end
                    
                    for j=1:numNewTimestamps
                        if stimulusTriggeredWaveformMode
                            scanWindow = scanWindowRelative + timestampOffsets_;
                        else
                            scanWindow = scanWindowRelative + timestampOffsets_(j);
                        end
                        
                        idxWindow = scanWindow - bufStartScanNum;
                        
                        if obj.debug && h==1
                            % ****************************** DEBUG *****************************
                            fprintf('numNewTimestamps: %d, j: %d, idxWindow: %s\n',numNewTimestamps,j,sprintf('%d ',idxWindow));
                            fprintf('j: %d, idxWindow: %s\n',j,sprintf('%d ',idxWindow));
                            fprintf('j: %d, idxWindow < 1: %s\n',j,sprintf('%d ',idxWindow < 1));
                            fprintf('j: %d, idxWindow >= 1: %s\n',j,sprintf('%d ',idxWindow >= 1));
                            % ****************************** DEBUG *****************************
                        end
                        
                        %Store waveform(s) for this channel, checking for edge cases
                        if find(idxWindow < 1) %'early' waveform in first timer chunk upon start of Spoke execution, before there is any history from prior chunks
                            tmpWindow = idxWindow + 1000 + abs(scanWindowRelative(1));
                            
                            offset = find(~idxWindow); % Find the offset in the waveform of the first non-negative value. This is the true "trigger level" crossing location. Don't lose this.
                            waveform = zeros(length(tmpWindow),1,'int16');
                            % Note: The following operations used to assume
                            % that the zero value in idxwindow/tmpwindow was
                            % the "trigger level" crossing. In actuality,
                            % this is relative to the relative scan window.
                            % So what we needed to do was renormalize the
                            % idxWindow to the relative scan window. Then we
                            waveform(tmpWindow < 1) = obj.fullDataBuffer(1,h); % Pad the waveform to the left of the zero crossing with whatever the first value in the waveform is.
                            waveform(tmpWindow >= 1) = obj.fullDataBuffer(tmpWindow >= 1,h); % Fill in the waveform to the right of the zero crossing with the actual waveform.
                            newwaveform = vertcat(zeros(offset,1),waveform); % Now, center the waveform on the original "trigger level" crossing index.
                            obj.reducedData{i}.waveforms{j} = newwaveform(1:length(tmpWindow));
                            
                            if obj.debug &&  h == 1
                                % ********************************** DEBUG ********************************
                                fprintf('yabba dabba doo\n');
                                fprintf('numNewTimestamps: %d, j: %d, tmpWindow: %s\n',numNewTimestamps, j,sprintf('%d ',tmpWindow));
                                fprintf('numNewTimestamps: %d, j: %d, waveform: %s\n',numNewTimestamps, j,sprintf('%d ',waveform));

                                if isempty(obj.diagrammm)
                                    obj.diagrammm = figure('name','temporary plot 2');
                                end
                                
                                figure(obj.diagrammm);
                                cla;
                                xlim manual;
                                xlim([0 600]);
                                ylim manual;
                                ylim([0, 17000]);
                                hold on;
                                
                                plot(waveform);
                                hold off;
                                % ********************************** DEBUG ********************************
                            end
                        elseif (idxWindowMax(end) > size(obj.fullDataBuffer,1)) %'late' waveform at end of a timer chunk, extending beyond available data
                            if stimulusTriggeredWaveformMode
                                if h==1 % EKK - Revisit this condition. Why does the first channel require different processing?
                                    obj.partialWaveformBuffer{end+1,i} = obj.fullDataBuffer(idxWindowMin(1):end,h);
                                    newWaveformWrapVal = idxWindowMax(end) - size(obj.fullDataBuffer,1); %TODO CHECK FOR +1 NECESSARY?
                                else
                                    obj.partialWaveformBuffer{end,i} = obj.fullDataBuffer(idxWindowMin(1):end,h);
                                end
                            else
                                %Warn the user if a 'late' waveform arises in spike-triggered waveform mode
                                fprintf('waveform #%d out of bounds, channel #%d, idxWindow min(1): %d, idxWindow min(end): %d, idxWindow max(1): %d, idxWindow max(end): %d\n', ...
                                    j, i, idxWindowMin(1), idxWindowMin(end), idxWindowMax(1), idxWindowMin(end) );
                                
                                %TODO:Determine if this can actually happen; if the partialWaveformBuffer/waveformWrap
                                %mechanism is needed/useful in spike-triggered waveform mode; and/or if the contract/augment buffer mechanism
                                %supercedes or is superceded by the partialWaveformBuffer/waveformWrap potential mechanism
                            end
                        else
                            obj.reducedData{i}.waveforms{j} = obj.fullDataBuffer(idxWindow,h);
                        end
                    end
                end
            end
            
            function znstDetectStimulus(bufStartScanNum)
                
                %Don't detect stimulus start within horizontal range (waveform or raster sweep) defined by last detected stimulus
                if isempty(obj.stimLastEventScanNumWindow)
                    reducedDataBufStartIdx = 1;
                elseif obj.stimLastEventScanNumWindow(2) >= bufStartScanNum
                    reducedDataBufStartIdx = obj.stimLastEventScanNumWindow(2) - bufStartScanNum + 1;
                else
                    reducedDataBufStartIdx = 1;
                end
                
                %Detect & record stimulus start and associated stimulus window
                %                stimIdx = find(diff(obj.fullDataBuffer(reducedDataBufStartIdx:end,obj.stimStartChannel + 1) > (obj.stimStartThreshold / obj.voltsPerBitNeural)) == 1, 1); %Should not have off-by-one error -- lowest possible value is fullDataBufferStartIdx+1 (if the second sample crosses threshold)
                stimChanFullDataIdx = find(obj.sglChanSubset==obj.stimStartChannel);
                % stimIdx = find(diff(obj.fullDataBuffer(reducedDataBufStartIdx:end,stimChanFullDataIdx) > (obj.stimStartThreshold / obj.voltsPerBitAux)) == 1, 1); %Should not have off-by-one error -- lowest possible value is fullDataBufferStartIdx+1 (if the second sample crosses threshold)
                %  stimIdx = find(diff(obj.fullDataBuffer(reducedDataBufStartIdx:end,stimChanFullDataIdx)) > (obj.stimStartThreshold / obj.voltsPerBitAux)) == 1; %Should not have off-by-one error -- lowest possible value is fullDataBufferStartIdx+1 (if the second sample crosses threshold)
                triggerThreshVal = obj.stimStartThreshold / obj.voltsPerBitAux;
                stimIdx = find((diff(obj.fullDataBuffer(reducedDataBufStartIdx:end,stimChanFullDataIdx)) > triggerThreshVal) == 1,1); %Fixed - Ed
                %dfprintf('stimChanFullDataIdx: %d min data: %g max data: %g\n', stimChanFullDataIdx, min(obj.fullDataBuffer(reducedDataBufStartIdx:end,stimChanFullDataIdx)), max(obj.fullDataBuffer(reducedDataBufStartIdx:end,stimChanFullDataIdx)));
                if ~isempty(stimIdx)
                    newStimScanNum = bufStartScanNum + stimIdx - 1;
                    obj.stimScanNums(end+1) = newStimScanNum;
                    
                    if rasterDisplayMode
                        obj.stimLastEventScanNumWindow = newStimScanNum + round(obj.horizontalRangeRaster/sampPeriod);
                    else
                        obj.stimLastEventScanNumWindow = newStimScanNum + round(obj.horizontalRange/sampPeriod);
                    end
                    
                    obj.stimWindowStartScanNums(end+1) = obj.stimLastEventScanNumWindow(1);
                    obj.stimWindowEndScanNums(end+1) = obj.stimLastEventScanNumWindow(2);
                    obj.stimEventTypeNames{end+1} = '';
                    
                    obj.stimTotalCount = obj.stimTotalCount + 1;
                    dfprintf('Detected stim! stimTotalCount: %d stimScanNum: %d stimWindowStartScanNum: %d stimWindowEndScanNum: %d bufStartScanNum: %d\n',...
                        obj.stimTotalCount,obj.stimScanNums(end),obj.stimWindowStartScanNums(end),obj.stimWindowEndScanNums(end),bufStartScanNum);
                end
            end
            
            
            function znstClassifyStimulus(bufStartScanNum)
                
                %Loop through detected stimuli that have not been classified, from the end - classifying them if possible.
                stimNum = length(obj.stimScanNums);
                while isempty(obj.stimEventTypeNames{stimNum})
                    
                    if isempty(obj.stimEventTypes)
                        obj.stimEventTypeNames{stimNum} = 'allstim';
                        obj.stimEventCount_.allstim = obj.stimEventCount_.allstim + 1;
                    elseif (obj.stimScanNums(stimNum) + obj.stimEventClassifyNumScans) <= (bufStartScanNum + length(obj.fullDataBuffer)) %Have sufficient scans to classify the event
                        
                        startIdx = obj.stimScanNums(stimNum) - bufStartScanNum + 1;
                        endIdx = obj.stimScanNums(stimNum) - bufStartScanNum + obj.stimEventClassifyNumScans ;
                        
                        %fprintf('bufStartScanNum: %d sizereducedDataBuf: %s startIdx: %d endIdx: %d \n',bufStartScanNum,mat2str(size(obj.fullDataBuffer)),startIdx,endIdx);
                        
                        eventType = feval(obj.stimEventClassifyFcn,obj.fullDataBuffer(startIdx:endIdx,obj.stimEventClassifyChannel+1));
                        if isempty(eventType)
                            fprintf(2,'WARNING: Failed to classify stimulus event type!\n');
                            obj.stimEventTypeNames{stimNum} = 'unknown';
                        else
                            obj.stimEventTypeNames{stimNum} = eventType;
                            obj.stimEventCount_.(eventType) = obj.stimEventCount_.(eventType) + 1;
                        end
                    end
                    stimNum = stimNum - 1;
                    if stimNum == 0
                        break;
                    end
                end
                
            end
            
            function chanNewSpikes = znstAssociateSpikesToStimuli()
                %For Raster mode, associate spikes to stimuli (and stimulus
                %event-types, if applicable). Update reduced data structure
                %with these associations. Return (just) the newly
                %associated spikes, to allow for a fast plot update.
                                
                for i=1:length(obj.stimEventTypes_)
                    chanNewSpikes.(obj.stimEventTypes_{i}) = zeros(numNeuralChans,1);
                end
                
                %numNewSpikes = 0;
                
                taggedSpikeIdxStructInit = cell2struct(repmat({[]},length(obj.stimEventTypes_),1),obj.stimEventTypes_);
                
                %                 for c=1:numNeuralChans
                %for b=1:numel(obj.sglChanSubset) % EDKANG
                for b=1:numel(obj.neuralChanAcqList)        % EDKANG
                    c=obj.sglChanSubset(b) + 1; % EDKANG
                    
                    taggedNewSpike = false;
                    spikesToClear = [];
                    taggedSpikeIdxsStruct = taggedSpikeIdxStructInit;
                    %                    fprintf('b: %d    ',b);
                    %Loop through spikes from most recent backwards, tagging event-type if possible
                    
                    %                    fprintf('length of spikedata: %d\n',length(obj.spikeData{c}.scanNums));
                    for spikeIdx = length(obj.reducedData{c}.scanNums):-1:1
                        %                        fprintf('hello!!!! %d\n',spikeIdx);
                        tmp1 =tic;
                        
                        %Reached previously-tagged spikes -- stop loop
                        if obj.reducedData{c}.stimNums(spikeIdx) > 0
                            break;
                        end
                        
                        %VI051112: Unlike previously planned - reducedData won't be deleted on plotting.
                        %assert(isempty(obj.reducedData{c}.stimEventTypes{spikeIdx}),'reducedData was stored longer than expected for an already stimulus event-tagged spike');
                        
                        %Find associated stim spike
                        spikeScanNum = obj.reducedData{c}.scanNums(spikeIdx);
                        
                        % fprintf('spikeScanNum: %d\n',spikeScanNum);
                        
                        stimIdx = find(spikeScanNum >= obj.stimWindowStartScanNums,1,'last');
                        
                        if ~isempty(stimIdx) && spikeScanNum <= obj.stimWindowEndScanNums(stimIdx)
                            %Found associated stimulus: store stim info & event tag
                            if isempty(obj.stimEventTypes)
                                eventTag = 'allstim'; %All stim-associated spikes get plotted
                            else
                                eventTag = obj.stimEventTypeNames{stimIdx}; %Note - this might still be empty. That's OK...will try again to tag it the next time around.
                            end
                            
                            if ~isempty(eventTag) && isempty(obj.reducedData{c}.stimEventTypes{spikeIdx})
                                obj.reducedData{c}.stimEventTypes{spikeIdx} = eventTag;
                                taggedSpikeIdxsStruct.(eventTag)(end+1) = spikeIdx; %Add spikeIdx to list of spikes associated with this eventTag
                                
                                obj.reducedData{c}.stimNums(spikeIdx) = stimIdx;
                                obj.reducedData{c}.stimRelScanNums(spikeIdx) = spikeScanNum - obj.stimScanNums(stimIdx);
                                
                                chanNewSpikes.(eventTag)(c) = chanNewSpikes.(eventTag)(c) + 1;
                                taggedNewSpike = true;
                            end
                            
                        elseif spikeScanNum < (obj.maxReadableScanNum + round(obj.horizontalRangeRaster(1)/sampPeriod))
                            %No hope of ever finding associated stimulus: mark for deletion
                            spikesToClear(end+1) = spikeIdx;
                        else
                            %fprintf('spikeScanNum: %d maxReadableScanNum: %d preStimTime: %d mostRecentHopefulScan: %d\n',spikeScanNum,obj.maxReadableScanNum, round(obj.horizontalRangeRaster(1)/sampPeriod),(obj.maxReadableScanNum - round(obj.horizontalRangeRaster(1)/sampPeriod)));
                            %Do nothing -- spike still has hope of finding associated stimulus
                        end
                    end
                    
                    tmp1 = tic;
                    %Maintain indices of stored spikes associated with each event, for per-event lookup %TODO: Determine if this speedup is actually apparent/important
                    if taggedNewSpike
                        for i=1:length(obj.stimEventTypes_)
                            obj.reducedData{c}.stimEventTypeStruct.(obj.stimEventTypes_{i}) =  [obj.reducedData{c}.stimEventTypeStruct.(obj.stimEventTypes_{i}) (taggedSpikeIdxsStruct.(obj.stimEventTypes_{i}) - length(spikesToClear))];
                        end
                    end
                    
                    
                    %Clear 'orphan' spikes with no hope of finding associated stimulus
                    tmp1 = tic;
                    if ~isempty(spikesToClear)
                        obj.reducedData{c}.scanNums(spikesToClear) = [];
                        obj.reducedData{c}.stimNums(spikesToClear) = [];
                        obj.reducedData{c}.stimRelScanNums(spikesToClear) = [];
                        obj.reducedData{c}.stimEventTypes(spikesToClear) = [];
                    end
                end
            end
            
            function znstUpdateBaselineStats(newSpikeScanNums,bufStartScanNum)
                %Extract fullDataBuffer data, excluding last spike-window post
                %time (not yet processed for spikes), and excluding just-detected
                %spike windows
                
                %newRmsData = cell(numDispChans,1);
                %batchLength = zeros(numDispChans,1);
                rmsDataIdxs = {1:(size(obj.fullDataBuffer,1)-obj.horizontalRangeScans(2))};
                rmsDataIdxs = repmat(rmsDataIdxs,numNeuralChans,1);
                
                firstPassMode = isempty(newSpikeScanNums); %Handle first pass at RMS detection, when there are no detected spikes yet
                
                if ~firstPassMode
                    %for i=1:numNeuralChans
                    %for i=1:numel(obj.sglChanSubset)
                    for i=1:numel(obj.neuralChanAcqList)
                        %newRmsData{i} = obj.fullDataBuffer(1:end-obj.horizontalRangeScans(2),i);
                        %batchLength(i) = length(newRmsData{i});
                        
                        %rmsDataIdxs{i} = 1:(size(obj.fullDataBuffer,1)-obj.horizontalRangeScans(2));
                        
                        spikeScanIdxs = newSpikeScanNums{i} - bufStartScanNum + 1;  %Convert scan numbers to spike-data-buffer indices
                        
                        %Exclude just-detected spike windows
                        badIdxs = [];
                        for j=1:length(spikeScanIdxs)
                            %newRmsData{i}(spikeScanIdxs:min(end,(spikeScanIdxs+obj.horizontalRangeScans(2)))) = [];
                            badIdxs = [badIdxs spikeScanIdxs(j):min(rmsDataIdxs{i}(end),spikeScanIdxs(j)+obj.horizontalRangeScans(2))];
                        end
                        
                        % remove all negative values
                        index = find(badIdxs>0);
                        badIdxs = badIdxs(index);
                        
                        rmsDataIdxs{i}(badIdxs) = [];
                    end
                end
                
                % Update mean & RMS computation for each pad channel
                warnNoData = false;
                %for i=1:numNeuralChans
                for i=1:numel(obj.neuralChanAcqList)
                    if isempty(rmsDataIdxs{i})
                        if ~warnNoData
                            warning('For at least one channel (%d), no data was available for RMS/mean calculations',i);
                            warnNoData = true;
                        end
                        continue;
                    end
                    
                    %fprintf('Computing mean & RMS for channel %d. Num idxs to average: %d\tSDB size: %s\n',i,length(rmsDataIdxs{i}),mat2str(size(obj.fullDataBuffer)));
                    dataLen = length(rmsDataIdxs{i});
                    if obj.filterWindow(1) > 0  %No need for per-channel mean subtraction if highpass-filtering
                        obj.baselineMean(i) = 0;
                        obj.baselineRMS(i) = sqrt(sum(double(obj.fullDataBuffer(rmsDataIdxs{i},i)).^2)/dataLen);
                    else %Use per-channel mean subtraction
                        obj.baselineMean(i) = sum(double(obj.fullDataBuffer(rmsDataIdxs{i},i)))/dataLen;
                        obj.baselineRMS(i) = sqrt(sum((double(obj.fullDataBuffer(rmsDataIdxs{i},i)) - obj.baselineMean(i)).^2)/dataLen);
                    end
                end
                
                %Redraw threshold lines, if it can change (only in case of 'mismatched' threshold type and display units)
                if strcmpi(obj.thresholdType,'rmsMultiple') && strcmpi(obj.waveformAmpUnits,'volts')
                    obj.zprvDrawThresholdLines();
                end
                
                %Update last scan at which RMS has been computed
                if ~firstPassMode
                    obj.baselineRMSLastScan = obj.bufScanNumEnd;
                end
                
            end
            
        end
        
        function newSpikeScanNums = zprvDetectNewSpikes(obj,bufStartScanNum)
            
            if strcmpi(obj.thresholdType,'rmsMultiple')
                if obj.baselineRMSLastScan == 0 %No RMS data has been computed yet
                    %newSpikeScanNums = cell(numDispChans,1);
                    threshVal = obj.baselineRMS * obj.INIT_RMS_THRESHOLD; %Use pre-set RMS multiplier for first buffer
                else
                    threshVal = obj.thresholdVal * obj.baselineRMS;
                end
                
                
                if obj.filterWindow(1) > 0
                    threshMean = 0; %Highpass filtering supercedes mean subtraction
                else
                    threshMean = obj.baselineMean;
                end
                
                [newSpikeScanNums, obj.maxNumWaveformsApplied] = ...
                    zlclDetectSpikes(obj.reducedData, ...
                                     obj.fullDataBuffer, ...
                                     bufStartScanNum, ...
                                     round(obj.spikeRefractoryPeriod * obj.sglParamCache.niSampRate), ...
                                     threshVal, ...
                                     obj.thresholdAbsolute, ...
                                     threshMean, ...
                                     obj.refreshPeriodMaxNumWaveforms, ...
                                     obj.sglChanSubset, ...
                                     obj.neuralChanAcqList, ...
                                     obj.horizontalRangeScans, ...
                                     obj.debug, ...
                                     obj.diagramm); %Detect spikes from beginning in all but the spike-window-post time, imposing a 'refractory' period of the spike-window-post time after each detected spike
                
                %             if maxNumWaveformsApplied && ~obj.maxNumWaveformsApplied
                %               fprintf(2,'WARNING: Exceeded maximum number of spikes (%d) on one or more channels; subsequent spikes were ignored.\n', obj.refreshPeriodMaxNumWaveforms);
                %             end
                %
                %             obj.maxNumWaveformsApplied = maxNumWaveformsApplied;
                
            else
                threshVal = obj.thresholdVal / obj.voltsPerBitNeural; %Convert to AD units
                threshMean = 0; %Don't do mean subtraction
                newSpikeScanNums = zlclDetectSpikes(obj.reducedData, ...
                                                    obj.fullDataBuffer, ...
                                                    bufStartScanNum, ...
                                                    round(obj.spikeRefractoryPeriod * obj.sglParamCache.niSampRate), ...
                                                    threshVal, ...
                                                    obj.thresholdAbsolute, ...
                                                    0, ...
                                                    obj.refreshPeriodMaxNumWaveforms, ...
                                                    obj.sglChanSubset, ...
                                                    obj.neuralChanAcqList, ...
                                                    obj.horizontalRangeScans, ...
                                                    obj.debug, ...
                                                    obj.diagramm); %Detect spikes from beginning in all but the spike-window-post time, imposing a 'refractory' period of the spike-window-post time after each detected spike
            end
            
        end
        
        function zprvUpdateRasterPlot(obj,chanNewSpikes)
            sampPeriod = 1 / obj.sglParamCache.niSampRate;
            colorOrder = get(0,'DefaultAxesColorOrder');
            
            plotAllSpikes = (nargin < 2);
            if plotAllSpikes
                obj.zprvResetStimNumsPlotted();
            end
            
            %Update vert axis limits if needed (i.e. if upper limit is set to Inf)
            ylim = obj.zprvverticalRangeRaster2YLim(obj.verticalRangeRaster);
            if ~isequal(get(obj.hRasters(1),'YLim'),obj.verticalRangeRaster)
                set(obj.hRasters,'YLim',ylim);
            end
            
            if isempty(obj.stimEventTypes)
                eventTypes = {'allstim'};
            else
                eventTypes = obj.stimEventTypesDisplayed;
            end
            
            if isscalar(eventTypes)
                eventType = eventTypes{1};
            end
            
            %Determine eventColorIdx up front, in cases where only one event type is being displayed
            if isempty(obj.stimEventTypes)
                eventColorIdx = 1;
            elseif isscalar(eventTypes)
                [tf,eventColorIdx] = ismember(eventType,obj.stimEventTypes);
                assert(tf);
            end
            
            %Add new line object for stim-associated spikes
            
            plotCount = 0;
            for c=obj.tabChanNumbers
                plotCount = plotCount + 1;
                % EKK - change channel number of tab to mapped channel
                % number.
                actualChannelNumber = obj.neuralChanDispOrder(c);

                
                
                %if isempty(obj.reducedData{j}) || chanNewSpikes(j) == 0 || isempty(obj.reducedData{j}.stimEventTypeStruct)
                if isempty(obj.reducedData{c}.scanNums) || (~plotAllSpikes && all(structfun(@(x)x(c)==0,chanNewSpikes)))
                    continue;
                end
                
                %Ensure only newly stim-associated spikes are plotted (unless refreshing whole plot)
                if isscalar(eventTypes)
                    plotSpikeIdxs = obj.reducedData{c}.stimEventTypeStruct.(eventType);
                    if ~plotAllSpikes
                        plotSpikeIdxs(1:(length(plotSpikeIdxs) - chanNewSpikes.(eventType)(c))) = [];
                    end
                else
                    %Find last-tagged spike for this channel
                    
                    for lastTaggedSpikeIdx=length(obj.reducedData{c}.stimEventTypes):-1:0
                        if lastTaggedSpikeIdx > 0 && ~isempty(obj.reducedData{c}.stimEventTypes{lastTaggedSpikeIdx})
                            break;
                        end
                    end
                    if lastTaggedSpikeIdx==0
                        break;
                    end
                    
                    if plotAllSpikes
                        plotSpikeIdxs = 1:lastTaggedSpikeIdx;
                    else
                        allEventTypes = obj.stimEventTypes;
                        numNewSpikes = 0;
                        for i=1:length(allEventTypes)
                            numNewSpikes = numNewSpikes + chanNewSpikes.(allEventTypes{i})(c);
                        end
                        
                        plotSpikeIdxs = (lastTaggedSpikeIdx-numNewSpikes+1):lastTaggedSpikeIdx;
                    end
                end
                
                
                %Number stims by their order within the event type(s) selected
                if isscalar(eventTypes)
                    stimNumsPlotted = obj.stimNumsPlotted(c).(eventType);
                else
                    stimNumsPlotted = obj.stimNumsPlotted{c};
                end
                
                [uniqueStimNums,~,orderedStims] = unique(obj.reducedData{c}.stimNums(plotSpikeIdxs),'sorted'); %orderedStims has same length as plotSpikeIdxs
                
                %TODO: Remove this when enough testing confirms this never happens
                haveUntaggedStims = ~isempty(uniqueStimNums) && any(uniqueStimNums == 0);
                assert(~haveUntaggedStims,'The plotSpikeIdxs identified contain one or more untagged spikes');
                
                %Determine ordering of stims in relation to previously plotted lines (determine stims that have been plotted, partially or fully before, for determining base count)
                repeatStims = [];
                idx = length(stimNumsPlotted);
                while idx > 0
                    if ismember(stimNumsPlotted(idx),uniqueStimNums)
                        repeatStims(end+1) = stimNumsPlotted(idx);
                        idx = idx - 1;
                    else
                        break;
                    end
                end
                
                orderedStimBaseNum = length(stimNumsPlotted) - length(repeatStims);
                orderedStims = orderedStims + orderedStimBaseNum;
                
                if isscalar(eventTypes)
                    %Add points to line of new stim-associated spikes --
                    %note 1) line can have points spanning more than one
                    %stim and 2) may not (typically won't) have all the
                    %spikes of a given stim
                    obj.hRasterLines{1}(plotCount).addpoints(obj.reducedData{c}.stimRelScanNums(plotSpikeIdxs) * sampPeriod, orderedStims);
                    
                else
                    startIdxs = [1; find(diff(orderedStims))+1];
                    endIdxs = [startIdxs(2:end) - 1; length(plotSpikeIdxs)];
                    
                    assert(length(startIdxs) == length(uniqueStimNums));
                    
                    if plotAllSpikes
                        %Determine event type for each uniqueStimNumber, encoded as an index into the eventTypes array
                        [~,eventClassOrder] = ismember(obj.stimEventTypeNames(uniqueStimNums),eventTypes);
                        
                        %Plot stimulus spike records for each event type, in stimulus order within each event type grouping
                        stimsPlotted = length(stimNumsPlotted);
                        for i=1:length(eventTypes)
                            stimNums = find(eventClassOrder == i);
                            
                            for j=1:length(stimNums)
                                stimNum = stimNums(j);
                                startIdx = startIdxs(stimNum);
                                endIdx = endIdxs(stimNum);
                                
                                %Plot line of new stim-associated spikes for a particular stimulus. Note line may not (typically won't) have all the spikes for the given stim
                                line('Parent',obj.hRasters(plotCount),'XData',obj.reducedData{c}.stimRelScanNums(plotSpikeIdxs(startIdx:endIdx)) * sampPeriod,'YData', (stimsPlotted+1) * ones(endIdx-startIdx+1,1),'Marker','d','MarkerSize',2,'Color',colorOrder(i,:),'LineStyle','none'); %'EraseMode','none');
                                stimsPlotted = stimsPlotted + 1;
                            end
                            
                        end
                        
                    else
                        %Plot stimulus spike records in temporal order, color-coded by event type
                        for i=1:length(uniqueStimNums)
                            
                            [~,eventColorIdx] = ismember(obj.stimEventTypeNames(uniqueStimNums(i)),obj.stimEventTypes);
                            %assert(tf);
                            
                            startIdx = startIdxs(i);
                            endIdx = endIdxs(i);
                            
                            %Plot line of new stim-associated spikes for a particular stimulus. Note line may not (typically won't) have all the spikes for the given stim
                            line('Parent',obj.hRasters(plotCount),'XData',obj.reducedData{c}.stimRelScanNums(plotSpikeIdxs(startIdx:endIdx)) * sampPeriod,'YData', orderedStims(startIdx:endIdx),'Marker','d','MarkerSize',2,'Color',colorOrder(eventColorIdx,:),'LineStyle','none'); %,'EraseMode','none');
                            
                        end
                    end
                end
                
                %Update stimNumsPlotted
                if isscalar(eventTypes)
                    obj.stimNumsPlotted(c).(eventType) = [stimNumsPlotted setdiff(uniqueStimNums,repeatStims)];
                else
                    obj.stimNumsPlotted{c} = [stimNumsPlotted setdiff(uniqueStimNums,repeatStims)];
                end
                
            end
            
        end
        
        function zprvUpdateWaveformPlot(obj)
            
            totalNewWaveforms = 0;
            
            %totalClearedSpikes = 0;
            for i=obj.tabChanNumbers
                if isempty(obj.reducedData{i})
                    continue;
                end
                plotIdx = mod(i-1,obj.PLOTS_PER_TAB) + 1;
                
                
                % EKK - change channel number of tab to mapped channel
                % number.
                actualChannelNumber = obj.neuralChanDispOrder(i);
                
                numNewWaveforms = length(obj.reducedData{i}.scanNums);
                totalNewWaveforms = totalNewWaveforms + numNewWaveforms;
                
                %Plot new spike lines
                horizontalRangeScansLength = diff(obj.horizontalRangeScans)+1;
                xData = linspace(obj.horizontalRange(1),obj.horizontalRange(2),horizontalRangeScansLength)';
                
                newWaveformCounts = obj.lastPlottedWaveformCount(i) + (1:numNewWaveforms);
                lineIdxs = mod(newWaveformCounts,obj.waveformsPerPlot) + 1; %The line object indices to use for these newly detected spikes
                
                % Clear spikes (if necessary)
                switch obj.waveformsPerPlotClearMode
                    case 'all'
                        if obj.lastPlottedWaveformCountSinceClear(i) + numNewWaveforms > obj.waveformsPerPlot
                            obj.hWaveforms(plotIdx).clearpoints();
                            obj.lastPlottedWaveformCountSinceClear(i) = 0;
                        end
                    case 'oldest'
                        
                    otherwise
                        disp('invalid plot clear mode');
                end
                
                % Draw new spikes
                if ~isempty(obj.waveformWrap)
                    j=0; % used only to pass assert in nested function below. J=0 in this case because we are referring to a previous waveform (as a product of a previously detected stimulus or spike)
                    znstPlotWaveform(obj.partialWaveformBuffer{1,i}); %plot oldest partial waveform stored
                else
                    for j=1:numNewWaveforms
                        znstPlotWaveform(obj.reducedData{i}.waveforms{j});
                    end
                end
            end
            
            %Clear oldest partial waveforms stored, now that they're plotted
            if ~isempty(obj.waveformWrap)
                obj.waveformWrap(1) = [];
                obj.partialWaveformBuffer(1,:) = [];
            end
            
            function znstPlotWaveform(waveform)
                if isempty(waveform)
                    return;
                end
                assert(length(waveform) == length(xData),'Waveform data for chan %d (%d), spike %d not of expected length (%d)\n',i,length(waveform),j,length(xData));             
                
                %Scale waveform from A/D units to target units, applying mean subtraction if thresholdType='rmsMultiple'
                switch obj.waveformAmpUnits
                    case 'volts'
                        if strcmpi(obj.thresholdType,'volts') %no mean subtraction...just show as is
                            waveform = double(waveform) * obj.voltsPerBitNeural;
                        else  %RMS-multiple threshold --> do mean subtraction
                            waveform = (double(waveform) - obj.baselineMean(i)) * obj.voltsPerBitNeural;
                        end
                    case 'rmsMultiple'
                        if obj.filterWindow(1) > 0
                            waveform = double(waveform) / obj.baselineRMS(i);
                        else %Use mean subtraction
                            waveform = (double(waveform) - obj.baselineMean(i)) / obj.baselineRMS(i);
                        end
                end
                %Update line object with waveform for current spike
                obj.hWaveforms(plotIdx).addpoints(vertcat(xData, NaN),vertcat(waveform, NaN)); % old
                obj.lastPlottedWaveformCount(i) = obj.lastPlottedWaveformCount(i) + 1;
                obj.lastPlottedWaveformCountSinceClear(i) = obj.lastPlottedWaveformCountSinceClear(i) + 1;
            end
            
        end
        
        function zprvSetAxesProps(obj,hAx)
            %Axes properties for spoke waveform grid axes
            set(hAx,'XTick',0,'YTick',0,'XGrid','on','YGrid','on','XTickLabel','','YTickLabel','');
        end
        
        
        function zprvShowDisplayFig(obj)
            currPosn = obj.gridFigPosition;
            displayType = obj.displayMode;
            set(obj.hFigs.(displayType),'Visible','on','Position',currPosn);
            switch displayType
                case 'waveform'
                    set([obj.hFigs.raster obj.hFigs.psth],'Visible','off');
                case 'raster'
                    set(obj.hFigs.waveform,'Visible','off');
            end
        end
        
        
        %     function zprvResetThreshold(obj)
        %
        %       %TODO(?): A smarter adjustment based on the last-cached RMS values, somehow handlign the variety across channels
        %         switch obj.waveformAmpUnits
        %           case 'volts'
        %               obj.thresholdVal = .1 * obj.hSpoke.x_fs;
        %           case 'rmsMultiple'
        %             obj.thresholdVal= 5;
        %         end
        %
        %         %Updates threshold value, and threshold lines in process
        %         obj.thresholdVal = newThreshVal;
        %     end
        %
        %     function zprvResetUnits(obj)
        %
        %
        %     end
        
        function zprvDrawThresholdLines(obj)
            
            numNeuralChans = numel(obj.neuralChansAvailable);
            
            %Clear existing threshold lines
            handlesToClear = [obj.hThresholdLines{1}(isgraphics(obj.hThresholdLines{1})); obj.hThresholdLines{2}(isgraphics(obj.hThresholdLines{2}))];
            %set(handlesToClear,'EraseMode','normal');
            delete(handlesToClear);
            
            %Compute all-channel threshold; determine lack of threshold val -- as applicable
            perChanThreshold = ~strcmpi(obj.thresholdType,obj.waveformAmpUnits);
            if perChanThreshold %RMS threshold with voltage units -- this is only mismatch type presently allowed
                if isempty(obj.baselineRMS)
                    obj.hThresholdLines = repmat({ones(numNeuralChans,1) * -1},2,1);
                    return; %nothing to draw
                end
            else %matched units/threshold-type
                threshold = obj.thresholdVal;
            end
            
            %Draw threshold for each channel, computing for each channel if needed
            xData = obj.horizontalRange;
            
            for i=obj.tabChanNumbers
                plotIdx = mod(i-1,obj.PLOTS_PER_TAB) + 1;
                
                if perChanThreshold
                    threshold = obj.thresholdVal * obj.baselineRMS(i) * obj.voltsPerBitNeural;
                end
                
                if ~isempty(threshold)
                    dataArgs = {'XData',xData,'YData',[threshold threshold]};
                    
                    if obj.thresholdAbsolute
                        dataArgs2 = {'XData',xData,'YData',[-threshold -threshold]};
                    end
                    
                    if numel(obj.hThresholdLines{1}) < plotIdx || ~isgraphics(obj.hThresholdLines{1}(plotIdx))
                        obj.hThresholdLines{1}(plotIdx) = line('Parent',obj.hPlots(plotIdx),'Color','r',dataArgs{:}); %'EraseMode','none',
                        
                        if obj.thresholdAbsolute
                            obj.hThresholdLines{2}(plotIdx) = line('Parent',obj.hPlots(plotIdx),'Color','r',dataArgs2{:}); %'EraseMode','none',
                        end
                    else
                        set(obj.hThresholdLines{1}(plotIdx),dataArgs{:});
                        if obj.thresholdAbsolute
                            set(obj.hThresholdLines{2}(plotIdx),dataArgs2{:});
                        end
                    end
                end
                
            end
            
        end
        
        function zprvResetAcquisition(obj,fileRollover)
            %Method to handle acquisition 'reset' cases, either hard (SpikeGL is stopped & restarted) or soft (file rollover)
            
            if nargin < 2
                fileRollover = false;
            end
            
            numNeuralChans = numel(obj.neuralChansAvailable);
            
            obj.fullDataBuffer = zeros(0,numel(obj.neuralChanDispList) + numel(obj.auxChanProcList));
            
            if ~fileRollover %&& strcmpi(obj.thresholdType,'rmsMultiple')
                obj.baselineRMS = zeros(numNeuralChans,1);
                obj.baselineMean = zeros(numNeuralChans,1);
                obj.baselineRMSLastScan = 0;
            end
            
            if ~isempty(obj.filterCoefficients)
                obj.filterCondition = zeros(max(cellfun(@length,obj.filterCoefficients)) - 1, 1);
            end
            
            obj.bufScanNumEnd = 0;
            obj.stimLastEventScanNumWindow = [];
            
            if ~fileRollover
            end
            
        end
        
        function zprvResetReducedData(obj)
            %Method to clear cached spike data; can be either on acquisition 'reset' or in some cases mid-acquisition
            
            % TODO: Does this work with channel subsets?
            
            numNeuralChans = numel(obj.neuralChansAvailable);
            
            obj.lastPlottedWaveformCount = zeros(numNeuralChans,1);
            obj.lastPlottedWaveformCountSinceClear = zeros(numNeuralChans,1);
            
            obj.reducedData = cell(numNeuralChans,1);
            for i=1:1:numNeuralChans
                %for i=1:1:numel(obj.sglChanSubset)
                if strcmpi(obj.displayMode,'waveform')
                    obj.reducedData{i} = struct('scanNums',[],'waveforms',{{}});
                else
                    obj.reducedData{i} = struct('scanNums',[],'stimRelScanNums',[],'stimNums',[],'stimEventTypes',{{}},'stimEventTypeStruct',struct());
                    %             obj.reducedData{i} = struct();
                    %
                    %             eventTypes = obj.stimEventTypes;
                    %             if isempty(eventTypes)
                    %               obj.reducedData{i}.plot = struct('scanNums',[],'stimRefScanNum',[],'stimNum',[]);
                    %             else
                    %               for j=1:length(eventTypes)
                    %                 obj.reducedData{i}.(eventTypes{j}) = struct('scanNums',[],'stimRefScanNum',[],'stimNum',[]);
                    %               end
                    %             end
                end
            end
        end
        
        function zprvResetStimNumsPlotted(obj)
            obj.stimNumsPlotted = struct(); %Clears existing struct data
            
            nca = numel(obj.neuralChansAvailable);
            if isempty(obj.stimEventTypes)
                obj.stimNumsPlotted(nca).allstim = [];
            elseif isscalar(obj.stimEventTypesDisplayed)
                for i=1:length(obj.stimEventTypes)
                    obj.stimNumsPlotted(nca).(obj.stimEventTypes{i}) = [];
                end
            else
                obj.stimNumsPlotted = cell(nca,1); %Cell array of empty arrays
            end
            
        end
        
        function zprvClearPlots(obj,displaysToClear,reuseThreshold)
            % displaysToClear: String or string cell array containing one or more of {'waveform' 'raster' 'psth'}. Used to signal this operation is specifically to clear plots for one of the display modes in particular.
            % reuseThreshold: <Default=true> If true, redraw the threshold lines based on the value that was prevailing before the plots are cleared
            
            if nargin < 3
                reuseThreshold = true;
            end
            
            %Don't bother clearing a plot that's invisible (or should be)
            if ischar(displaysToClear)
                displaysToClear = {displaysToClear};
            end
            
            shouldBeCleared = setdiff(displaysToClear,{'psth' obj.displayMode});
            if ~isempty(shouldBeCleared)
                cellfun(@(x)set(obj.hFigs.(x),'Visible','off'),shouldBeCleared); %Should already be off
            end
            
            if isempty(setdiff(displaysToClear,shouldBeCleared))
                return;
            end
            
            if nargin < 3
                refreshThreshold = true;
            end
            
            redrawThresholdLines = false;
            
            for i=1:length(displaysToClear)
                %for j=1:min(obj.PLOTS_PER_TAB,numel(obj.hThresholdLines{1}))
                for j=1:min(obj.PLOTS_PER_TAB)
                    
                    displayToClear = displaysToClear{i};
                    
                    switch displayToClear
                        case 'waveform'
                            
                            %Clear out graphics
                            reuseThreshold = isgraphics(obj.hThresholdLines{1}(j)) && reuseThreshold;
                            if reuseThreshold
                                threshold = unique(get(obj.hThresholdLines{1}(j),'YData'));
                                
                                if isgraphics(obj.hThresholdLines{2}(j))
                                    threshold = [threshold -threshold];
                                end
                            end
                            
                            cla(obj.hPlots(j));
                            obj.zprvSetAxesProps(obj.hPlots(j));
                            
                            if reuseThreshold
                                for k=1:length(threshold)
                                    obj.hThresholdLines{k}(j) = line('Parent',obj.hPlots(j),'Color','r','XData',obj.horizontalRange,'YData',[threshold(k) threshold(k)]); %'EraseMode','none',
                                end
                            else %Compute/draw threshold lines from scratch
                                redrawThresholdLines = true;
                            end
                            
                        case 'raster'
                            for k = 1:length(obj.hRasterLines)
                                obj.hRasterLines{k}(j).clearpoints();
                            end
                            obj.zprvSetAxesProps(obj.hRasters(j));
                            
                        case 'psth'
                            cla(obj.hPSTHs(j));
                            obj.zprvSetAxesProps(obj.hPSTHs(j));
                    end
                    
                    
                    %preallocate animated lines for spike waveforms
                    %obj.hWaveforms(j) = animatedline('Parent',obj.hPlots(j),'MaximumNumPoints',Inf,'Marker','.','MarkerSize',3,'LineStyle','none');
                    obj.hWaveforms(j) = animatedline('Parent',obj.hPlots(j),'MaximumNumPoints',obj.maxPointsPerAnimatedLine,'Marker','.','MarkerSize',3,'LineStyle','-');
                end
            end
            
            if redrawThresholdLines
                obj.zprvDrawThresholdLines();
            end
            
        end
        
        function ylim = zprvverticalRangeRaster2YLim(obj,val)
            ylim = val + [-1 0]; %Ensure the first element is seen
            if isinf(ylim(2))
                numStimsToDisplay = max(obj.stimEventCount - ylim(1) + 1,1);
                
                numIncrementFactor = numStimsToDisplay/obj.verticalRangeRasterInfIncrement;
                numIncrementsNeeded = ceil(numIncrementFactor);
                lastIncrementFraction = 1 - (numIncrementsNeeded - numIncrementFactor);
                
                if lastIncrementFraction > obj.RASTER_DISP_STIM_RANGE_INCREMENT_FRACTION
                    numIncrementsNeeded = numIncrementsNeeded + 1;
                end
                
                ylim = [ylim(1) numIncrementsNeeded * obj.verticalRangeRasterInfIncrement];
            end
        end
        
        function [neural,analogmux,analogsolo,digwords] = zprvGetAvailAcqChans(obj)
            %Determine from teh SpikeGL assignments the set of acquisition
            %channel numbers, for each of the channel type groups
            %
            %Current limitations wrt SpikeGLX NI configuration
            % * Only 1 DAQ device supported
            % * Required that NI configuration satisfies (All MN) < (All MA) < (All MX)
            %     (In other words, an AUX channel interleaved between 2 banks of neural channels, is not allowed)
            
            muxFactor = obj.sglParamCache.niMuxFactor;
            
            neural = [];
            analogmux = [];
            analogsolo = [];
            nextchan = 0;
            
            %Extract the 4 types of chans supported through IMEC phase 2
            mn = str2num(num2str(obj.sglParamCache.niMNChans1)); %#ok<ST2NM>
            ma = str2num(num2str(obj.sglParamCache.niMAChans1)); %#ok<ST2NM>
            xa = str2num(num2str(obj.sglParamCache.niXAChans1)); %#ok<ST2NM>
            dw = str2num(num2str(obj.sglParamCache.niXDChans1)); %#ok<NASGU,ST2NM>
            
            %Determine the acquisition channel numbers
            for i=1:length(mn)
                neural = [neural ((i-1) *muxFactor) + (0:(muxFactor-1))]; %#ok<AGROW>
            end
            nextchan = neural(end) + 1;
            
            if ~isempty(ma)
                for i=1:length(ma)
                    analogmux = [analogmux nextchan + (i-1)*muxFactor + (0:(muxFactor-1))]; %#ok<AGROW>
                end
                nextchan = analogmux(end) + 1;
            end
            
            
            for i = 1:length(xa)
                analogsolo = [analogsolo nextchan + (i-1)]; %#ok<AGROW>
            end
            
            digwords = []; %Not supported (or used, anecdotally) at this time. need to understand line to channel mapping rules.
            
        end
        
        function zprvAssertAvailChansConstant(obj)
            [neural,analogmux,analogsolo] = zprvGetAvailAcqChans(obj);
            
            assert(isequal([neural analogmux analogsolo],[obj.neuralChansAvailable, obj.analogMuxChansAvailable, obj.analogSoloChansAvailable]),...
                'A change in the available neural and/or auxiliary channels has been detected. Changes to the SpikeGLX NI channel configuration is not allowed currently. Close & restart Spoke to employ new configuration.');
            
        end
        
        function zprvApplyChanOrderAndSubset(obj)
            
            if isempty(obj.sglParamCache.snsNiChanMapFile)
                fprintf('Channel Remapping: no snsNiChanMapFile defined in SpikeGLX, defaulting to standard mapping...\n');
                obj.neuralChanDispOrder = obj.neuralChansAvailable;
            else
                fprintf('Channel Remapping: using snsNiChanMapFile located at: %s\n',obj.sglParamCache.snsNiChanMapFile);
                obj.neuralChanDispOrder = parseChanMapFile(obj,obj.sglParamCache.snsNiChanMapFile);                
            end
            
            obj.sglChanSubset = GetSaveChansNi(obj.hSGL); %channel subset as specified in SpikeGLX. Wierd - this /has/ to be done here, outside of zprvZpplyChanOrderAndSubset() to avoid a hang.
            
            obj.neuralChanAcqList = intersect(obj.sglChanSubset,obj.neuralChansAvailable);
            obj.neuralChanDispList = obj.neuralChanDispOrder; %TODO: Apply subsetting to neuralChanDispOrder - ed: I thought subsets were displayed in strict channel order - what would change in the display due to subsetting?
            obj.auxChanProcList = [obj.analogMuxChansAvailable obj.analogSoloChansAvailable];
        end
        
        function zprvInitializeRasterGridLines(obj)
            colorOrder = get(0,'DefaultAxesColorOrder');
            
            obj.hRasterLines = {};
            
            for i = 1:numel(obj.stimEventTypes_)
                ppt = obj.PLOTS_PER_TAB;
                obj.hRasterLines{i} = gobjects(ppt,1);
                for j=1:ppt
                    obj.hRasterLines{i}(j) = animatedline('Parent',obj.hRasters(j),'LineStyle','none','Marker','d','MarkerSize',2,'LineStyle','none','Color',colorOrder(i,:));
                end
            end
        end

        function zprvSnapshot(obj)
            numPlots = obj.PLOTS_PER_TAB;
            
            for i=1:numPlots
                %Iterate through the plots and save to disk.
                %obj.hPlots(i) = axes('Parent',obj.hPanels.waveform(i),'Position',[0 0 1 1],'XLim',obj.horizontalRange); %,'YLim',obj.verticalRange);
                %obj.hRasters(i) = axes('Parent',obj.hPanels.raster(i),'Position',[0 0 1 1],'XLim',obj.horizontalRangeRaster);
                %obj.hPSTHs(i) = axes('Parent',obj.hPanels.psth(i),'Position',[0 0 1 1]);
                print(obj.hPlots(i),'BarPlot','-dpng')
            end
            disp('Saved plots to location: ');
            
        end
        
        function zprvPause(obj)
            disp('Paused Operation');
        end
    end
    
    
    %% ABSTRACT PROPERTY REALIZATIONS (most.Model)
    properties (Hidden,SetAccess=protected)
        mdlPropAttributes = zlclInitPropAttributes();
        mdlHeaderExcludeProps = {};
    end
    
    
end


%% LOCAL FUNCTIONS
function [newSpikeScanNums, maxNumWaveformsApplied] = zlclDetectSpikes(reducedData,fullDataBuffer,bufStartScanNum,postSpikeNumScans,thresholdVal,thresholdAbsolute,baselineMean,maxNumSpikes,sglChanSubset,chanSubset, horizontalRangeScans, debug, diagramm)
%Detect spikes from beginning in all but the spike-window-post time, imposing a 'refractory' period of the spike-window-post time after each detected spike
%
% reducedData: Cell array, one element per channel, containing data for each detected spike (from earlier timer callback period(s))
% fullDataBuffer: Array of recently acquired scans (samples)
% bufStartScanNum: Scan number of first element in fullDataBuffer
% postSpikeNumScans: Number of scans following each detected spike to exclude from spike detection (the spike detection 'refractory period')
% thresholdVal: May be a scalar, or a vector with one value per channel
% thresholdAbsolute: Logical. If true, both crossings above thresholdVal or below -thresholdVal are considered spikes.
% baselineMean: Mean value to subtract from data before detecting threshold crossings.
% maxNumSpikes: Scalar, indicating max number of spikes to detect per channel (from the start of the fullDataBuffer)
% chanSubset: The defined channel subset, which is specified by the user in spikeGLX. This is used to account for situations where the user-defined channels are a subset of the total number of channels.
% horizontalRangeScans: The negative and positive offsets that define the horizontal range of the display in samples.
% debug: set to true if debugging is activated. This should be equal to obj.debug.
%
% NOTES:
%  VI050812: Not clear that recentSpikeScanNums can ever be non-empty -- might be able to get rid of this logic (and reducedData argument) altogether?
maxNumWaveformsApplied = false;
numNeuralChans = length(reducedData);
newSpikeScanNums = cell(numNeuralChans,1);

spikesFoundPerChan = zeros(numNeuralChans,1);

if isscalar(thresholdVal)
    thresholdVal = repmat(thresholdVal,numNeuralChans,1);
end

if isscalar(baselineMean)
    baselineMean = repmat(baselineMean,numNeuralChans,1);
end

localspikes = cell(numel(chanSubset),1);

%for i=1:numNeuralChans
for h=1:numel(chanSubset)
    i = sglChanSubset(h)+1;
    
    %Determine recent (already detected) spike scan numbers to exclude from spike search
    lastSpikeScanNumIdx = find(reducedData{i}.scanNums < bufStartScanNum,1,'last');
    if isempty(lastSpikeScanNumIdx)
        recentSpikeScanNums = reducedData{i}.scanNums;
    else
        recentSpikeScanNums = reducedData{i}.scanNums(lastSpikeScanNumIdx + 1:end);
    end
    
    %Find new spikes one at a time, imposing refractory period
    spikesFound = 0;
    scansToSearch = size(fullDataBuffer,1) - postSpikeNumScans;
    
    %maxIdx = bufStartScanNum + scansToSearch;
    %Since the fullDatabuffer contains historic data, start detection at
    %the first point at which a theoretical spike at the beginning of the
    %detection window would have a full negative horizontal range's worth
    %of data. If this value is too small, then there won't be any data in a
    %portion of the user's display.
    
    currIdx = abs(horizontalRangeScans(1)); %Index into fullDataBuffer. 
    % NOTE: If currIdx is initialized to 1, it will begin spike detection at the very start of the fullDataBuffer.
    % Doing because it will detect as many spikes as possible. The downside is that it may detect a spike with no 
    % waveform data in the native part of the time window.
        
    while currIdx < scansToSearch
        %fprintf('currIdx: %d scansToSearch: %d postSpikeNumScans: %d\n',currIdx,scansToSearch,postSpikeNumScans);
        %Find at most one spike (threshold crossing) in the fullDataBuffer
        if thresholdAbsolute %Find crossings above or below absolute threshold level
            nextSpikeIdx = currIdx + find(diff(abs(fullDataBuffer(currIdx:scansToSearch,h) - baselineMean(i)) > abs(thresholdVal(i))) == 1,1);
        else
            if thresholdVal >= 0 %Find crossings above threshold level
                nextSpikeIdx = currIdx + find(diff((fullDataBuffer(currIdx:scansToSearch,h) - baselineMean(i)) > thresholdVal(i)) == 1,1); %Find at most one spike
            else %Find crossings below threshold level
                nextSpikeIdx = currIdx + find(diff((fullDataBuffer(currIdx:scansToSearch,h) - baselineMean(i)) < thresholdVal(i)) == 1,1); %Find at most one spike
            end
        end
        
        if isempty(nextSpikeIdx) %no spikes found in whole remainder of fullDataBuffer
            break;
        else
            spikesFound = spikesFound + 1;

            if spikesFound > maxNumSpikes
                maxNumWaveformsApplied = true;
                break;
            end
        end
        
        nextSpikeScanNum = bufStartScanNum + nextSpikeIdx - 1;

        %Add new spike, if not added already
        if ~ismember(nextSpikeScanNum,recentSpikeScanNums)
            %Note: fullDataBuffer is in "local" space. This means that it
            %is zero referenced to the beginning of the current "block" of
            %samples.
            %
            %bufStartScanNum is in spikeglx's "global space". This means
            %that a bufStartScanNum of 0 is the very beginning of the
            %spikeglx acquisition.
            %
            %This is why nextSpikeScanNum adds the nextSpikeIdx to the
            %bufStartScanNum. Think of bufStartScanNum as being the offset
            %in "Global Space" that equals zero in the "local space".
            newSpikeScanNums{i}(end+1) = nextSpikeScanNum;
            localspikes{h}(end+1) = nextSpikeIdx;
            %fprintf('bufStartScanNum: %d, nextSpikeScanNum: %d, size of fullDataBuffer: %s \n',bufStartScanNum, nextSpikeScanNum, sprintf('%d ',size(fullDataBuffer)));
            spikesFoundPerChan(i) = spikesFoundPerChan(i) + 1;
            
        end

        %Impose refractory period
        currIdx = nextSpikeIdx + postSpikeNumScans; %Will start with final scan of the post-spike-window...to use as first scan for next diff operation (first element never selected)
    end
    
end 
    if debug
        % *********************************************************************
        % ****************************** DEBUG ********************************
        % *********************************************************************
        buffersize = size(fullDataBuffer,1);

        if ~isempty(localspikes{1})
            if max(fullDataBuffer(1:buffersize,1)) > thresholdVal(1)
                tempbuffer = zeros(1,buffersize);
                figure(diagramm);
                cla;
                % *********************************************************************
                % boundschecking
                % *********************************************************************
                beginplot = (localspikes{1}(1)-500);
                endplot = (localspikes{1}(end)+500);
                for iter = 1:numel(localspikes{1})
                    spikelocation(iter) = localspikes{1}(iter);
                end
                fprintf('beginplot: %d, endplot: %d, spikeLocation: %d, localspikes: %s \n',beginplot,endplot,spikelocation,sprintf('%d ',localspikes{1}));

                if beginplot < 1
                    beginplot = 1;
                end
                if endplot > buffersize
                    endplot = buffersize;
                end
                % *********************************************************************
                % boundschecking
                % *********************************************************************

                tempbuffer(beginplot:endplot) = fullDataBuffer(beginplot:endplot,1);

                xlim manual;
                xlim([0 buffersize]);
                ylim manual;
                ylim([0,thresholdVal(1)+500]);
                hold on;

                text(10,4000,sprintf('beginplot: %d',beginplot));
                text(10,3300,sprintf('endplot: %d',endplot));
                text(10,2600,sprintf('spikelocation: %d',spikelocation));
                text(10,1900,sprintf('localspikes: %s',sprintf('%d ',localspikes{1})));
                plot(tempbuffer);
                plot(thresholdVal(1)*ones(1,buffersize));

                for iter = 1:numel(localspikes{1})
                    line([spikelocation(iter) spikelocation(iter)],get(gca, 'ylim'),'Color','red','LineStyle','--')
                end
                
                hold off;
            end
        end    
    end
    % *********************************************************************
    % ****************************** DEBUG ********************************    
    % *********************************************************************
end

function s = zlclInitPropAttributes()

s.running = struct('Classes','binaryflex','Attributes','scalar');
s.stimStartChannel = struct('Attributes',{{'integer' 'finite' 'nonnegative'}},'AllowEmpty',1);

s.displayMode   = struct('Options',{{'waveform' 'raster'}});
s.tabDisplayed  = struct('Attributes',{{'scalar' 'finite' 'positive' 'integer'}});

s.thresholdType = struct('Options',{{'volts' 'rmsMultiple'}});
s.thresholdVal = struct('Attributes',{{'scalar' 'nonempty' 'finite'}});
s.thresholdAbsolute = struct('Classes','binaryflex','Attributes','scalar');
s.baselineStatsRefreshPeriod = struct('Attributes',{{'scalar' 'positive' 'finite'}});

s.horizontalRange = struct('Attributes',{{'numel' 2 'finite'}});
s.verticalRange = struct('Attributes',{{'numel' 2 'finite'}});

s.waveformAmpUnits = struct('Options',{{'volts' 'rmsMultiple'}});

s.waveformsPerPlot = struct('Attributes',{{'scalar' 'finite' 'positive'}});
s.waveformsPerPlotClearMode = struct('Options',{{'all' 'oldest'}});
s.spikeRefractoryPeriod = struct('Attributes',{{'scalar' 'finite' 'nonnegative'}});

s.dataReadMode = struct('Options',{{'file' 'spikeGL'}});

s.filterWindow = struct('Attributes',{{'nonnegative' 'numel' 2}});
s.globalMeanSubtraction = struct('Classes','binaryflex','Attributes','scalar');

s.stimStartThreshold = struct('Attributes',{{'finite' 'scalar'}});
s.stimEventTypesDisplayed = struct();
s.horizontalRangeRaster = struct('Attributes',{{'numel' 2 'finite'}});
s.verticalRangeRaster = struct('Attributes',{{'numel' 2 'nonnegative'}});
s.verticalRangeRasterInfIncrement = struct('Attributes',{{'positive' 'scalar' 'finite'}});
s.stimEventClassifyFcn = struct();

s.refreshRate = struct('Attributes',{{'finite' 'positive' 'scalar'}});
s.refreshPeriodMaxWaveformRate = struct('Attributes',{{'scalar' 'positive'}});

s.psthTimeBin = struct('Attributes',{{'nonnegative' 'scalar' 'finite'}});
s.psthAmpRange = struct('Attributes',{{'nonnegative' 'finite' 'numel' 2}});
s.psthTimerPeriod = struct('Attributes',{{'positive' 'scalar'}});

end

function dfprintf(varargin)
fprintf(varargin{:});
end

%% Developer notes
%
% UNITS OF PROCESSING: Processing is done in volts/bit through most of
% pipeline. Conversion to volts done only in spots where needed: threshold
% determination (if in volts units) & plotting (if in volts units).

