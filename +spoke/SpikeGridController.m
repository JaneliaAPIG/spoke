classdef SpikeGridController < most.Controller;
  %SPIKEGRIDCONTROLLER Summary of this class goes here
  %   Detailed explanation goes here
  
  %% ABSTRACT PROPERTY REALIZATIONS (most.Controller)
  properties (SetAccess=protected)
    propBindings = lclInitPropBindings();
  end
  
  
  %% VISIBLE PROPERTIES
  properties (AbortSet)
    filterEnable; %Logical specifying whether 
    filterWindowSpec = [300 8e3]; %Specification of filterWindow to use when filterEnable=true. When filterEnable=false, property is unchanged, while model filterWindow = [0 inf].             
    
  end
  
  properties (SetAccess=protected)
    stimEventClassifyFcnPath = ''; %Path of stimEventClassifyFcn last-specified via changeStimEventClassifyFcn
  end
  
  %% HIDDEN PROPERTIES
  

  
  properties (Hidden,Dependent)
    stimEventDisplayTypes; %Same as model stimEventTypes with 'all' appended to end
    
    stimEventDisplayAllTypes;
  end
  
  
  %% CONSTRUCTOR/DESTRUCTOR
  methods
    
    function obj = SpikeGridController(hModel)
      obj = obj@most.Controller(hModel,{'SpikeGrid'},{});  
      
      set(obj.hGUIs.SpikeGrid,'Name','Spoke Grid Control');      
            
    end     
    
    function initialize(obj)
      initialize@most.Controller(obj);
      
      figure(obj.hGUIs.SpikeGrid);
    end
    
  end
  
  %% PROPERTY ACCESS 
  
  methods
    
    function set.filterEnable(obj,val)
      validateattributes(val,{'logical'},{'scalar'});
      obj.filterEnable = val;
      
      %Side-effects
      obj.changeFilterSpec();
    end
    
    
    function set.filterWindowSpec(obj,val)
      validateattributes(val,{'numeric'},{'nonnegative' 'numel' 2});
      assert(val(2) > val(1), 'Property filterWindowSpec of class %s must be a 2 element array arranged in ascending order',mfilename('class'));
      obj.filterWindowSpec = val;      
      
      %Side-effects
      obj.changeFilterSpec();    
    end
    
      
    
    function val = get.stimEventDisplayTypes(obj)      
      
      if isempty(obj.hModel.stimEventTypes)
        val = {};
      else
        val = [obj.hModel.stimEventTypes {'all'}];
      end      
    end
    
        
    function val = get.stimEventDisplayAllTypes(obj)            
      val = iscell(obj.hModel.stimEventTypeDisplayed) && length(obj.hModel.stimEventTypeDisplayed) == length(obj.hModel.stimEventTypes);
    end
    
  end
  
  
  %% HIDDEN METHODS
  methods (Hidden)
    
    
    function changedDisplayMode(obj,~,~)      
            
      rasterControls = {'pbStimClassify' 'pmStimEventTypeDisplayed'};
      
      switch obj.hModel.displayMode
        case 'waveform'
          set(obj.hGUIData.SpikeGrid.rbWaveformDisplay,'Value',1);
          set(obj.hGUIData.SpikeGrid.rbRasterDisplay,'Value',0);
          
          set(obj.hGUIData.SpikeGrid.pmSpikeDisplayUnits,'Enable','on');
          set(obj.hGUIData.SpikeGrid.pmStimEventTypeDisplayed,'Enable','off');
          
        case 'raster'
          set(obj.hGUIData.SpikeGrid.rbWaveformDisplay,'Value',0);
          set(obj.hGUIData.SpikeGrid.rbRasterDisplay,'Value',1);
                    
          set(obj.hGUIData.SpikeGrid.pmSpikeDisplayUnits,'Enable','off');          
          
          %Enable StimEvenTypeDisplayed menu based on whether event types are being used or not
          obj.changedStimEventClassifyFcn();
      end      
      
      %Update pbStimClassify state
      obj.pbStimClassifyUpdate();
      
    end
    
    function changeFilterSpec(obj)
      
      if obj.filterEnable
        obj.hModel.filterWindow = obj.filterWindowSpec;
      else
        obj.hModel.filterWindow = [0 inf];
      end        
      
      set(obj.hGUIData.SpikeGrid.etFilterWindow,'String',mat2str(obj.filterWindowSpec));
      set(obj.hGUIData.SpikeGrid.cbFilterData,'Value',obj.filterEnable);      
    end
    
    function changedFilterWindow(obj,~,~)
      
      if isequal(obj.hModel.filterWindow(:),[0;inf])
        obj.filterEnable = false;
      else
        obj.filterWindowSpec = obj.hModel.filterWindow;
        obj.filterEnable = true;
      end                    
    end
    
    function changedMaxNumSpikesApplied(obj,~,~)
      
      if obj.hModel.maxNumSpikesApplied
        set(obj.hGUIData.SpikeGrid.etThresholdVal,'ForegroundColor',[.8 0 0]);
      else
        set(obj.hGUIData.SpikeGrid.etThresholdVal,'ForegroundColor',[0 0 0]);
      end
      
    end
    
    function changeTabDisplayed(obj,newVal)       
        try
           obj.hModel.tabDisplayed = newVal; 
        catch ME
          obj.changedTabDisplayed();
        end        
    end
    
    function changedTabDisplayed(obj,~,~)      
            
      for i=1:obj.hModel.MAX_NUM_TABS
        if i == obj.hModel.tabDisplayed
          set(obj.hGUIData.SpikeGrid.(sprintf('tbTab%d',i)),'Value',1);
        else
          set(obj.hGUIData.SpikeGrid.(sprintf('tbTab%d',i)),'Value',0);
        end
      end      
      
    end            
    
    function changeRunning(obj)      
      if obj.hModel.running
        obj.hModel.stop();
      else
        obj.hModel.start();
      end              
    end
    
    function changedRunning(obj,~,~)      
      
      if obj.hModel.running
        set(obj.hGUIData.SpikeGrid.pbStartOrStop,'String','Stop','BackgroundColor',[.8 0 0]);        
      else
        set(obj.hGUIData.SpikeGrid.pbStartOrStop,'String','Start','BackgroundColor',[0 .8 0]);                       
      end
      
      obj.pbStimClassifyUpdate();
      
    end
    
    function changeStimEventTypeDisplayed(obj,src)
      
      assert(~isempty(obj.hModel.stimEventClassifyFcn));
      
      idx = get(src,'Value');          
      eventType = obj.stimEventDisplayTypes{idx};
      if strcmpi(eventType,'all')      
        obj.hModel.stimEventTypeDisplayed = obj.hModel.stimEventTypes;
      else
        obj.hModel.stimEventTypeDisplayed = eventType;
      end
      
      %Update pbStimClassify state
      obj.pbStimClassifyUpdate();
      
    end
    
    function changedStimEventTypeDisplayed(obj,~,~)      
      
      if isempty(obj.hModel.stimEventTypeDisplayed)
        assert(isempty(obj.hModel.stimEventTypes));
        set(obj.hGUIData.SpikeGrid.pmStimEventTypeDisplayed,'Value',1); %should already be the case, but do it here anyway
      else                
        [tf,loc] = ismember(obj.hModel.stimEventTypeDisplayed,obj.stimEventDisplayTypes);                        
        assert(all(tf),'Unexpected value for SpikeGrid ''stimEventTypeDisplayed'' property'); %Not really needed, but done anyway...
        
        if isscalar(tf)       
          set(obj.hGUIData.SpikeGrid.pmStimEventTypeDisplayed,'Value',loc);
        else
          numEventTypes =  length(obj.hModel.stimEventTypes);
          assert(length(tf) == numEventTypes); %Not really needed, but done anyway..
          set(obj.hGUIData.SpikeGrid.pmStimEventTypeDisplayed,'Value',numEventTypes + 1); %Select 'all' option
        end
      end            
    end
    
    function pbStimClassifyPressed(obj)     
        
      %If button is available while running, must be in rastergram display-all-types mode
      if obj.hModel.running 
        obj.hModel.sortRastergram();
      else
        %Update model's stimEventClassifyFcn        
        
        [f,p] = uigetfile('*.m','Select Stim Classification Function',obj.stimEventClassifyFcnPath);
        
        if isnumeric(f) && f==0
          obj.hModel.stimEventClassifyFcn = [];
          obj.stimEventClassifyFcnPath = '';
        else
          currDir = cd();
          
          exception = [];
          try
            cd(p);
            [~,f] = fileparts(f);
            obj.hModel.stimEventClassifyFcn = str2func(f);
            obj.stimEventClassifyFcnPath = p;
          catch ME
            exception = ME;
          end
          
          cd(currDir);
          if ~isempty(exception)
            obj.hModel.stimEventClassifyFcn = [];
            obj.stimEventClassifyFcnPath = '';
            throw(exception);
          end
        end
      end
    end
    
    function pbStimClassifyUpdate(obj)
      
      rasterMode = strcmpi(obj.hModel.displayMode,'raster');
      hButton = obj.hGUIData.SpikeGrid.pbStimClassify;
      
      if obj.hModel.running        
        if rasterMode
          if obj.stimEventDisplayAllTypes
            set(hButton,'String','Stim Sort','Enable','on');
          else
            set(hButton,'String','Stim Classify...','Enable','off');
          end
        end
      else
        set(hButton,'String','Stim Classify...');
        if rasterMode
          set(hButton,'Enable','on');
        else
          set(hButton,'Enable','off');
        end
      end      
    end
    
    
    function changedStimEventClassifyFcn(obj,~,~)      
      if isempty(obj.hModel.stimEventTypes)
        set(obj.hGUIData.SpikeGrid.pmStimEventTypeDisplayed,'Enable','off','String',' ','Value',1);
      else
        set(obj.hGUIData.SpikeGrid.pmStimEventTypeDisplayed,'Enable','on','String',obj.stimEventDisplayTypes);
      end  
      obj.stimEventClassifyFcnPath = ''; %Don't store path generally -- only when stimEventClassifyFcn is set via this Controller class
    end
    
    
    function stepThresholdVal(obj,direction)
      %Direction: 1 to increment, -1 to decrement
      
      
      switch obj.hModel.thresholdType
        case 'volts' %step in millivolts
          if direction == 1
            obj.hModel.thresholdVal = obj.hModel.thresholdVal + .001; 
          else
            obj.hModel.thresholdVal = obj.hModel.thresholdVal - .001;
          end
          
        case 'rmsMultiple' %stop in multiple units
          if direction == 1
            obj.hModel.thresholdVal = obj.hModel.thresholdVal + 1;
          else
            obj.hModel.thresholdVal = obj.hModel.thresholdVal - 1;
          end      
      end
      
    end
      
    

      
    
    
    
  end
  
end

function s = lclInitPropBindings()
  
  s = struct();
  
  s.thresholdVal = struct('GuiIDs',{{'SpikeGrid','etThresholdVal'}});  
  s.globalMeanSubtraction = struct('GuiIDs',{{'SpikeGrid','cbGlobalMeanSubtract'}});  
  
  s.displayMode = struct('Callback','changedDisplayMode');
  s.filterWindow = struct('Callback','changedFilterWindow');  
  s.maxNumSpikesApplied = struct('Callback','changedMaxNumSpikesApplied');
  s.tabDisplayed = struct('Callback','changedTabDisplayed');    
  s.running = struct('Callback','changedRunning');
  
  s.stimEventClassifyFcn = struct('Callback','changedStimEventClassifyFcn');
  s.stimEventTypeDisplayed = struct('Callback','changedStimEventTypeDisplayed');
  
  s.thresholdType = struct('GuiIDs',{{'SpikeGrid','pmThresholdUnits'}});
  s.spikeAmpUnits = struct('GuiIDs',{{'SpikeGrid','pmSpikeDisplayUnits'}});
  
  s.psthTimeBin = struct('GuiIDs',{{'SpikeGrid','etTimeBinMs'}},'ViewScaling',1000);  
  s.psthAmpRange = struct('GuiIDs',{{'SpikeGrid','etPSTHAmpRange'}});

    
  %frequently-adjusted property table bindings
  s.spikeAmpWindow = struct('GuiIDs',{{'SpikeGrid','pcFreqAdjustProps'}},'PropControlData',struct('format','char'));
  s.spikeTimeWindow = struct('GuiIDs',{{'SpikeGrid','pcFreqAdjustProps'}},'PropControlData',struct('format','char'));
  s.spikesPerPlot = struct('GuiIDs',{{'SpikeGrid','pcFreqAdjustProps'}},'PropControlData',struct('format','numeric'));

  s.stimNumDisplayRange = struct('GuiIDs',{{'SpikeGrid','pcFreqAdjustProps'}},'PropControlData',struct('format','numeric'));

  s.refreshPeriodMaxSpikeRate = struct('GuiIDs',{{'SpikeGrid','pcFreqAdjustProps'}},'PropControlData',struct('format','numeric'));
  s.refreshRate = struct('GuiIDs',{{'SpikeGrid','pcFreqAdjustProps'}},'PropControlData',struct('format','numeric'));


  
  %rarely-adjusted property table bindings
  s.numAuxChans = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));
  s.gatingChannel = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));
  s.stimStartChannel = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));

  s.gatingDuration = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));
  s.gatingThreshold = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));
  
  s.stimStartThreshold = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));
  s.stimTimeWindow = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));
  s.stimNumDisplayRangeInfIncrement = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));

  s.spikesPerPlotClearMode = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','char'));
  
  s.thresholdRMSRefreshPeriod = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));
  s.thresholdRMSRefreshOnRetrigger =  struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));
  s.thresholdAbsolute = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','logical'));

  s.psthTimerPeriod = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','numeric'));


  %s.spikeAmpUnits = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','char')); 
  %s.stimEventTypeDisplayed = struct('GuiIDs',{{'SpikeGrid','pcRareAdjustProps'}},'PropControlData',struct('format','char'));

  
end

