%Using this function on 2/24/12 with a quad dual-core Xeon, observed the following:
% useMatlabAPI = false
%   128 chans, display paused: readTime mean=15ms, std=21ms; getScanNumTime mean=.05ms, std=.02ms (but get varying readTime readings on successive 150 read measures)
%   128 chans, display active: mean=13ms, std=16ms (i.e. display does not slow down reading in any way)
%
%
%
% useMatlabAPI = true, normal drive
%   128 chans, display paused: mean=1034ms, std=56ms
%
% useMatlabAPI = true, RAM drive (64MB, temp-file size = 16MB)
%   128 chans, display paused: readTime mean=70.9ms, std=95ms; getScanNumTime mean=69ms, std=81ms  (typical time is 20ms, but periodically much longer)



function testDataRead(useMatlabAPI)
  
  if nargin < 1
    useMatlabAPI = false;
  end
  
  if useMatlabAPI && isempty(strfind(getenv('tmp'),'r:'))
    numReads = 5;
  else
    numReads = 150;
  end
  
  readTimeSecs = .2;
  bytesPerSample = 2;
  dispTimeSecs = 1;
  
  dispTimeFactor = ceil(dispTimeSecs/readTimeSecs);
  
  hSGL = SpikeGL();
  
  %Read needed params
  params = GetParams(hSGL);
  scanRate = params.srate;
  numChans = length(GetChannelSubset(hSGL));
  
  numScansRead = 0;
  
  scansPerPeriod = scanRate * readTimeSecs;
  
  if ~useMatlabAPI
    
    assert(IsSaving(hSGL) > 0,'SpikeGL must be saving data to disk if not using Matlab data API');
    
    fileName = params.outputFile;
    
    fid = fopen(fileName,'rb');
    assert(fid > 0,'Unable to open SpikeGL data file');
    
  end
  
  readTimes = nan(numReads,1);
  getScanNumTimes = nan(numReads,1);
  
  delete(timerfindall);
  hTimer = timer('Name','testSpikeGLMatlabDataAPI','ExecutionMode','fixedRate','TimerFcn',@timerFcn,'StopFcn',@stopFcn,'Period',readTimeSecs,'StartDelay',readTimeSecs,'TasksToExecute',numReads);
  start(hTimer);
  
  function timerFcn(src,evnt)
    
    readIdx = get(src,'TasksExecuted');
    
    
    if useMatlabAPI
      
      tic;
      totalNumScans = GetScanCount(hSGL);
      getScanNumTimes(readIdx) = toc();
      
    else
      tic;
      d = dir(fileName);
      numBytes = d.bytes;
      totalNumScans = floor(numBytes/numChans/bytesPerSample);
      getScanNumTimes(readIdx) = toc();
    end
    
    scansAvailable = totalNumScans - numScansRead;
    assert(scansAvailable >= scansPerPeriod,'Expected at least %d scans would be available to read, but only found %d scans\n',scansPerPeriod,scansAvailable);
    
    if useMatlabAPI
      tic;
      newData = GetDAQData(hSGL,totalNumScans-scansPerPeriod,scansPerPeriod); %Read most recent scans
      readTimes(readIdx) = toc();
    else
      tic;
      fseek(fid,numScansRead * bytesPerSample * numChans,'bof');
      [newData,cnt] = fread(fid,scansPerPeriod * numChans,'int16=>int16');
      readTimes(readIdx) = toc();
      
      assert(cnt == scansPerPeriod * numChans,'Expected to read %d samples, but instead read %d',scansPerPeriod * numChans,cnt);
    end
    
    if mod(readIdx,dispTimeFactor) == 0
      fprintf('Read %d scans...\n',scansPerPeriod);
    end
    %newData = reshape(newData,numChans,[])'; %Create into matrix with numChans columns, and rows representing scans
    
    
  end
  
  function stopFcn(src,evnt)
    
    readTimes(isnan(readTimes)) = [];
    getScanNumTimes(isnan(getScanNumTimes)) = [];
    
    fprintf('Read time statistics -- mean: %g ms\tstd: %g ms\n',1000*mean(readTimes),1000*std(readTimes));
    fprintf('Get scan num time statistics -- mean: %g ms\tstd: %g ms\n',1000*mean(getScanNumTimes),1000*std(getScanNumTimes));
    
  end
  
end