%Using this script on 2/24/12 with a quad dual-core Xeon, observed the following:
% No Acquisition: mean=17.6ms, std=0.2ms
% Acquiring 128 chans: mean=438ms, std=202ms
% Acquiring 128 chans, paused: mean=34ms, std=41ms
% Acquiring 128 chans,  running with logging disabled: mean=115ms, std=32ms
%
% Here's a dump of the GetParams result in the Acquiring 128 chan case:
% ans =
%
%                      acqMode: 4
%                    acqPDChan: 4
%             acqPDOffStopTime: 0.5
%          acqPDPassthruChanAO: -1
%                  acqPDThresh: 30405
%                 acqPDThreshW: 5
%              acqStartEndMode: 0
%        acqStartTimedDuration: 60
%           acqStartTimedImmed: 'false'
%           acqStartTimedIndef: 'false'
%            acqStartTimedTime: 0
%     aiBufferSizeCentiSeconds: 9
%                     aiString: '0:7'
%                 aiTermConfig: 10106
%                        aoDev: 'Dev2'
%                   aoPassthru: 'false'
%             aoPassthruString: '0=1,1=2'
%                   aoRangeMax: 10
%                   aoRangeMin: -10
%                      auxGain: 1
%                          dev: 'Dev1'
%                    doCtlChan: 0
%      doPreJuly2011IntanDemux: 'false'
%                     extClock: 'true'
%             fastSettleTimeMS: 15
%                   lowLatency: 'false'
%                   outputFile: 'C:/Documents and Settings/labadmin/My Documents/data_1.'
%                     rangeMax: 1.25
%                     rangeMin: -1.25
%              silenceBeforePD: 0.01
%                        srate: 20000
%             stimGlTrigResave: 'false'
%                 subsetString: 'ALL'
%               suppressGraphs: 'false'


hSGL = SpikeGL();

if IsAcquiring(hSGL)
  numCalls = 50;
else
  numCalls = 1000;
end

getTimes = zeros(numCalls,1);
for i=1:numCalls
  t0 = tic;
  param = GetParams(hSGL);
  getTimes(i) = toc(t0);
  
  if mod(i,50) == 0
    fprintf('Completed %d calls...\n',i);
  end
end

fprintf('GetParams time -- mean: %g ms\tstd: %g ms\n',1000*mean(getTimes),1000*std(getTimes));

%clear hSGL;