function event = stimEventClassifyFcn_model(stimData)
  % Model stimulus event-classification function
  % Function processes supplied stimData containing scans on the stimClassificationChannel immediately following stimulus detection. Processing classifies stimulus into one of the pre-defined event types.
  % When function is called with no arguments, it must return a structure with fields 'stimEventTypes' and 'stimEventy
  
  
  if nargin == 0
    event = struct();
    event.stimEventTypes = {'a' 'b'}; % (REQUIRED) String cell array of event names
    event.stimEventClassifyNumScans = 200; % (REQUIRED) Scalar integer specifying number of samples, following stimulus detection, used for stimulus event classification
    event.stimEventClassifyChannel = 58; % (OPTIONAL - Default=[]) Specify SpikeGL (demuxed) channel number of data to use for stimulus event classification. Must be an 'auxiliary' channel number. This channel's data is passed to function as stimData upon stimulus detection. If empty, the stimStartChannel specified in Spoke is used.
    return;
  end 
  
  if rand() > 0.5
    event = 'a';    
  else
    event = 'b';
  end
  
  
  
  
end

