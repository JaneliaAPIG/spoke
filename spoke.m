function spoke()
% This is spoke, an interface for acquiring data via SpikeGL

% make the spoke main controller object

if nargin == 0
   [~,s] = system('ipconfig /all');
   
   out = regexp(s,'IPv4 Address[^0-9]*([.0-9]*)','tokens');
   assert(numel(out) > 0, 'Unable to determine the IP address for this machine');
   ipAddress = out{1}{1};
   sprintf('Detected local IP Address %s. Connecting to the SpikeGL Remote Connection server.',ipAddress);
end    
   
hGrid = spoke.SpikeGrid(ipAddress);
hGridCtl = spoke.SpikeGridController(hGrid);

assignin('base','hGrid',hGrid);
assignin('base','hGridCtl',hGridCtl);

hGrid.initialize();
