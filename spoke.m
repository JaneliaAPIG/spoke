function spoke(ipAddress)
% This is spoke, an interface for acquiring data via SpikeGL

% Fix for freezing of MATLAB when running in hardware OpenGL mode.
% TODO: Fix Spoke to work with Hardware Rendering!!
opengl('save','software');

% make the spoke main controller object

if nargin == 0
   [~,s] = system('ipconfig /all');
   
   % Check IPv4 address of this machine. Note, this only works on newer machines.
   out = regexp(s,'IPv4 Address[^0-9]*([.0-9]*)','tokens');
   if not(numel(out) > 0)
       % In case IPv4 address check failed, try a different regexp compatible with older machines.
       out = regexp(s,'IP Address[^0-9]*([.0-9]*)','tokens');
   end
   assert(numel(out) > 0, 'Unable to determine the IP address for this machine');
   ipAddress = out{1}{1};
   sprintf('Detected local IP Address %s. Connecting to the SpikeGL Remote Connection server.',ipAddress);
end    
   
hGrid = spoke.SpikeGrid(ipAddress);
hGridCtl = spoke.SpikeGridController(hGrid);

assignin('base','hGrid',hGrid);
assignin('base','hGridCtl',hGridCtl);

set(hGridCtl.hGUIsArray,'CloseRequestFcn',@(src,evnt)hGrid.delete); % For the momnet GUIs owned by Controller are never really killed

hGrid.initialize();
