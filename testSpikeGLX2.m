function [param, a, mat] = testSpikeGLX()

hSGL = SpikeGL('10.102.20.125');

if IsRunning(hSGL)
  numCalls = 100;
else
  numCalls = 1;
end

getTimes = zeros(numCalls,1);


subset = GetChannelSubset( hSGL );
subset = subset(1:128);

chans = [10,66,192,193,194,195];

%Initialize output args
param = struct();
a = 0;
mat = [];

for i=1:numCalls
  t0 = tic;

%   prm = struct();
%   prm.bandPass = 'true';
%   prm.autoStart = 'true';
%   prm.device = 'Dev2';
%   prm.rangeMin = -10;
%   prm.rangeMax = 10;
%   prm.clock = 'PFI2';
%   prm.channels = '0=192';
%   SetAOParams( hSGL, prm );
%   SetAOEnable( hSGL, 0 );
%   SetTrgEnable( hSGL, 1 );

  %param = GetParams(hSGL)
  %SetParams(hSGL, param)
  %count = GetScanCount(hSGL)
  %file = GetRunDir(hSGL);
  %IsConsoleHidden(hSGL)
  %Par2(hSGL,'v','G:/SGL_DATA/foo_g0_t0.bin');
  %VerifySha1( hSGL, 'G:/SGL_DATA/foo_g0_t0.bin' )
  %StopRun( hSGL );
  %StartRun( hSGL );
  %SetDigOut( hSGL, 0, 'Dev6/port0/line0' )
  %SetTrgEnable( hSGL, 0 );
  %IsSaving(hSGL)
  %SetRunName( hSGL, 'bar' )
  %file = EnumRunDir(hSGL)
  %a = GetChannelSubset( hSGL )
  %a = GetAcqChanCounts( hSGL );

  mat = GetLastNDAQData(hSGL, 4000, subset);
  getTimes(i) = toc(t0);

%  pause(10);
%  ConsoleHide(hSGL);
%  pause(10);
%  ConsoleShow(hSGL);

  if mod(i,50) == 0
    fprintf('Completed %d calls...Just extracted data matrix of size %s & class %s \n',i,mat2str(size(mat)),class(mat));
  end

   showdata( mat );

end

fprintf('Get time -- mean: %g ms\tstd: %g ms\n',1000*mean(getTimes),1000*std(getTimes));

%clear hSGL;

end

function showdata( mat )

showChanNum = 16;

x = 1:4000;
y = mat(:,showChanNum);
figure(1);
%set(gcf,'doublebuffer','on');
p = plot(x,y);
set(p,'XData',x,'YData',y);
drawnow;

end





