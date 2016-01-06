function spoke()
% This is spoke, an interface for acquiring data via SpikeGL

% make the spoke main controller object

hGrid = spoke.SpikeGrid('10.102.20.125');
hGridCtl = spoke.SpikeGridController(hGrid);

assignin('base','hGrid',hGrid);
assignin('base','hGridCtl',hGridCtl);

hGrid.initialize();
