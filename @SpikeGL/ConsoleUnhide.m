% myobj = ConsoleUnhide( myobj )
%
%     Show SpikeGL console window.
%
function [s] = ConsoleUnhide( s )

    s = DoSimpleCmd( s, 'CONSOLEUNHIDE' );
end
