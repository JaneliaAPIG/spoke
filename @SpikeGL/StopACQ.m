% myobj = StopACQ( myobj )
%
%     Unconditionally stop the currently running acquisition (if
%     any), closing all data files immediately and returning the
%     app to the idle state. See IsAcquiring.m to determine if
%     an acquisition is running.
%
function [s] = StopACQ( s )

    DoSimpleCmd( s, 'STOPACQ' );
end
