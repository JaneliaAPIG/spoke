% boolval = IsAcquiring( myobj )
%
%     Returns 1 if SpikeGL is currently acquiring data.
%
function [ret] = IsAcquiring( s )

    ret = sscanf( DoQueryCmd( s, 'ISACQ' ), '%d' );
end
