% myobj = StartACQ( myobj )
% myobj = StartACQ( myobj, params )
% myobj = StartACQ( myobj, outfile )
%
%     Start an acquisition. The optional second argument,
%     params, is a struct of acquisition parameters as returned
%     from GetParams.m. Alternatively, the second argument can
%     be a string in which case it will denote the filename to
%     save data to, and the acquisition will implicitly use the
%     last set of parameters it has seen. If the second argument
%     is omitted then the last active parameters are used. If
%     acquisition was already in progress StartACQ will *not*
%     restart it. Rather, an error is thrown. If the params are
%     invalid, an error is thrown.
%
function [s] = StartACQ( varargin )

    s      = varargin{1};
    fname  = [];
    params = [];

    if( nargin > 1 )

        arg2 = varargin{2};

        if( ischar( arg2 ) )
            fname = arg2;
        elseif( isstruct( arg2 ) )
            params = arg2;
        else
            error( 'Invalid argument type for second argument. Must be a string or struct.' );
        end
    end

    if( IsAcquiring( s ) )
        error( 'Acquisition is already running.' );
    end

    if( isempty( params ) )
        params = GetParams( s );
    end

    if( ~isempty( fname ) )
        params.lastOutFile = fname;
    end

    % NB: Always send params as that does some important setup!
    % BK: Not sure that's true!!

    SetParams( s, params );

    DoSimpleCmd( s, sprintf( 'STARTACQ' ) );
end
