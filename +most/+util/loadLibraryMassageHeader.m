function loadLibraryMassageHeader(massageFlags,infile,outfile)
%LOADLIBRARYMASSAGEHEADER Massage C/C++ header files into format required for use with loadlibrary function
%
%   massageFlags: String cell array of one or more 'flags' specifying massage options
%   infile: <OPTIONAL> Specify header file to be massaged. If none specified, user is prompted to select file.
%   outfile: <OPTIONAL> Specify name of massaged header file to generate. If none specified, the input file name is appended with '_MOD' and saved to same folder.
%
% Massage flags:
%   'includeHeader',XXX: 
%   'cppRemoveComments': 
%   'cppPassByRefToPointer':
%   'cppTypedefEnums':
%   'preproRemoveSpaces':

% NOTES
%   Many of the massage options are needed in particular for use under 64-bit Matlab;
%   Under 64-bits, Matlab's loadlibrary requires a strictly pure C header file
%
%   


assert(nargin >= 1 && iscellstr(massageFlags));

flagNames = massageFlags;

%Handle flags with args (includeHeader,...)
includeHeaderIdxs = find(strcmpi('includeHeader',massageFlags));
headersToInclude = {};
for i=1:length(includeHeaderIdxs)
    idx = includeHeaderIdxs(i);
    assert(length(massageFlags) > idx, 'Header filename must be provided after the ''includeHeader'' flag is used');
    headersToInclude{end+1} = massageFlags{idx+1};

    massageFlags(idx:idx+1) = []; %Will pass index only 
end
if ~isempty(includeHeaderIdxs)
    massageFlags{end+1} = 'includeHeader';
end

%Open source header file; read into a string
if nargin < 2 || isempty(infile)
    infile = uigetfile('*.h','Select Header file');    
else
    assert(ischar(infile) && exist(infile,'file'));           
end

inFid = fopen(infile,'r');
headerString = '';
while ~feof(inFid)
    headerString = [headerString fgets(inFid)]; %#ok<AGROW>
end
fclose(inFid);

%Generate massaged header string
for i=1:length(massageFlags)
    
    switch massageFlags{i}
        
        case 'includeHeader'
            for j=1:length(headersToInclude)                
                %TODO: Handle quoted and bracketed header specs
                
                %Prepend #include statement to mod header file
                headerString = sprintf('#include <%s>\n%s',headersToInclude{j},headerString);
            end
            
        case 'cppRemoveComments'
            %Removes all C++ sytle comments in mod header file
            headerString = regexprep(headerString,'/\*.*?\*/','');
            
        case 'cppPassByRefToPointer'
            %TODO: Maybe make this more restrictive, to avoid
            %possible false positives
            headerString = strrep(headerString,'&','*');
            
        case 'cppTypedefEnums'
            %Adds a #typedef after every enum definition
            headerString = regexprep(headerString,'enum\s*(?<enumName>\w*)(?<enumBlock>\s*\{.*?\};)','enum $<enumName>$<enumBlock>\ntypedef enum $<enumName> $<enumName>;');
            
        case 'preproRemoveSpaces'
            %TODO: Remove leading whitespace before c preprocessor macros -- the Matlab thunk compiler rejects this
            
        otherwise
            error('Unrecognized flag ''%s''', massageFlags{i});
    end
    
end

%Create output header file
if nargin < 3 || isempty(outfile)
    [p,f,e] = fileparts(infile);
    outfile = fullfile(p,[f '_MOD' e]);
else
    assert(ischar(outfile) && exist(outfile,'file'));
end

outFid = fopen(outfile,'w');
fprintf(outFid,'%s\n',headerString);
fclose(outFid);

end

