function modpath()
  % need to add the SpikeGL mex files to the path
  % this version of spoke is intended for use with SpikeGL v20111103
  % this version of spoke is intended for use in Matlab R2011b
  if ispc
    spikegl_matlab_path= ...
      ['C:/Documents and Settings/labadmin/My Documents/' ...
       'projects_janelia/spoke_related/SpikeGL_20111103/Matlab'];
  else
    spikegl_matlab_path= ...
      '~/projects_janelia/spoke_related/SpikeGL_2011103/Matlab';
  end    
  addpath(genpath(spikegl_matlab_path));
  % add toolbox to path
  %addpath(desvn(genpath('~/packages/taylor_matlab_toolbox/tags/release_1.10')));
  addpath(desvn(genpath('toolbox')));
end



function path_no_svn=desvn(path_raw)
  % eliminate .svn directories from a path string
  path_raw_as_array=split_path(path_raw);
  path_no_svn_as_array=cell(0,1);
  for i=1:length(path_raw_as_array)
    k=strfind(path_raw_as_array{i},'.svn');
    if isempty(k)
      path_no_svn_as_array{end+1}=path_raw_as_array{i};
    end
  end
  path_no_svn=combine_path(path_no_svn_as_array);
end

function path_as_array=split_path(path)
  % split a path on pathsep into a cell array of single dir names
  i_pathsep=strfind(path,pathsep);
  n=length(i_pathsep)+1;
  path_as_array=cell(n,1);
  if n>0
    if n==1
      path_as_array{1}=path;
    else
      % if here, n>=2
      path_as_array{1}=path(1:i_pathsep(1)-1);
      for i=2:(n-1)
        path_as_array{i}=path(i_pathsep(i-1):i_pathsep(i)-1);
      end
      path_as_array{n}=path(i_pathsep(n-1)+1:end);
    end
  end
end

function path=combine_path(path_as_array)
  % combine a cell array of dir names into a single path string
  n=length(path_as_array);
  if n>0
    path=path_as_array{1};
    for i=2:n
      path=[path pathsep path_as_array{i}];
    end
  end
end
