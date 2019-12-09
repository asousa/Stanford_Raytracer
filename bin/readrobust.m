function x = readrobust(fid)
% Older versions of matlab don't handle inf or nan on an fscanf.  This
% works around it.
  x = [];
  while 1
    [xnew cnt] = fscanf(fid, '%g');  % read data
    x(end+1:end+cnt,1) = xnew;       % append
    if( feof(fid))                   % break out if end of file
      break
    end         
    str = fscanf(fid, '%s', 1);      % read one string
    if(strcmp(lower(str), 'inity'))  % Check for newer version of matlab
      % Newer versions of matlab will properly handle +-inf but not
      % +-infinity.  Do nothing and skip past it.
    elseif(strcmp(lower(str), 'nan'))    % if is a NaN
      x(end+1,1) = NaN;              %   append
    elseif(strcmp(lower(str), 'inf'))
      x(end+1,1) = inf;              %   append
    elseif(strcmp(lower(str), 'infinity'))
      x(end+1,1) = inf;              %   append
    elseif(strcmp(lower(str), '+inf'))
      x(end+1,1) = inf;              %   append
    elseif(strcmp(lower(str), '+infinity'))
      x(end+1,1) = inf;              %   append
    elseif(strcmp(lower(str), '-inf'))
      x(end+1,1) = -inf;              %   append
    elseif(strcmp(lower(str), '-infinity'))
      x(end+1,1) = -inf;              %   append
    else
      % UNRECOGNIZED
      disp(['Unrecognized string "' str '" in readrobust.m.  Replacing with NaN']);
      x(end+1,1) = NaN;
    end
  end
  
