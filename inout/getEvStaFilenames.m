function files = getEvStaFilenames( thisdir, suff1, suff2, suff3 )
%
% Find and match the three components of the filenames
%
% files = getEvStaFilenames( thisdir, suff1, suff2, suff3 )

% start with the first suffix
try,
  [fnames1, nf1] = getFilenames( thisdir, suff1 );
catch ME
    error('No files with suffix %s', suff1)
end
 
try,
  [fnames2, nf2] = getFilenames( thisdir, suff2 );
catch
      error('No files with suffix %s', suff2)
end

try,
  [fnames3, nf3] = getFilenames( thisdir, suff3 );
catch
  error('No files with suffix %s', suff3)
end

% print a warning if varying number of files
if( nf1 ~= nf2 | nf1 ~= nf3 | nf2 ~= nf3 ),
  fprintf('Warning: \n\tNumber of *%s files: ', suff1 );
  fprintf('%i\n', nf1);
  fprintf('\tNumber of *%s files: ', suff2 ); fprintf('%i\n', nf2);  
  fprintf('\tNumber of *%s files: ', suff3 ); fprintf('%i\n', nf3);  
  fprintf('This function will only use the prefixes for *%s files\n\n', suff1 );
end

% check files exist
if( nf1 == 0 ),
  error('No files with suffix %s', suff1)
end

% loop through all files
for i = 1:nf1,
  
  % get the prefix
  m = regexp( fnames1(i,:), suff1, 'split');
  fullpref = m{1}{1};
  
  % check the other files exist
  if( exist([fullpref,suff2], 'file') &&  exist([fullpref,suff3], 'file') ),
    files(i).name1 = fnames1{i};
    files(i).name2 = [fullpref,suff2];
    files(i).name3 = [fullpref,suff3];   
  end
end

return