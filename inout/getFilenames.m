function [fnames, nf] = getFilenames( DIR, SUFFIX)
% From a directory and filename suffix, get a list of files.  Uses
% unix find
%
% [fnames nf] = getFilenames( DIR, SUFFIX)
%
% IN:
% DIR = (string) directory to search
% SUFFIX = (string) suffix to look for in DIR
%
% OUT:
% fnames = (cell array) list of filenames such that fnames{1} gives the first one
% nf = number of filenames
%

% use the unix command to find the files
[a, fnames]=unix(['find ', DIR, ' -name \*', SUFFIX]);

% split based on new line
fnames = strread(fnames, '%s', 'delimiter', sprintf('\n'));

% get the number
nf = size(fnames,1);

return 
