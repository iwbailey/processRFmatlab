function tarr = getArrTime( arrname, times, labels )

% Get the time of an arrival with label arrname from the arrival times and label arrays
%
% tarr = getarrtime( arrname, times, labels )
% 
% IN: 
% arrname = label of phase looking for (e.g., 'P')
% times = array of all times to look through
% labels = labels corresponding to times
%
% OUT:
% tarr = time of desired arrival
%
tarr = NaN;

for i = 1:length( times ),
  % match the correct label 
  if(strcmp( strtrim(labels(i,:)), arrname) == 1 ),
    % set to corresponding time
    tarr = times(i);
  end
end

% Check if found
if( isnan( tarr ) ),
  for i = 1:length( times ),
    disp( labels(i,:) )
  end
  error( ['ERROR: time of ',arrname,' arrival not defined'] );
end

return


%---------------------------------------------------------
