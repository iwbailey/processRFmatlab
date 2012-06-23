function checkPRfn
% Open all SAC files containing a receiver functions in a
% directory. Check whether it falls into a number of criteria, then
% decide whether or not to keep it.

%--- checkPRfn.m ---
%
% Filename: checkPRfn.m
%
% Author: Iain W Bailey
% Maintainer:
% Created: Wed Oct 20 15:16:39 2010 (-0700)
% Version: 1
% Last-Updated: Mon Jun 18 22:29:14 2012 (-0400)
%           By: Iain W. Bailey
%     Update #: 140
%
%- Change Log:
%
% 18 June 2012 Removed adding paths in the first part of the
% script. This is now done by a script in the root directory. /
% Made some commands more windows friendly by taking out unix stuff
% /
%
%- Code:


DIR = fullfile('prfns','prfns_iter_2.50'); % directory for in/output data

isAuto = false;

if( isAuto ),
    fprintf(['Automatically checking receiver function quality in directory ' ...
        '%s\n'], DIR);
else
    fprintf(['Manually checking receiver function quality in directory ' ...
        '%s\n'], DIR);
end

% make a log of all the outputs
logfile = fullfile( DIR,'prf_check.log' );
if( exist( logfile, 'file' ) == 2 ),
  delete( logfile );
end
fprintf('Saving terminal output to %s\n', logfile)
diary(logfile);

% PARAMETERS
MAXRMS = 0.3 % Don't allow any RMS values greater than this
PSEARCH = 5 % Largest spike must be +ve within this many seconds of zero
MAXSPIKE2 = 0.5 % 2nd max spike cannot be > this much times max spike
MAXPOWT = 0.004 % maximum power of transverse receiver function

STATION_PREFIX='TA';
PSUFFIX='PRF.sac';
TSUFFIX='TRF.sac';

% directories for unwanted RFs
REJDIR = fullfile( DIR, 'auto_rejected/');
if( ~exist(REJDIR,'dir') ),
  mkdir(REJDIR);
end
fprintf('Putting bad receiver functions in %s\n',REJDIR);

% get station directories
dirlist = dir( fullfile( DIR, [STATION_PREFIX,'_*']) );
fprintf('Number of station directories: %i\n', numel(dirlist));

% loop through directories
for di = 1:numel(dirlist),

    fprintf('DIR # = %i\n',di);

    % P Rfn file list
    files1 = dir( fullfile( DIR, dirlist(di).name,['*.',PSUFFIX] ));

    % T Rfn file list
    files2 = dir( fullfile( DIR, dirlist(di).name,['*.',TSUFFIX] ));

    fprintf('Number of PRF files: %i\n', numel(files1));
    fprintf('Number of TRF files: %i\n', numel(files2));

    % loop through files in directory
    for fi = 1:numel(files1),

        fprintf('File # = %i\n',fi);

        % Open files
        try
            filename = fullfile( DIR, dirlist(di).name, files1(fi).name );
            [t, seis, hdr] = sac2mat(filename);
        catch ME
            error(ME.message)
            continue;
        end

        % Flag for checking whether we keep the receiver function
        isOk = true;
        
        % Get and check the RMS from the deconvolution
        rms = hdr.user.data(7);
        fprintf('\tDeconvolution RMS = %f\n',rms);
        if( rms > MAXRMS ),
            fprintf( ['\tLarger than defined threshold, indicating ',...
                'deconvolution didnt work\n'] );
            isOk=false;            
            if( isAuto ), fprintf('\t\tREJECTED\n'); end
        else

        end

        % Get and check the maximum peak
        [ a, ~, tmax ] = getAbsMaxSeisTime( seis, t );
        fprintf('\tMax abs peak at t=%.3f\n',tmax);
        if( abs(tmax) > PSEARCH ),
            fprintf('\tGreater than threshold indicating no strong p-wave signal\n')
            fprintf( 'Note this logic only works for R/Z receiver functions\n');
            isOk=false;
            if( isAuto ), fprintf('\t\tREJECTED\n'); end
        end

        % check size of maximum peak isn't too large
        printf('\tMax amplitude = %f\n',a);
        if( a > 1 ),
            fprintf('\tGreater than 1, indicating deconvolution didn''t work\n');
            isOk=false;
            if( isAuto ), fprintf('\t\tREJECTED\n'); end            
        end

        % get size of second largest spike
        a2 = max( getAbsMaxSeisTime( pRfn.seis, pRfn.time , PSEARCH ), ...
                  getAbsMaxSeisTime( pRfn.seis, pRfn.time , pRfn.time(1), ...
                  -PSEARCH) );
        if( isAuto && abs(a2)/a > MAXSPIKE2 ),
            fprintf(['\tRejected because 2nd largest spike > %f x max spike', ...
                ' indicating another phase or noise is contaminating the ',...
                'signal\n'], ...
                    MAXSPIKE2 )
            isOk=false;
        else
            fprintf('\tSize of the second largest spike: %f\n',a2)
        end

        % Plot the receiver function
            clf;
    p1 = plot( rftime, rfseis, '-b', 'linewidth', 2 ); hold on;

        if( ~isAuto ),

            % Review the plot
            % plot colour according to above criteria
            if( isok )
                plot2Rfn( pRfn.time, pRfn.seis, tRfn.seis, ...
                          [pRfn.kstnm,' P'], 'T' )
            else
                plot2Rfn( pRfn.time, pRfn.seis, tRfn.seis, ...
                          [pRfn.kstnm,' P'], 'T' ,'-r', '-b' )
            end

            % alter axis dimensions to weed out excessive RFs
            axis([ min(pRfn.time), max(pRfn.time), -0.3, 0.5 ])

            % display parameters to screen
            fprintf('p = %f\nbaz=%f\n',[pRfn.rayp,pRfn.baz])

            % get user input
            display('Keep with L click.  Middle click to go back. Reject with R click...');
            [xtmp,ytmp,keepoption] = ginput(1);
            
            % message
            if( keepoption == 1), disp('...keeping');
            else disp('...rejecting');
            end

        else

            % keep or reject based on above criteria
      if( isok  ),
	plot2Rfn( pRfn.time, pRfn.seis, tRfn.seis, ...
		  [pRfn.kstnm,' P'], 'T' )
	keepoption = 1;
      	disp('...keeping')
	if( isPause ), pause(0.5); end
      else
	plot2Rfn( pRfn.time, pRfn.seis, tRfn.seis, ...
		  [pRfn.kstnm,' P'], 'T' ,'-r', '-b' )
	keepoption = 3;
	disp('...rejecting');
	if( isPause ), pause(0.5); end
      end
    end

    % move file out of directory, this won't always work
    if( keepoption == 3 ),
      REJDIR2 = [ REJDIR,'/',dirlist(di).name];
      if( exist(REJDIR2,'dir') ~=7 ),
	unix(['mkdir ',REJDIR2]);
      end
      unix(['mv ',filename,' ',REJDIR2]);
      unix(['mv ',filename2,' ',REJDIR2]);
     end

    if( keepoption == 2 ),
      disp('...Reversing')
      fi = fi - 2
      if( fi < 0 ),
	di = di - 2
	fi = length(files1) + fi
	if( di < 1 ),
	  di = 1
	  fi = 0
	end
      end
    end

    fi=fi+1;

  end % END file loop
  di=di+1;
end % end directory loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- checkPRfn.m ends here
