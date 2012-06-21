function testHP
%
%
% Try out a range of different High pass filters to see which is  best RF
%
%
clear;
format compact;
clf;


% Different high passes to investigate
FHP=[ 0.01, 0.02, 0.03, 0.04, 0.05]; 

% fixed parameters for processing
TBEFOREP=20; % how much of seismogram we require
TAFTERP=115;
TAPERW=0.05; % width of taper
DTOUT=0.1; % desired sampling
T0 = -5; % time limits for receiver function
T1 = 50;
TSHIFT = 5; % time to add on to beginning of rf
FLP=2; % low pass frequency (Hz)
F0 = 1.5; % gaussian width for filtering rfns
ITERMAX = 200;
MINDERR=1e-2; % stop iterating for decon at this change in error
PSEARCH=5; % number of seconds P can be earlier than reported
FORDER = 3; % order of butterworth filter 


fileprefix='TA.Q20A';

% Read in the data
try
  pSeis0 = readPseisSacfiles( ['./seismograms/',fileprefix], ...
				   TBEFOREP, TAFTERP );
  %clf; plotSeis(pSeis); 
catch
  disp('****Problem Reading file****')
  le=lasterror;
  disp(le.message)
end

% Loop through the frequencies
for hp = FHP,
  
  % Process the seismogram without plotting steps
  pSeis = processSeis( pSeis0, TAPERW, hp, FLP, FORDER, DTOUT, false);

  % Use this option for modified ligorria and ammon method
  [RFI3, RMS3] = makeRFitdecon(pSeis.seis(:,2), pSeis.seis(:,3) , ...
			       pSeis.delta, size(pSeis.seis,1) ,...
			       T0, T1, TSHIFT,...
			       F0, ITERMAX, MINDERR , 0);
  
  % Generate the Rfn structure
  idx = find(RMS3~=0);   % get final value of RMS

  pRfn = struct( 'seis', RFI3,...
		 'f0',F0,...
		 'rms',RMS3(idx(end)),...
		 'nit',length(idx),...
		 'time',T0-TSHIFT + pSeis.delta*(0:1:(length(RFI3)-1)) );

  % store the rf
  if( hp == FHP(1) ),
    pRfns = [ pRfn.seis ];
  else
    pRfns = [ pRfns; pRfn.seis ];
  end

end

% plot results

clf;
plot( pRfn.time, pRfns ); hold on;
axis tight;
legend( [num2str(FHP(1)),' Hz'], ...
	[num2str(FHP(2)),' Hz'], ...
	[num2str(FHP(3)),' Hz'],...
	[num2str(FHP(4)),' Hz'],...
	[num2str(FHP(5)),' Hz'] )
xlabel('Time (s)')  
ylabel('Amplitude (1/s)')
