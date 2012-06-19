function compare_water
%
% compare the water level deconvolution methods using the example in Ammon's code
format compact;

addpath '../'  % decon functions
addpath '../../ioFunctions/' % for file in 
addpath '../../plotFunctions/' % for plotting input data 
addpath '../../sigprocFunctions/' % for plotting input data 

% read data
[Zin,hdrZ]=sacsun2mat('test_data/uln_1995_113_b_057_d_045.z');
[Ein,hdrE]=sacsun2mat('test_data/uln_1995_113_b_057_d_045.e');
[Nin,hdrN]=sacsun2mat('test_data/uln_1995_113_b_057_d_045.n');

% get time axis
t0 = hdrZ.times.b;
dt = hdrZ.times.delta;
nt = hdrZ.trcLen;
baz = hdrZ.evsta.baz;
time = t0 + dt*(0:1:(nt-1));
% hdrZ.times

% % plot data and wait for user input
% figure(1);
% clf;
% plot3seis( time, Zin, time, Ein, time, Nin )

% tmp=input('prompt after plotting components.');

% rotate to ZRT coordinates
seis = rotateSeisENZtoTRZ( [Ein, Nin, Zin] , baz );
rseis = seis(:,2);
zseis = seis(:,3);

% Receiver function parameters from Ammon's readme file
tdel = 30; 
f0 = 1.0; % pulse width
wlevel = 0.05; % can't tell from readme, chosen based on best fit

%--------------------
% Make receiver function
fprintf('\nAmmon method\n')
[rf1,rms] = makeRFwater_ammon( rseis, zseis, tdel, dt, nt, wlevel, f0, true);
% get time for RF
time = - tdel  + dt*(0:1:numel(rf1)-1);

fprintf('RMS: %.2f\n', rms)

t=cputime;
for i=1:1e2,
  makeRFwater_ammon( rseis, zseis, tdel, dt, nt, wlevel, f0, false); 
end
timeTakenFor100 = cputime-t


%--------------------
% levander method
fprintf('\nLevander method\n')
wvt = 0;%wavelet = 0 for Gaussian, 1 for Ricker wavelet
[rf2,tmp] = makeRFwater_levander( rseis, zseis, -tdel, dt, nt, wlevel, f0, wvt, true);

t=cputime;
for i=1:1e3,
  makeRFwater_levander( rseis, zseis, -tdel, dt, nt, wlevel, f0, wvt, false);
end
timeTakenFor100 = cputime-t

fprintf('RMS: %.2f\n', tmp)

%--------------------
% Adjusted Ammon method
fprintf('\nAdjusted Ammon method\n')
[rf3,rms] = makeRFwater( rseis, zseis, tdel, dt, nt, wlevel, f0, true);

t=cputime;
for i=1:1e2,
  makeRFwater( rseis, zseis, tdel, dt, nt, wlevel, f0, false); 
end
timeTakenFor100 = cputime-t

fprintf('RMS: %.2f\n', rms)

%--------------------
% plot results
clf;
h1 = plot(time,rf1,'k'); hold on;
h2 = plot(time,rf2,'--m'); hold on;
h3 = plot(time,rf3,'--c'); hold on;

legend([ h1, h2, h3 ], ...
       'Ammon method',...
       'Levander method',...
       'Adjusted Ammon method')

disp('The difference in the result is because the ammon method uses more of the signal')
disp('Levander sets anything that could wraparound to zero.')
disp('Adjusted ammon method pads with extra zeros to avoid wraparound.')

axis tight

xlabel('Time (s)')

