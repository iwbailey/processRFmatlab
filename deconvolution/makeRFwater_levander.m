function [RF,RMS]=makeRFwater_levander(U,W,t0,dt,nt,wlevel,f0,wavelet,vb)

% Water level deconvolution
%
% [RF,RMS]=makeRFwater(U,W,t0,dt,nt,wlevel,f0,wavelet)
%
% Water level deconvolution to compute receiver function
% Uses ammons definition of the gaussian shaping filter
%
% OUT:
% RF = receiver function for correct number of samples
%
% IN:
% U = horizontal component (column form)
% W = vertical component  (column vector)
% t0 = time 
% dt = sample interval
% nt = number of samples
% wlevel = water level (proportion of max abs value of denomonator)
% f0 = gaussian width
% wavelet = 0 for Gaussian, 1 for Ricker wavelet
% vb: verbose

if( nargin < 9 ), vb = true; end


% 
% calculate receiver function
RF=zeros(nt,1);

% Fourier transform parameters
%
nfft = 2^nextpow2(nt);
df=1.0/(nfft*dt);
fny=1.0/(2.0*dt);
nfpts=nfft/2 + 1;

if( vb ),
  fprintf('nt: %i, nfpts: %i, nft: %i, fny: %.2f, df: %.3f\n',...
	  [nt, nfpts, nfft, fny, df]);
end

% fft
% uu = U(1:nt)';
% ww = W(1:nt)';
uu = U(1:nfpts)';
ww = W(1:nfpts)';
UFFT=fft(uu,nfft);
WFFT=fft(ww,nfft);

if wavelet == 1
    [wr2,C,tv,fv]=f_Ricker0(dt,nfft,f0);
elseif wavelet == 0
  C = gaussFilter( dt, nfft, f0 );
end

% Do water level operation on denominator
%
denom=WFFT.*conj(WFFT);
wlevel = max(abs(denom))*wlevel;
nwl = numel(denom( denom<wlevel ));
denom( denom<wlevel ) = wlevel;

if( vb ),
  fprintf('Dmax: %.3e, W. level: %.3e, Number corrected: %i/%i\n',...
	  [max(abs(denom)), wlevel, nwl,numel(denom)]);
end

% Numerator
numer = UFFT.*conj(WFFT);

% get receiver function
RFF=(C).*numer./denom;


% Get the phase
f=fftshift(-fny:df:fny-df);
phase = exp(-1i.*2.*pi.*f.*t0).';

% ifft
RFF = RFF.*phase';
RF2 = (real(ifft( RFF ,nfft))).';
RF(1:nt)= RF2(1:nt);

% rms fit
RMS = 1; % TODO
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wr2,C,tv,fv]=f_Ricker0(dt,nfft,f0)
%
% convolve with a source wavelet : Ricker
%
% Construct a 0-phase Ricker wavelet using the
% 
%
df=1.0/(nfft*dt);
fny=1.0/(2.0*dt);
%
% Ricker spectrum
%
C=zeros(1,nfft);

for k=2:nfft/2-1
    f=(k-1)*df;
    C(k)=(f/f0)^2*exp(1.-(f/f0)^2);
    C(nfft-k+2)=C(k);
end

%w1=fftshift(ifft(C));
%
%
tv= -nfft/2*dt:dt:(nfft-1)/2*dt;
fv=-fny:df:fny-df;
%
%
% wavelet
%
w2=ifft(C);
wr2=real(fftshift(w2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C]=f_Gaussian(dt,nfft,f0)

% [C]=f_Gaussian(dt,nfft,f0)
%
%  makes a gaussian spectrum of width f0
%  
% Langston 1979, JGR, uses G(omega)=exp(-omega^2/(4.*a^2))  
%  G(f0)=exp(-1/2)=.606
% 
%   a=sqrt(2)*pi*f0
%
% the formulation used here is somewhat broader band: rather than use
% exp(-1) for f0, use exp(-.5)
%  
%

df=1.0/(nfft*dt);
%a=sqrt(2)*pi*f0;
%a = 4*f0^2;
%a=1.5
    
nfft2=nfft/2;
    
f=zeros(1,nfft);
for ifreq=1:nfft2
  f(ifreq)=(ifreq-1)*df;
  f(nfft-ifreq+1)=-(ifreq*df);
end
        
C=zeros(1,nfft);
for ifreq=1:nfft2
  C(ifreq)=exp(-(2*pi*f(ifreq))^2/(4.*f0^2));
  C(nfft-ifreq+1)=exp(-(2*pi*f(nfft-ifreq+1))^2/(4.*f0^2));
end

C = C/dt;  % added to compare to iterative results

% clf;
% plot( C ) ; hold on;
% 
% w = 2*pi*f;
% plot( exp(-w.*w./(4.*f0*f0) )/dt, '--r');
% pause