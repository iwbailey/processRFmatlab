function [RFI,RMS]=makeRFitdecon_levander(UIN,WIN,t0,dt,nt,f0,itermax)
% Levanders slightly different code for calculating the iterative decon
%
% [RFI,RMS]=makeRFitdecon(UIN,WIN,t0,dt,nt,f0,itermax)
%
% UIN = radial component 
% WIN = vertical component
% t0 = start time
% dt = sample rate
% nt = number of samples
% f0 = width of filter
% itermax = max # iterations
%
% RFI = receiver function
% RMS = Root mean square error after each iteration
%

%RFI=zeros(nt);

nfft = 2^nextpow2(nt);

%[nfft,f,fny,df]=fftparms2(nt,dt);
nfft2=nfft/2;

i_t0=round(abs(t0/dt))+1;

%t=t0:dt:(t0+(nt-1)*dt);

% make shaping filter positioned at t0
ntw = 2*round(1/(dt*f0));
nt2 = round(1/(dt*f0));
C = f_GaussianA(dt,nfft,f0); % filter in freq domain

wr2 = ifft(C,nfft); % filter in time domain
wr2 = fftshift(wr2); % shift zero component to 0

% shift zero to correct location
g = zeros(nt,1);
g(i_t0+1-ntw:i_t0+1+ntw) = wr2(nfft2+1-ntw:nfft2+1+ntw);

% normalize so that integral g is 1: definition of the delta
g = g./sum(g);  

% %U0=UIN;
% U=UIN;
% Z=WIN;
    
% put z component in freq domain and filter, then back in time domain
IZ=ifft( C.*fft(WIN,nfft) , nfft);
Z=IZ(1:nt);

% same for horizontal
IU=ifft( C.*fft(UIN,nfft), nfft);
U=IU(1:nt);

U0=U;

% IWB: Get power in numerator for error scaling
powerU = sum(U.^2);

% use U as the first estimate of R
R=U;

RMS=zeros(itermax-1,1);

for iter=1:itermax
%       clf;
%      subplot(3,1,1); plot(R,'-k');
%      subplot(3,1,2); plot(Z,'-k');
%      subplot(3,1,3); plot(U,'-k');
%      tmp = input('prompt');
    
    RZt=conv(R,Z); % compute predicted horiz component
    
    RZ(1:nt,1)=RZt(i_t0+1:nt+i_t0); % get part after t0

    %  calculate error between RZ and U0
    %RMS(iter)=sqrt(mean((U0-RZ).^2));
    if( iter > 1 ),
      RMS(iter-1)=sum((U0-RZ).^2)/powerU; % change to ligorria & ammon definition
    end
    % Cross correlate the radial with the vertical
    UZ=xcorr(U,Z);
    

    % get the lag between two biggest spikes
    [m1,i1]=max(abs(UZ(nt:length(UZ))));
    %    [m1,i1]=max(abs(URZ));
    if( i1==1 && iter > 1),
        [m1,i1]=max(abs(UZ(nt+nt2:length(UZ))));
        i1=i1+nt2;
    end

    i1=i1+nt-1;
    
    m1=sign(UZ(i1))*m1; % get absolute value, IWB: added /dt
    
    % Get location and time of max spike in the z component
    ZZ=xcorr(Z,Z);
    [m2,i2]=max(ZZ);

%     clf;
%     subplot(4,1,1); plot(U,'-k');
%     subplot(4,1,2); plot(UZ,'-k'); hold on;
%     plot(i1,m1,'rx');
%     subplot(4,1,3); plot(Z,'-k');
%     subplot(4,1,4); plot(ZZ,'-k'); hold on;
%     plot(i2,m2,'rx');
%     tmp = input('prompt');
   
    M1=m1/m2; % get amplitude of the spike
    ishift=i1-i2+1;  % lag of spike

    % Shift spike according to the lag
    gshift=zeros(nt,1);
    gshift(ishift:nt)=g(1:(nt-ishift+1));

    % Update Receiver function
    if( iter==1 ),
        R=M1*gshift;
    else
        R=M1*gshift+R;
    end
    
    % Remove the spike from the horizontal
    %Utemp=conv(R,Z); % IWB: added *dt
    Utemp=conv(R,WIN(1:nt)); % IWB: added *dt
    Uest(1:nt,1)=Utemp(i_t0+1:nt+i_t0);
    U=U0-Uest;
      
end

RFI=R/dt;

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C=f_GaussianA(dt,nfft,f0)
%
%  makes a gaussian spectrum of width f0
%  
% Liguria & Ammon, use G(omega)=exp(-omega^2/(4.*L^2))  
%  G(f0)=exp(-1/2)=.606
% 
%   L=f0^2*4
%
% the formulation used here is somewhat broader band: rather than use
% exp(-1) for f0, use exp(-.5)
%  
%

df = 1./(nfft*dt);
a = 4*f0^2;
    
nfft2 = 0.5*nfft;

%f = zeros(nfft,1);
f = ones(nfft,1);

for ifreq=1:nfft2,
    f(ifreq)=(ifreq-1)*df;
    f(nfft-ifreq+1)=-(ifreq*df);
end
        
C=zeros(nfft,1);
    
for ifreq=1:nfft2
    C(ifreq)=exp(-(2*pi*f(ifreq))^2/a);
    C(nfft-ifreq+1)=exp(-(2*pi*f(nfft-ifreq+1))^2/a);
end
    
C(nfft2+1)=exp(-(2*pi*f(nfft2+1))^2/a);

%C = C./dt; % added for comparison with iterative results