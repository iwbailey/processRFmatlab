function plotDepthSRfn

clear;
format compact;
clf;

addpath ../ioFunctions/
addpath ../plotFunctions/
addpath ../depthFunctions/

% read a receiver function, plot the depth map

% define parameters
MAXZ = 400; % max depth in km
NZ = 401; % number of depth points to plot
ISPRF = false; % look at P Rfns 

% load a velocity model
vmod = '../velmodels1D/TNA.vsmod';
[hdr, x]=hdrload(vmod);
z = x(:,1);
vs = x(:,2);

% compute vp 
vp = 1.8.*vs;

% get reg vel model
[regz, regVp, regVs] = regVelmodel(  NZ, MAXZ, z, vp, vs  );

% load a receiver function
DIR='./srfns/';
dirlist=dir([DIR,'/TA*']);
files=dir( [DIR,'/',dirlist(1).name,'/*SRF.sac'] );
filename=[DIR, dirlist(1).name, '/',files(1).name];
disp(filename)

try 
  srfn = readRFsacfile(filename);
catch ME
  error(ME.message)
end

% check the slowness units
if srfn.rayp > 1 , 
  fprintf( 'Changing slowness from %f s/rad to %f s/km\n' , ...
	   [ srfn.rayp, srfn.rayp/6371] );
  srfn.rayp = srfn.rayp/6371;
end

% plot the ray path
figure(1);clf;
plotRayPaths( srfn.rayp, regz, regVp, regVs );

% get the depth of the receiver function
[ z, rf2 ] = mapRF2depth( srfn.time, srfn.seis, srfn.rayp,...
			  regz, regVp, regVs , ISPRF );

figure(2); clf;
subplot(1,2,1); plot(srfn.seis,srfn.time ); hold on;
set(gca, 'YDir', 'reverse');
set(gca, 'XTick',[]);
axis tight;
ylabel('Time (s)')

subplot(1,2,2); plot(rf2,z);
set(gca, 'YDir', 'reverse');
axis tight;
ylabel('Depth (km)');
set(gca, 'XTick',[]);

% get longitude and latitude
[lon,lat,z] = pRayPath( srfn.rayp, regz, regVp, ...
			srfn.stlo, srfn.stla, srfn.baz );

fprintf('Back azimuth: %f degrees\n', srfn.baz );

% plot results
figure(3); clf;
plot3( lon+5*rf2, lat , z,'b'); hold on;
plot3( srfn.stlo, srfn.stla, 0 , 'rv' , 'MarkerSize', 10, ...
       'MarkerFaceColor', 'r' );
view(10,30);
daspect( [ 1 1 10])
grid on;
set(gca, 'ZDir', 'reverse');
