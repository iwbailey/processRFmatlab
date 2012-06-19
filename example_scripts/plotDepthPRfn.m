function plotDepthPRfn

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
ISPRF = true; % look at P Rfns 

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
DIR='./prfns/';
dirlist=dir([DIR,'/TA*']);
files=dir( [DIR,'/',dirlist(1).name,'/*PRF.sac'] );
filename=[DIR, dirlist(1).name, '/',files(1).name];
disp(filename)

try 
  pRfn = readRFsacfile(filename);
catch ME
  error(ME.message)
end

% check the slowness units
if pRfn.rayp > 1 , 
  fprintf( 'Changing slowness from %f s/rad to %f s/km\n' , ...
	   [ pRfn.rayp, pRfn.rayp/6371] );
  pRfn.rayp = pRfn.rayp/6371;
end

% plot the ray path
figure(1);clf;
plotRayPaths( pRfn.rayp, regz, regVp, regVs );

% get the depth of the receiver function
[ z, rf2 ] = mapRF2depth( pRfn.time, pRfn.seis, pRfn.rayp,...
			  regz, regVp, regVs , ISPRF );

figure(2); clf;
subplot(1,2,1); plot(pRfn.seis,pRfn.time ); hold on;
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
[lon,lat,z] = sRayPath( pRfn.rayp, regz, regVs, ...
		      pRfn.stlo, pRfn.stla, pRfn.baz );

fprintf('Back azimuth: %f degrees\n', pRfn.baz );

% plot results
figure(3); clf;
plot3( lon+5*rf2, lat , z,'b'); hold on;
plot3( pRfn.stlo, pRfn.stla, 0 , 'rv' , 'MarkerSize', 10, ...
       'MarkerFaceColor', 'r' );
view(10,30);
daspect( [ 1 1 10])
grid on;
set(gca, 'ZDir', 'reverse');
