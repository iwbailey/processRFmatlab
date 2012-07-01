% function plotDepthPRfn

% read a receiver function, plot the depth map

% define parameters
MAXZ = 400; % max depth in km
DZ = 1;
ISPRF = true; % look at P Rfns
SUFFIX='PRF.sac';

% load a velocity model
[z, vs] = tna();

% compute vp
vp = 1.8.*vs;

% Remove discontinuities for later interpolation
z(diff(z)==0) = z(diff(z)==0) - 1e-8;

% Get at regular intervals
zReg = 0:DZ:MAXZ;
vsReg = interp1( z, vs, zReg , 'linear');
vpReg = interp1( z, vp, zReg , 'linear');

% load a receiver function
DIR = 'prfns/prfns_iter_2.50';
dirlist = dir(fullfile( DIR,'TA*'));
files = dir( fullfile( DIR, dirlist(1).name, ['*',SUFFIX]) );
filename = fullfile( DIR, dirlist(1).name, files(1).name );
fprintf('filename: %s\n',filename)

try
  [t,prf,SAChdr] = sac2mat(filename);
catch ME
  error(ME.message)
end

% check the slowness units
rayp = SAChdr.user(1).data;
if rayp > 1 ,
  fprintf( 'Changing slowness from %f s/rad to %f s/km\n' , ...
	   [ rayp, rayp/6371] );
  rayp = rayp/6371;
end

%% plot the ray path
figure(1);clf;
plotRayPaths( rayp, zReg, vpReg, vsReg );

%% get the depth of the receiver function
rf2  = mapRF2depth( t, prf, rayp,...
    zReg, vpReg, vsReg, ISPRF );
 
%%
figure(2); clf;
subplot(1,2,1); plot(prf,t ); hold on;
set(gca, 'YDir', 'reverse');
set(gca, 'XTick',[]);
axis tight;
ylabel('Time (s)')

subplot(1,2,2); plot(rf2,zReg);
set(gca, 'YDir', 'reverse');
axis tight;
ylabel('Depth (km)');
set(gca, 'XTick',[]);

%% Get longitude and latitude
[epos, npos, zpos ] = sRaypath_1d( rayp, SAChdr.evsta.baz, ...
    diff(zReg(1:2)), MAXZ, zReg, vsReg );

fprintf('Back azimuth: %f degrees\n', SAChdr.evsta.baz);

% plot results
figure(3); clf;
plot3( epos+5*rf2, npos , zpos,'b'); hold on;
plot3( 0, 0, 0 , 'rv' , 'MarkerSize', 10, ...
       'MarkerFaceColor', 'r' );
view(10,30);
daspect( [ 1 1 5])
grid on;
set(gca, 'ZDir', 'reverse');
