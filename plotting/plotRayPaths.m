function plotRayPaths( PP, REGZ, REGVP , REGVS )

% plotRayPaths( PP, REGZ, REGVP , REGVS )
%
% For a ray parameter and a velocity model compute the P and S ray path
%
% IN:
% PP = horizontal slowness (s/km)
% REGZ = regularly spaced depths for midpoints in layer const
%        velocity model (km)
% REGVP = corresponding vp of reg spaced vel model (km/s)
% REGVS = corresponding vs of reg spaced vel model (km/s)
%

% get associated vertical slowness
qb = sqrt(1./REGVS.^2 - PP^2);
qa = sqrt(1./REGVP.^2 - PP^2);

% compute horizontal travel for each 
dz = REGZ(2) - REGZ(1); % distance of layers
dx = PP * dz ./ qb;

% convert to horizontal coord
xs = [ 0 ; cumsum(dx) ];
dx = PP * dz ./ qa;
xp = [ 0 ; cumsum(dx) ];

% depth of interfaces
zintf = [ 0 ; REGZ + 0.5*dz ];

% Make Figure 
clf;
subplot(1,2,1)

% plot s ray path
plot( xs, zintf ,'-k'); hold on;

% compare to straight line
plot([xs(1), xs(end)],[REGZ(1), REGZ(end)],'--','Color',[.8 .8 .8])

% plot p ray path
plot( xp, zintf ,'-r'); hold on;

% compare to straight line
plot([xp(1), xp(end)], [REGZ(1), REGZ(end)],'--','Color',[.8 0 0])
xlabel('x [km]');
ylabel('z [km]')
set(gca, 'YDir', 'reverse');
set(gca,'XAxisLocation','top')
grid on;
axis equal;

% plot the velocity model
subplot(1,2,2)
plot( REGVS, REGZ ,'k'); hold on;
plot( REGVP, REGZ ,'r'); hold on;
xlabel('V [km/s]');
axis( [ 0, max(REGVP), 0, REGZ(end)] )
set(gca, 'YDir', 'reverse');
set(gca, 'YTickLabel',[]);
set(gca,'XAxisLocation','top');
grid on;
