%% Prepare Pin Trap Spline
%
% Prepare a matrix for splining the pin trap potential.
f = importdata('Broad322DynamicsEnergyt.dat',' ',9);
x = f.data(:,1);
y = f.data(:,2);
z = f.data(:,3);
u = f.data(:,4);

xs = sort(uniquetol(x,1e-6,'DataScale',1));
ys = sort(uniquetol(y,1e-6,'DataScale',1));
zs = sort(uniquetol(z,1e-6,'DataScale',1));

% get the datapoint spacing, used for taking derivatives later.
xsp = mode(diff(xs));
ysp = mode(diff(ys));
zsp = mode(diff(zs));

% the size of the 3D data matrices to be filled
fullsize = [length(xs) length(ys) length(zs)];

% for each x,y,z value in the COMSOL loaded x,y,z columns, check which
% linear index this corresponds to in a 3D matrix of size fullsize.
x2i = @(x) int16(1+x/xsp);
y2i = @(y) int16(1+y/ysp);
z2i = @(z) int16(1+z/zsp);
locs = sub2ind(fullsize,x2i(x),y2i(y),z2i(z));

% Now fill the matrices:
xx = zeros(fullsize);
yy = zeros(fullsize);
zz = zeros(fullsize);
uu = zeros(fullsize);
xx(locs) = x;
yy(locs) = y;
zz(locs) = z;
uu(locs) = u;



%x = -3:0.1:3;
%y = -1:0.1:1;
%z = -2.5:0.1:2.5;
%[xxx, yyy, zzz] = ndgrid(xs,ys,zs);

%uu = interpn(xx,yy,zz,vv,xxx,yyy,zzz);
uu = uu - min(uu(:));
uu = uu/10-1;

save('pininterpB.mat','uu','xs','ys','zs')