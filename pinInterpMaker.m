%% Prepare Pin Trap Spline
%
% Prepare a matrix for splining the pin trap potential.
f = load('magpinbroadmatrix.mat');
x = -3:0.1:3;
y = -1:0.1:1;
z = -2.5:0.1:2.5;
[xx, yy, zz] = ndgrid(x,y,z);

uu = interpn(f.xx,f.yy,f.zz,f.vv,xx,yy,zz);
uu = uu - max(uu(:));
uu = -uu./(min(uu(:)));

save('pininterpB.mat','uu')