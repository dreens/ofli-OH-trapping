%% Prepare Pin Trap Spline
%
% Prepare a matrix for splining the pin trap potential.
f = load('magpinbroadmatrix.mat');
x = -2:0.1:2;
y = x;
z = y;
[xx, yy, zz] = ndgrid(x,y,z);

uu = interpn(f.xx,f.yy,f.zz,f.vv,xx,yy,zz);
uu = uu - max(uu(:));
uu = -uu./(min(uu(:)));

save('pininterp.mat','uu')