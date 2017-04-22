%% Repeat Crossed Dipole Calculations but with a cubic spline!
%
% This function combines the capabilities of several previous functions.
% getOFLIYY2_ODE_System could produce an arbitrary propagator for checking
% higher order trajectory perturbations given an analytic geometry, but
% this function will take points in and produces a cubic spline.

%% First we start with crossed dipole potential for comparison
%
% U = -exp(-2x^2-2z^2)/2-exp(-2y^2-2z^2)/2
u = @(x,y,z) -exp(-2*x.^2-2*z.^2)/2-exp(-2*y.^2-2*z.^2)/2;

lim = 1;
dense = 20;
x = -lim:lim/dense:lim;
y = x;
z = x;
[xx,yy,zz] = ndgrid(x,y,z);
uu = u(xx,yy,zz);

% see what it looks like:
% figure; fnplt(uu);

%pcs for potential cubic spline. It's inverted in advance: ?(-u)/?x = +f
pcs = csapi({x,y,z},-uu);

% Take directional derivatives to get the forces.
% a is a three valued spline; a = (ax ay az).
a = fndir(pcs,eye(3));

% % Setup an ODE and solve it.
% f = @(t,y) [y(4); y(5); y(6); fnval(a,y(1:3))] ;
% t = linspace(0,10,100);
% y0 = [0 0 0 0 0 1.2];
% options = odeset('RelTol',2e-9,'AbsTol',2e-10);
% ode45(f,t,y0,options);

% These will be used for the "df/dy" term:
da = fndir(a,eye(3));

% make da output 3x3 matrices instead of 9x1 vector. Here is it's output:
%                     ( da_x/dx  da_x/dy  da_x/dz )
% fnval(da,[x;y;z]) = ( da_y/dx  da_y/dy  da_y/dz )
%                     ( da_z/dx  da_z/dy  da_z/dz )
da = fnchg(da,'dim',[3 3]); 

% here I take second derivatives. In principle this is a single tensor, 
% but I've split it up for simplicity.
d2ax = fnchg(fndir(da,[1;0;0]),'dim',[3 3]);
d2ay = fnchg(fndir(da,[0;1;0]),'dim',[3 3]);
d2az = fnchg(fndir(da,[0;0;1]),'dim',[3 3]);

% here is the actual ODE stepping function. The variables are ordered in
% one giant length 18 column vector: 
% [x y z vx vy vz dx dy dz dvx dvy dvz d2x d2y d2z d2vx d2vy d2vz]'
% This function gives what the derivative should be for each of these in
% the coupled differential equation:
% y' = f(y), dy' = (?f/?y)*dy, d2y = (?f/?y)*d2y+(?^2f/?^2y)ij*dyi*dyj
f = @(t,y) ...
[y(4:6); fnval(a,y(1:3)); y(10:12); fnval(da,y(1:3))*y(7:9); y(16:18); ...
 fnval(da,y(1:3))*y(13:15) + ...
(fnval(d2ax,y(1:3))*y(7) + ...
 fnval(d2ay,y(1:3))*y(8) + ...
 fnval(d2az,y(1:3))*y(9)       ) * y(7:9)   ]; 

%% 
% Now I try to repeat figure 3 of "nonlinear dynamics of atoms in a crossed
% optical dipole trap" González-Férez et al, PRE 2014.
% Get vz with total energy -0.6:
xi = 430e-3; yi = 430e-3; zi = 0;
vz = sqrt(-.6*2-2*u(xi,yi,zi));

% Initial position.
y0 = [xi yi zi 0 0 vz];

% Initial dy. Should be along the "flow" and a unit vector:
dy0 = f(0,[y0 zeros(1,12)]');
dy0 = dy0(1:6);
dy0 = dy0'/sqrt(sum(dy0.^2));

% Initial d2y=0.
d2y0 = zeros(1,6);

% Solve the ODE and keep track of how long it takes.
tic
options = odeset('RelTol',1e-6,'AbsTol',1e-7);%,'OutputFcn','odeplot');
sol = ode45(f,[0 1000],[y0 dy0 d2y0],options);
toc

%
% Unpack the solutions
y = sol.y(1:6,:);
dy = sol.y(7:12,:);
d2y = sol.y(13:18,:);

% Get flows for OFLI calculation
flowy = y;
flowd2y = y;
for i=1:max(size(sol.y))
    temp = f(0,sol.y(:,i));
    flowy(1:6,i) = temp(1:6);
    flowd2y(1:6,i) = temp(13:18);
end

% Define projector
projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;

% Finalize OFLI terms
ofli1 = dy - projection(dy,flowy);
fli2 = (dy+.5*d2y);
ofli2 = fli2 - projection(fli2,flowy);
ofli2n = sqrt(sum(ofli2.^2));

%figure(23); hold on;
%plot(sol.x,log10(ofli2n))

figure(111111)
hold on
plot(sol.x,ofli2n)
%ylim([-10 25])
%xlim([0 1000])


%
% Now I try to repeat figure 3 of "nonlinear dynamics of atoms in a crossed
% optical dipole trap" González-Férez et al, PRE 2014.
% Get vz with total energy -0.6:
xi = 400e-3; yi = 0e-3; zi = 0;
vz = sqrt(-.6*2-2*u(xi,yi,zi));

% Initial position.
y0 = [xi yi zi 0 0 vz];

% Initial dy. Should be along the "flow" and a unit vector:
dy0 = f(0,[y0 zeros(1,12)]');
dy0 = dy0(1:6);
dy0 = dy0'/sqrt(sum(dy0.^2));

% Initial d2y=0.
d2y0 = zeros(1,6);

% Solve the ODE and keep track of how long it takes.
tic
options = odeset('RelTol',1e-6,'AbsTol',1e-7);%,'OutputFcn','odeplot');
sol = ode45(f,[0 1000],[y0 dy0 d2y0],options);
toc

%
% Unpack the solutions
y = sol.y(1:6,:);
dy = sol.y(7:12,:);
d2y = sol.y(13:18,:);

% Get flows for OFLI calculation
flowy = y;
flowd2y = y;
for i=1:max(size(sol.y))
    temp = f(0,sol.y(:,i));
    flowy(1:6,i) = temp(1:6);
    flowd2y(1:6,i) = temp(13:18);
end

% Define projector
projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;

% Finalize OFLI terms
ofli1 = dy - projection(dy,flowy);
fli2 = (dy+.5*d2y);
ofli2 = fli2 - projection(fli2,flowy);
ofli2n = sqrt(sum(ofli2.^2));

%figure(23); hold on;
%plot(sol.x,log10(ofli2n))

figure(111111)
hold on
plot(sol.x,ofli2n)
%ylim([-10 25])
%xlim([0 1000])



%
% Now I try to repeat figure 3 of "nonlinear dynamics of atoms in a crossed
% optical dipole trap" González-Férez et al, PRE 2014.
% Get vz with total energy -0.6:
xi = 0e-3; yi = 0e-3; zi = 0;
vz = sqrt(-.6*2-2*u(xi,yi,zi));

% Initial position.
y0 = [xi yi zi 0 0 vz];

% Initial dy. Should be along the "flow" and a unit vector:
dy0 = f(0,[y0 zeros(1,12)]');
dy0 = dy0(1:6);
dy0 = dy0'/sqrt(sum(dy0.^2));

% Initial d2y=0.
d2y0 = zeros(1,6);

% Solve the ODE and keep track of how long it takes.
tic
options = odeset('RelTol',1e-6,'AbsTol',1e-7);%,'OutputFcn','odeplot');
sol = ode45(f,[0 1000],[y0 dy0 d2y0],options);
toc

%
% Unpack the solutions
y = sol.y(1:6,:);
dy = sol.y(7:12,:);
d2y = sol.y(13:18,:);

% Get flows for OFLI calculation
flowy = y;
flowd2y = y;
for i=1:max(size(sol.y))
    temp = f(0,sol.y(:,i));
    flowy(1:6,i) = temp(1:6);
    flowd2y(1:6,i) = temp(13:18);
end

% Define projector
projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;

% Finalize OFLI terms
ofli1 = dy - projection(dy,flowy);
fli2 = (dy+.5*d2y);
ofli2 = fli2 - projection(fli2,flowy);
ofli2n = sqrt(sum(ofli2.^2));

%figure(23); hold on;
%plot(sol.x,log10(ofli2n))

figure(111111)
hold on
plot(sol.x,ofli2n)
%ylim([-10 25])
%xlim([0 1000])




