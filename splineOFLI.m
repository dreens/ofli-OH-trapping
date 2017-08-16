%% Repeat Crossed Dipole Calculations but with a cubic spline!
%
% This function combines the capabilities of several previous functions.
% getOFLIYY2_ODE_System could produce an arbitrary propagator for checking
% higher order trajectory perturbations given an analytic geometry, but
% this function will take points in and produces a cubic spline.
%
% First repeated different types of behavior- periodic, quasiperiodic, and
% chaotic- on April 22, 2017. See spline_success.fig.
%
% Rebranded as a function that can evaluate OFLI at a single point. I
% basically took the modifications done to get a frame running on the
% cluster and brought them back to some extent into this earlier edition.

function out = splineOFLI(E,X,trap)

%% First we start with crossed dipole potential for comparison
%

lim = 1;
dense = 20;
x = -lim:lim/dense:lim;
y = x;
z = x;
[xx,yy,zz] = ndgrid(x,y,z);


if strcmp(trap,'xdip')
    u = @(x,y,z) -exp(-2*x.^2-2*z.^2)/2-exp(-2*y.^2-2*z.^2)/2;
    uu = u(xx,yy,zz);
elseif strcmp(trap,'pin')
    load('pininterp.mat','uu');
end


% see what it looks like:
% figure; fnplt(uu);

%potential spline. It's inverted in advance: ?(-u)/?x = +f
o = 5; %order
splinepot = spapi({o+1,o+1,o+1},{x,x,z},-uu);

% Take directional derivatives to get the forces.
% a is a three valued spline; a = (ax ay az).
a = fndir(splinepot,eye(3));

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
xi = 100e-3; yi = 100e-3; zi = 0;
vz = sqrt(E*2+2*fnval(splinepot,[xi; yi; 0]));

% Initial position.
y0 = [xi yi zi 0 0 vz];

% Initial dy. Should be along the "flow" and a unit vector:
dy0 = f(0,[y0 zeros(1,12)]');
dy0 = dy0(1:6);
dy0 = dy0'/sqrt(sum(dy0.^2));

% Initial d2y=0.
d2y0 = zeros(1,6);

% Output function trajectories only
function stop = trimplot(t,y,flag,varargin)
    if size(y,1)>=6
        stop = odeplot(t,y(13:18,:),flag,varargin);
    else
        stop = false;
    end
end

% Try plotting OFLI live
function stop = oflilive(t,y,flag,varargin)
    if size(y,1)>=6
        y1 = y(1:6,:);
        y2 = y(7:12,:);
        y3 = y(13:18,:);
        fy = zeros(size(y));
        for ii=1:size(y,2)
            fy(:,ii) = f(0,y(:,ii));
        end
        fy = fy(1:6,:);
        fl = (y2+0.5*y3);
        pj = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;
        ofl = fl - pj(fl,fy);
        ofl = sqrt(sum(ofl.^2));

        stop = odeplot(t,ofl,flag,varargin);
    else
        stop = false;
    end
end

% "Event" Function to exit if the molecule flies out.
function [vals, terms, dirs] = escape(~,y) 
    y2 = y(7:12);
    y3 = y(13:18);
    fy = f(0,y);
    fy = fy(1:6);
    fl = (y2+0.5*y3);
    pj = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;
    ofl = fl - pj(fl,fy);
    ofl = sqrt(sum(ofl.^2));
    
    val1 = max(abs(y(1:6)))-10;
    val2 = ofl-10^1;
    vals = [val1 val2];
    terms = [1 1];
    dirs = [0 0];
end

% Solve the ODE and keep track of how long it takes.
options = odeset('RelTol',2e-8,'AbsTol',2e-8,'OutputFcn',@oflilive,'Events',@escape);
sol = ode45(f,[0 10 50 100],[y0 dy0 d2y0],options);

% Check for escape
if sol.x(end)<10^X
    out = -2;
else

    % Unpack the solutions
    y = sol.y(1:6,:);
    dy = sol.y(7:12,:);
    d2y = sol.y(13:18,:);

    % Get flows for OFLI calculation
    flowy = y;
    for i=1:max(size(sol.y))
        temp = f(0,sol.y(:,i));
        flowy(1:6,i) = temp(1:6);
    end

    % Define projector
    projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;

    % Finalize OFLI terms
    fli2 = (dy+0.5*d2y);
    ofli2p = fli2 - projection(fli2,flowy);
    ofli2n = sqrt(sum(ofli2p.^2));


    out = sol;
end

end

% check energy
%figure;hold on
%plot(0.7-fnval(splinepot,sol.y(1:3,:))+.5*sum(sol.y(4:6,:).^2))
% With rel 1e-6, abs 1e-7, 1.6e-8 energy lost in time = 100;

%figure(111111)
%hold on
%plot(sol.x,ofli2n)
%grid on
%set(gca,'YScale','log')
%set(gca,'XScale','log')
%xlim([1 100])
%ylim([0.1 100])