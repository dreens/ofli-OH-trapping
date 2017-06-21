%% Frame version
%
% This is a merge of splineOFLI.m and run_crossed_dipole_frame.m. The goal
% is to see how well the spline version does the full panel from their
% paper. 
%
% N = number of pixels length and width of frame.
% E = energy of the frame, negative is trapped.
% X = exponent of time to run to. 3 is ideal, but slow.
function ofli2 = splineOFLIframeArb(N,E,X,trap)


% setup the grid to use for defining the spline
lim = 1;
dense = 20;
x = -lim:lim/dense:lim;
y = x;
z = x;
[xx,yy,zz] = ndgrid(x,y,z);

% function giving unitless potential energy
if strcmp(trap,'xdip')
    u = @(x,y,z) -exp(-2*x.^2-2*z.^2)/2-exp(-2*y.^2-2*z.^2)/2;
    uu = u(xx,yy,zz);
elseif strcmp(trap,'pin')
    load('pininterp.mat','uu');
end

%potential spline. It's inverted in advance: ?(-u)/?x = +f
o = 5; %order
splinepot = spapi({o+1,o+1,o+1},{x,x,z},-uu);

% Take directional derivatives to get the forces.
% a is a three valued spline; a = (ax ay az).
a = fndir(splinepot,eye(3));

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

% Setup plane of trajectories to investigate. All begin in x-y plane with
% velocity orthogonal to the plane. All have their velocity chosen
% according to their initial potential energy so that total energy is the
% same.
x = linspace(-1,1,N);
y = linspace(-1,1,N);
[xs, ys] = meshgrid(x,y);
vz = sqrt(E*2+2*splinepot(xs,ys,0));

% All ofli results in the plane will be stored in this variable
ofli2 = zeros(N,N);

% Defining this once instead of every time:
options = odeset('RelTol',2e-6,'AbsTol',2e-7);

parfor i=1:N
    for j=1:N
        % Given the specified energy, some startpoints will be out of
        % reach. We can tell if vz is imaginary:
        if ~isreal(vz(i,j))
            ofli2(i,j) = -1;
        else
            % Initialize the ODE solver
            y0 = [xs(i,j) ys(i,j) 0 0 0 vz(i,j)];
            dy0 = f(0,[y0 zeros(1,12)]');
            dy0 = dy0(1:6);
            dy0 = dy0'/sqrt(sum(dy0.^2));
            
            % Solve the ODE
            sol = ode45(f,[0 logspace(0,X,100)],[y0 dy0 zeros(1,6)],options);

            % Unpack the solutions
            y = sol.y(1:6,:);
            dy = sol.y(7:12,:);
            d2y = sol.y(13:18,:);
            
            % Note if orbit escaped the trap (or at least the spline)
            if max(max(abs(y(:,1:3))))>lim
                ofli2(i,j) = 100;
            else

                % Get flows for OFLI calculation
                flowy = y;
                for k=1:max(size(sol.y))
                    temp = f(0,sol.y(:,k));
                    flowy(1:6,k) = temp(1:6);
                end

                % Define projector
                projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;

                % Finalize OFLI terms
                fli2 = (dy+0.5*d2y);
                ofli2p = fli2 - projection(fli2,flowy);
                ofli2(i,j) = log(max(sqrt(sum(ofli2p.^2))));
            end
        end
    end
end
end