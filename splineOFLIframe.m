%% Frame version
%
% This is a merge of splineOFLI.m and run_crossed_dipole_frame.m. The goal
% is to see how well the spline version does the full panel from their
% paper. 
%
% N = number of pixels length and width of frame.
% E = energy of the frame, negative is trapped.
% X = exponent of time to run to. 3 is ideal, but slow.
% trap = "xdip" for crossed dipole or "pin" for pintrap.
% Assumes X=0, Y=0 mirror symmetry.
function ofli2 = splineOFLIframe(N,E,X,trap,P)

% lets collect info for a few different integration times:
T = 20;
times = logspace(X-2,X,20);


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

if ~exist('//data/ye/dare4983/splines/pin.mat','file')
    

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
    ff = @(t,y) ...
    [y(4:6); fnval(a,y(1:3)); y(10:12); fnval(da,y(1:3))*y(7:9); y(16:18); ...
     fnval(da,y(1:3))*y(13:15) + ...
    (fnval(d2ax,y(1:3))*y(7) + ...
     fnval(d2ay,y(1:3))*y(8) + ...
     fnval(d2az,y(1:3))*y(9)       ) * y(7:9)   ]; 
 
    % Put the spline on the data drive maybe?
    save('//data/ye/dare4983/splines/pin.mat','ff','-v7.3')

end

% Setup plane of trajectories to investigate. All begin in x-y plane with
% velocity orthogonal to the plane. All have their velocity chosen
% according to their initial potential energy so that total energy is the
% same.
x = linspace(0,1,N);
y = linspace(0,1,N);
[xs, ys] = meshgrid(x,y);
vz = zeros(size(xs));
vz(:) = sqrt(E*2+2*fnval(splinepot,[xs(:)';ys(:)';zeros(size(xs(:)'))]));


% This function checks for escaping the trap, and is registered as an
% "event" indicating termination by the ode solver.
% function [vals, terms, dirs] = escbase(~,y,fff) 
%     y2 = y(7:12);
%     y3 = y(13:18);
%     fy = fff(0,y);
%     fy = fy(1:6);
%     fl = (y2+0.5*y3);
%     pj = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;
%     ofl = fl - pj(fl,fy);
%     ofl = log10(sqrt(sum(ofl.^2)));
%     
%     val1 = max(abs(y(1:6)))-10;
%     val2 = ofl-20;
%     vals = [val1 val2];
%     terms = [1 1];
%     dirs = [0 0];
% end

% All ofli results in the plane will be stored in this variable
ofli2 = zeros(N,N,T);

parfor i=1:N
    
    % Now get the spline back out now that we're on a compute node:
    f = load('//data/ye/dare4983/splines/pin.mat');
    
    % Now define options for the ODE solver. Would have done this outside
    % the parfor, but it needs the spline that we need to separately load
    % on the compute nodes to reduce information transfer.
    escape = @(t,y) escbase(t,y,f.ff);
    options = odeset('RelTol',10^-P,'AbsTol',10^-P,'Events',escape);

    
    fprintf('%d\n',i)
    for j=1:N
        % Given the specified energy, some startpoints will be out of
        % reach. We can tell if vz is imaginary:
        ofli2row = zeros(1,T);
        if ~isreal(vz(i,j))
            ofli2row = ofli2row - 1;
        else
            % Initialize the ODE solver
            y0 = [xs(i,j) ys(i,j) 0 0 0 vz(i,j)];
            dy0 = f.ff(0,[y0 zeros(1,12)]');
            dy0 = dy0(1:6);
            dy0 = dy0'/sqrt(sum(dy0.^2));

            % Solve the ODE
            sol = ode45(f.ff,times,[y0 dy0 zeros(1,6)],options);

            % Check for erors
            if ~isempty(sol.ie)
                ofli2row = ofli2row + sol.ie;
            end
            
            % Unpack the solutions
            ntimes = times(times<=max(sol.x));
            yall = deval(sol,ntimes);
            y   =  yall(1:6,:);
            dy  =  yall(7:12,:);
            d2y =  yall(13:18,:);

            % Get flows for OFLI calculation
            flowy = y;
            for l=1:length(ntimes)
                temp = f.ff(0,yall(:,l));
                flowy(1:6,l) = temp(1:6);
            end
            
            % Define projector
            projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;

            % Finalize OFLI terms
            fli2 = (dy+0.5*d2y);
            ofli2p = fli2 - projection(fli2,flowy);
            ofli2row(1:length(ntimes)) = log10(sqrt(sum(ofli2p.^2)));
        end
        for k=1:T
            ofli2(i,j,k) = ofli2row(k);
        end
    end
end
end