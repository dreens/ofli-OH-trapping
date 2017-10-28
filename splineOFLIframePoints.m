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
function ofli2 = splineOFLIframePoints(N,E,X,trap,P,plane)

% lets collect info for a few different integration times:
T = 200;
times = linspace(0,10^X,T);



% function giving unitless potential energy
if strcmp(trap,'xdip')
    lim = 1;
    dense = 20;
    xs = -lim:lim/dense:lim;
    ys = xs ; zs = xs;
    [xx,yy,zz] = ndgrid(xs,ys,zs);
    u = @(x,y,z) -exp(-2*x.^2-2*z.^2)/2-exp(-2*y.^2-2*z.^2)/2;
    uu = u(xx,yy,zz);
elseif strcmp(trap,'pin')
    load('pininterpB.mat','uu','xs','ys','zs');
else
    error('splineOFLI:input','Trap type ''%s'' not recognized.',trap)
end

%potential spline. It's inverted in advance: ?(-u)/?x = +f
o = 4; %order
%sk = @(x) x;
sk = @(x) x(1:2:end,1:2:end,1:2:end);
splinepot = spapi({o+1,o+1,o+1},{sk(xs),sk(ys),sk(zs)},-sk(uu));

% Let's define some directories
datadir = '//data/ye/dare4983/splines/';
thisdir = sprintf('N%d_E%.1f_X%d_%s_P%d_%s',N,E,X,trap,P,plane);
mkdir(datadir,thisdir);
for mi=0:999
    mkdir([datadir thisdir],num2str(mi))
end


if ~exist([datadir 'pinC.mat'],'file')
    

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
    
    % This is useful for tuning the sign of da evaluated in the positive
    % octant:
    o123 = [1:3;1:3;1:3];

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
    [y(4:6); ...
    fnval(a,abs(y(1:3))).*sign(y(1:3)); ...
    y(10:12); ...
    fnval(da,abs(y(1:3))).*sign(y(o123)).*sign(y(o123'))*y(7:9); ...
    y(16:18); ...
     fnval(da,abs(y(1:3))).*sign(y(o123)).*sign(y(o123'))*y(13:15) + ...
    (fnval(d2ax,abs(y(1:3)))*sign(y(1)).*sign(y(o123)).*sign(y(o123'))*y(7) + ...
     fnval(d2ay,abs(y(1:3)))*sign(y(2)).*sign(y(o123)).*sign(y(o123'))*y(8) + ...
     fnval(d2az,abs(y(1:3)))*sign(y(3)).*sign(y(o123)).*sign(y(o123'))*y(9)       ) * y(7:9)   ]; 
 
    % Put the spline on the data drive maybe?
    save([datadir 'pinC.mat'],'ff','-v7.3')

else
    
    load([datadir 'pinC.mat'],'ff')
    
end

% Setup plane of trajectories to investigate. All begin in x-y plane with
% velocity orthogonal to the plane. All have their velocity chosen
% according to their initial potential energy so that total energy is the
% same.
switch(plane)
    case 'xy'
        n1 = 3; n2 = 1;
        a = linspace(0,n1,n1*N);
        va = sqrt(2*E+2*fnval(splinepot,[a;zeros(2,length(a))]));
        a = a(~imag(va));
        b = linspace(0,n2,n2*N);
        vb = sqrt(2*E+2*fnval(splinepot,[zeros(size(b));b;zeros(size(b))]));
        b = b(~imag(vb));
        [as, bs] = ndgrid(a,b);
        vc = zeros(size(as));
        cs = vc;
        vc(:) = sqrt(E*2+2*fnval(splinepot,[as(:)';bs(:)';cs(:)']));
    case 'yz'
        n1 = 1; n2 = 2;
        a = linspace(0,n1,n1*N);
        va = sqrt(2*E+2*fnval(splinepot,[zeros(size(a));a;zeros(size(a))]));
        a = a(~imag(va));
        b = linspace(0,n2,n2*N);
        vb = sqrt(2*E+2*fnval(splinepot,[zeros(2,length(b));b]));
        b = b(~imag(vb));
        [as, bs] = ndgrid(a,b);
        vc = zeros(size(as));
        cs = vc;
        vc(:) = sqrt(E*2+2*fnval(splinepot,[cs(:)';as(:)';bs(:)']));
    case 'xz'
        n1 = 3; n2 = 2;
        a = linspace(0,n1,n1*N);
        va = sqrt(2*E+2*fnval(splinepot,[a;zeros(2,length(a))]));
        a = a(~imag(va));
        b = linspace(0,n2,n2*N);
        vb = sqrt(2*E+2*fnval(splinepot,[zeros(2,length(b));b]));
        b = b(~imag(vb));
        [as, bs] = ndgrid(a,b);
        vc = zeros(size(as));
        cs = vc;
        vc(:) = sqrt(E*2+2*fnval(splinepot,[as(:)';cs(:)';bs(:)']));
    otherwise
        error('splineOFLI:input',...
       'Plane input must be among the following:\n xy\n yz\n xz\n.')
end
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
N1 = length(a); N2 = length(b); NT = N1*N2;
ofli2 = zeros(NT,T);
save([datadir thisdir '/dims_' num2str(N1) 'x' num2str(N2) '.mat'],'N1','N2');

% Randomize parfor ordering
rindex = randperm(NT);

vcp = vc; vcp(:) = vc(rindex);
asp = as; asp(:) = as(rindex);
bsp = bs; bsp(:) = bs(rindex);

parfor iii=1:NT
    
    i = rindex(iii);
    
    % Declare this slice of the 3D variable:
    oflirow = zeros(1,T);

    
    % Look to see if this step in the parfar loop has been completed by a
    % previously attempted job submission. If so, get the results and reuse
    % them directly.
    filestash = [datadir thisdir '/' num2str(mod(i,1000)) '/i=' num2str(i) '.mat'];
    if exist(filestash,'file')
        d = load(filestash);
        for k=1:T
            oflirow(k) = d.stash(k);
        end
    else

        % Now get the spline back out now that we're on a compute node:
        % f = load([datadir 'pinB.mat']);

        % Now define options for the ODE solver. Would have done this outside
        % the parfor, but it needs the spline that we need to separately load
        % on the compute nodes to reduce information transfer.
        escape = @(t,y) escbase(t,y,ff);
        options = odeset('RelTol',10^-P,'AbsTol',10^-P,'Events',escape);
        
        % Given the specified energy, some startpoints will be out of
        % reach. We can tell if vc is imaginary:
        oflirow = zeros(1,T);
        if ~isreal(vcp(iii))
            oflirow = oflirow - 1;
        else
            % Initialize the ODE solver
            switch(plane)
                case 'xy'
                    y0 = [asp(iii) bsp(iii) 0 0 0 vcp(iii)];
                case 'yz'
                    y0 = [0 asp(iii) bsp(iii) vcp(iii) 0 0];
                case 'xz'
                    y0 = [asp(iii) 0 bsp(iii) 0 vcp(iii) 0];
            end
            dy0 = ff(0,[y0 zeros(1,12)]');
            dy0 = dy0(1:6);
            dy0 = dy0'/sqrt(sum(dy0.^2));

            % Solve the ODE
            sol = ode113(ff,times,[y0 dy0 zeros(1,6)],options);

            % Check for errors
            if ~isempty(sol.ie)
                oflirow = oflirow + sol.ie;
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
                temp = ff(0,yall(:,l));
                flowy(1:6,l) = temp(1:6);
            end

            % Define projector
            projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;

            % Finalize OFLI terms
            fli2 = (dy+0.5*d2y);
            ofli2p = fli2 - projection(fli2,flowy);
            oflirow(1:length(ntimes)) = log10(sqrt(sum(ofli2p.^2)));
        end
        
        % Save the results of this step of the parfor loop directly so that
        % if the job is killed or runs out of time or jiliac crashes, we
        % don't start from scratch.
        saveparfor(filestash,oflirow)
        
    end % end if exist filestash
    ofli2(iii,:) = oflirow;
end % end parfor i=1:N

ofli2(rindex,:) = ofli2(:,:);

ofli2 = reshape(ofli2,[N1 N2 T]);

end % end function splineOFLIframe

function saveparfor(file,stash)
    save(file,'stash','-v7.3')
end