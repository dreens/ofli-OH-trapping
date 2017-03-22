%% Crossed Dipole Chaos
%
% Here I actually try to reproduce Rosario et. al. 's results.
%
% Successfully saw distinction between chaotic, quasiperiodic, periodic
% after correcting several typos on July 31, 2016. -D.R.
%
function ofli2 = run_crossed_dipole_frame(N)
potentialenergy = @(x,y,z) -.5*(exp(-2*(z.^2+y.^2))+exp(-2*z.^2-2*x.^2));
%potentialenergy = @(x,y,z) x^2+y^2+z^2;

%odefunc = getOFLITT2_ODE_System(potentialenergy);%,'File','crossedDipole');
%N = 3;

x = linspace(-1,1,N);
y = linspace(-1,1,N);
[xs, ys] = meshgrid(x,y);
vz = sqrt(-.55*2-2*potentialenergy(xs,ys,0));
d2y0 = zeros(1,6);
options = odeset('RelTol',2e-6,'AbsTol',2e-7);
projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;
ofli2 = zeros(N,N);

parfor i=1:N
    for j=1:N
        if ~isreal(vz(i,j))
            ofli2(i,j) = -1;
        else
            y0 = [xs(i,j) ys(i,j) 0 0 0 vz(i,j)];
            dy0 = crossedDipole(0,[y0 zeros(1,12)]);
            dy0 = dy0(1:6);
            dy0 = dy0'/sqrt(sum(dy0.^2));
            [~, yy] = ode45(@crossedDipole,[0 logspace(0,3,100)],[y0 dy0 d2y0],options);

            y = yy(:,1:6)';
            dy = yy(:,7:12)';
            d2y = yy(:,13:18)';
            
            if max(max(abs(y(:,1:3))))>5
                ofli2(i,j) = 100;
            else

                flowy = y;
                for k=1:max(size(y))
                    temp = crossedDipole(0,yy(k,:));
                    flowy(1:6,k) = temp(1:6);
                end

                fli2 = (dy+.5*d2y);
                ofli2p = fli2 - projection(fli2,flowy);
                ofli2(i,j) = log(max(sqrt(sum(ofli2p.^2))));
            end
        end
    end
end
end