%% Crossed Dipole Chaos
%
% Here I actually try to reproduce Rosario et. al. 's results.
%
% Successfully saw distinction between chaotic, quasiperiodic, periodic
% after correcting several typos on July 31, 2016. -D.R.
%

potentialenergy = @(x,y,z) -.5*(exp(-2*(z^2+y^2))+exp(-2*z*z-2*x*x));
%potentialenergy = @(x,y,z) x^2+y^2+z^2;

odefunc = getOFLITT2_ODE_System(potentialenergy);%,'File','crossedDipole');


vz = sqrt(-.6*2-2*potentialenergy(500e-3,0e-3,0));
y0 = [500e-3 0e-3 0 0 0 vz];
dy0 = odefunc(0,[y0 zeros(1,12)]);
dy0 = dy0(1:6);
dy0 = dy0'/sqrt(sum(dy0.^2));
d2y0 = zeros(1,6);
%dy0 = [1 0 0 0 0 0];
tic
options = odeset('RelTol',2e-6,'AbsTol',2e-7);
options = odeset('RelTol',1e-7,'AbsTol',1e-8);
sol = ode45(@crossedDipole,[0 1000],[y0 dy0 d2y0],options);
toc
if false
figure;
plot(sol.y(13:18,:)')
end
y = sol.y(1:6,:);
dy = sol.y(7:12,:);
d2y = sol.y(13:18,:);
% y = sol.y(1:3,:);
% vy = sol.y(4:6,:);
% dy = sol.y(7:9,:);
% vdy = sol.y(10:12,:);
% d2y = sol.y(13:15,:);
% vd2y = sol.y(16:18,:);

flowy = y;
flowd2y = y;
for i=1:max(size(sol.y))
    temp = odefunc(0,sol.y(:,i)');
    flowy(1:6,i) = temp(1:6);
    flowd2y(1:6,i) = temp(13:18);
end
%flowy = vy;

projection = @(a,b) repmat(sum(a.*b)./sum(b.^2),size(a,1),1).*b;

ofli1 = dy - projection(dy,flowy);
fli2 = (dy+.5*d2y);
ofli2 = fli2 - projection(fli2,flowy);
ofli2n = sqrt(sum(ofli2.^2));

%figure(23); hold on;
%plot(sol.x,log10(ofli2n))

figure(111111)
hold on
plot(sol.x,ofli2n)
ylim([-10 25])
xlim([0 1000])