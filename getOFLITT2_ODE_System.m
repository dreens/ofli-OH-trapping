function odefunc = getOFLITT2_ODE_System(potentialfunction,varargin)

if ~isempty(varargin)
    if ~strcmp(varargin{1},'File')
        error(['Additional Arguments to getOFLITT2_ODE_System must be a '
            'file, name pair passed to matlabFunction']);
    end
end

syms x y z vx vy vz 
syms dx dy dz vdx vdy vdz 
syms d2x d2y d2z vd2x vd2y vd2z
yy = [x y z vx vy vz];
dyy = [dx dy dz vdx vdy vdz];
d2yy = [d2x d2y d2z vd2x vd2y vd2z];
assume(yy,'real')
assume(dyy,'real')
assume(d2yy,'real')

%mass = 17*1.67e-27;
vv = potentialfunction(x,y,z);
a = -gradient(vv,[x y z]);
vs = [vx vy vz]';
xyzmat = [zeros(3) diag(vs) ; diag(a) zeros(3)];
f = [vx;vy;vz;a];

dfdy = arrayfun(@(xx) diff(f,xx),yy,'UniformOutput',false);
dfdy = [dfdy{:}];
dxdydzmat = dfdy;
d2fdy2 = diff(dfdy,x)*dx+diff(dfdy,y)*dy+diff(dfdy,z)*dz+...
    diff(dfdy,vx)*vdx+diff(dfdy,vy)*vdy+diff(dfdy,vz)*vdz;
fullMat = [xyzmat   zeros(6)  zeros(6)  ; 
          zeros(6) dxdydzmat zeros(6)  ;
          zeros(6) d2fdy2    dxdydzmat ] ;

tempfunc = matlabFunction(fullMat * [ones(6,1) ; dyy' ; d2yy'],...
    'Vars',[yy dyy d2yy],varargin{:});

odefunc = @(t,y) tempfunc(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),y(12),y(13),y(14),y(15),y(16),y(17),y(18));
end