%%
% David Reens, sometime last year.
% Commentary added 3/27/17:
%
% This function creates another function, and takes a function as an
% argument. Lovely. The argument function gives the potential energy of a
% molecule in a trap as a function of position in 3D. The generated
% function is the update function for a system of linear differential
% equations representing the motion of the particle and two orders of
% variation on that motion. In 3D, (or 6D phase space), this is an
% 18-variable system.
%
% The way that getOFLITT2_ODE_System achieves this is by symbolically
% expressing the potential function, analytically computing the derivative
% terms associated with the variations using the symbolic toolbox, and then
% exporting the ODE update function (i.e. the "f" in y' = f(y,t) ) either
% back to the caller as an anonymous function or saving in an optimized
% file by the matlabFunction converter in the symbolic toolbox.
function odefunc = getOFLITT2_ODE_System(potentialfunction,varargin)

% additional arguments are related to optimized function export
if ~isempty(varargin)
    if ~strcmp(varargin{1},'File')
        error(['Additional Arguments to getOFLITT2_ODE_System must be a '
            'file, name pair passed to matlabFunction']);
    end
end

% Create all 18 symbolic variables, group some as 6D phase space vectors.
syms x y z vx vy vz 
syms dx dy dz vdx vdy vdz 
syms d2x d2y d2z vd2x vd2y vd2z
yy = [x y z vx vy vz];
dyy = [dx dy dz vdx vdy vdz];
d2yy = [d2x d2y d2z vd2x vd2y vd2z];
assume(yy,'real')
assume(dyy,'real')
assume(d2yy,'real')

% Thinking of the target ODE function as a matrix that multiplies the input
% column vector, we can start with the submatrix that converts the
% non-variational component of phase space to its derivative. This is the
% familiar equation of motion v' = a and x' = v:
vv = potentialfunction(x,y,z);
a = -gradient(vv,[x y z]);
vs = [vx vy vz]';
xyzmat = [zeros(3) diag(vs) ; diag(a) zeros(3)];

% Here we generate the more complex parts of the variation. You can read
% about them in "Nonlinear dynamics of atoms in a crossed optical dipole
% trap" from Rosario et al in PRE 2014, Volume 90, Issue 6.
%
% Actually that paper will just refer you to Roberto Barrio's work:
% "Painting Chaos: a Gallery of Sensitivity Plots of Classical Problems"
% published in International Journal of Bifurcation and Chaos, 2006, Volume
% 16, Issue 10.
%
% I found it hard to read but sufficient nonetheless. It's also helpful to
% look at "Sensitivity tools vs. Poincaré sections", also by Roberto
% Barrio, in Chaos, Solitons & Fractals, year 2005, Volume 25, Issue 3.
f = [vx;vy;vz;a];
dfdy = arrayfun(@(xx) diff(f,xx),yy,'UniformOutput',false);
dfdy = [dfdy{:}];
dxdydzmat = dfdy;
d2fdy2 = diff(dfdy,x)*dx+diff(dfdy,y)*dy+diff(dfdy,z)*dz+...
    diff(dfdy,vx)*vdx+diff(dfdy,vy)*vdy+diff(dfdy,vz)*vdz;

% Here you can see how the big 18x18 matrix is built in blocks:
fullMat = [xyzmat   zeros(6)  zeros(6)  ; 
          zeros(6) dxdydzmat zeros(6)  ;
          zeros(6) d2fdy2    dxdydzmat ] ;

% Here is the conversion from symbolic function back to matlab function.
tempfunc = matlabFunction(fullMat * [ones(6,1) ; dyy' ; d2yy'],...
    'Vars',[yy dyy d2yy],varargin{:});

% The matlabFunction insists on generating a function that asks for 18
% separate inputs, so this anonymous function unpacks that. This can be
% done with arraylists, but it is super slow to create and then unpack an
% arraylist every time the ODE update function is called upon.
odefunc = @(t,y) tempfunc(y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),y(12),y(13),y(14),y(15),y(16),y(17),y(18));
end