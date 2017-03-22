function energy = pintrapfunction(pos,varargin)
% Pintrapfunction uses an analytic formula to approximate the potential
% energy of the OH magnetic pin trap. It uses a collection of weighting
% coefficients for a 3D fourier series computed with fitfourier3D.m. 

persistent Acoeffs A0 N

if ~isempty(varargin)
    Acoeffs = varargin{1};
    A0 = num2cell(varargin{2});
    N = varargin{3};     end

if isempty(Acoeffs)
    error(['A coefficients not specified. ',...
        'Provide them as a second argument to pintrapfunction ',...
        'the first time it is called. ',...
        'Provide the scales A0 and the order N ',...
        'as 3rd and 4th inputs as well.']);          end
[a,b,c,d] = A0{:};
x = abs(pos(:,1));
y = abs(pos(:,2));
z = pos(:,3);
R = d*ones(length(x),8*N^3+1);

%parpool(6);

i = 2;
for j=1:N
for k=1:N
for l=1:N
    R(:,i)=  d*sin(j*pi.*x/(2*a)).*sin(k*pi.*y/(2*b)).*sin(k*pi.*z/(2*c));
    R(:,i+1)=d*cos(j*pi.*x/(2*a)).*sin(k*pi.*y/(2*b)).*sin(k*pi.*z/(2*c));
    R(:,i+2)=d*sin(j*pi.*x/(2*a)).*cos(k*pi.*y/(2*b)).*sin(k*pi.*z/(2*c));
    R(:,i+3)=d*cos(j*pi.*x/(2*a)).*cos(k*pi.*y/(2*b)).*sin(k*pi.*z/(2*c));
    R(:,i+4)=d*sin(j*pi.*x/(2*a)).*sin(k*pi.*y/(2*b)).*cos(k*pi.*z/(2*c));
    R(:,i+5)=d*cos(j*pi.*x/(2*a)).*sin(k*pi.*y/(2*b)).*cos(k*pi.*z/(2*c));
    R(:,i+6)=d*sin(j*pi.*x/(2*a)).*cos(k*pi.*y/(2*b)).*cos(k*pi.*z/(2*c));
    R(:,i+7)=d*cos(j*pi.*x/(2*a)).*cos(k*pi.*y/(2*b)).*cos(k*pi.*z/(2*c));
    i = i + 8;
end
end
end

energy=R*Acoeffs;
        