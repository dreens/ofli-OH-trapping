function energyfunc = pintrapforce(A,A0,N)
% Pintrapfunction uses an analytic formula to approximate the potential
% energy of the OH magnetic pin trap. It uses a collection of weighting
% coefficients for a 3D fourier series computed with fitfourier3D.m. 
A0 = num2cell(A0);
[a,b,c,d] = A0{:};

bs = '@(x,y,z)';
bs = [bs num2str(A(1)*d)];
i = 2;
for j=1:N
for k=1:N
for l=1:N
    if A(i)~=0
        bs = [bs '+' num2str(d*A(i)) '*sin(x*' num2str(j*pi/(2*a)) ').*sin(y*' num2str(k*pi/(2*b)) ').*sin(z*' num2str(k*pi/(2*c)) ')']; end
    if A(i+1)~=0
        bs = [bs '+' num2str(d*A(i+1)) '*cos(x*' num2str(j*pi/(2*a)) ').*sin(y*' num2str(k*pi/(2*b)) ').*sin(z*' num2str(k*pi/(2*c)) ')']; end
    if A(i+2)~=0
        bs = [bs '+' num2str(d*A(i+2)) '*sin(x*' num2str(j*pi/(2*a)) ').*cos(y*' num2str(k*pi/(2*b)) ').*sin(z*' num2str(k*pi/(2*c)) ')']; end
    if A(i+3)~=0
        bs = [bs '+' num2str(d*A(i+3)) '*cos(x*' num2str(j*pi/(2*a)) ').*cos(y*' num2str(k*pi/(2*b)) ').*sin(z*' num2str(k*pi/(2*c)) ')']; end
    if A(i+4)~=0
        bs = [bs '+' num2str(d*A(i+4)) '*sin(x*' num2str(j*pi/(2*a)) ').*sin(y*' num2str(k*pi/(2*b)) ').*cos(z*' num2str(k*pi/(2*c)) ')']; end
    if A(i+5)~=0
        bs = [bs '+' num2str(d*A(i+5)) '*cos(x*' num2str(j*pi/(2*a)) ').*sin(y*' num2str(k*pi/(2*b)) ').*cos(z*' num2str(k*pi/(2*c)) ')']; end
    if A(i+6)~=0
        bs = [bs '+' num2str(d*A(i+6)) '*sin(x*' num2str(j*pi/(2*a)) ').*cos(y*' num2str(k*pi/(2*b)) ').*cos(z*' num2str(k*pi/(2*c)) ')']; end
    if A(i+7)~=0
        bs = [bs '+' num2str(d*A(i+7)) '*cos(x*' num2str(j*pi/(2*a)) ').*cos(y*' num2str(k*pi/(2*b)) ').*cos(z*' num2str(k*pi/(2*c)) ')']; end
    i = i + 8;
end
end
end

energyfunc = eval(bs);        