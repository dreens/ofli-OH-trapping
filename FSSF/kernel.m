function output=kernel(data,n,series,bcx_0,bcx_end,bcy_0,bcy_end)

% fitfourier(data,n,bcx_0,bcx_end,bcy_0,bcy_end)
% Function for fitting a 2-dimensional Fourier series
% to an unstructured 3-dimensional point map.
% Fourier series formulation in the fashion which appears
% in the linear elastic plate bending theory.
% Fitting is performed by means of LSQ
% data :: matrix (of at least n x 3 dimensions) containing the original data
% n :: order of the Fourier series
% bcx_0 :: boundary condition for the x-start --> 1 for supported, 2 for free
% bcx_end :: boundary condition for the x-end --> 1 for supported, 2 for free   
% bcy_0 :: boundary condition for the y-start --> 1 for supported, 2 for free
% bcy_end :: boundary condition for the y-end --> 1 for supported, 2 for free
% function returns the vector of Aij coefficients and plots the resulting surface with deviation data
% additional columns of the input matrix are to be ignored
%(c) Miltiades Salvanos 2003

datasize=size(data);
if (datasize(2))<3
    error('insufficient data, breaking...');
else
    x=data(:,1); y=data(:,2); z=data(:,3);
end
switch series;
case 'SIN*SIN'
    a=abs(max(x)-min(x)); b=abs(max(y)-min(y));
    j=1; k=1; R=[];
    for j=1:2:n
        for k=1:2:n
            Rjk=max(abs(z))*sin(j*pi.*x/(bcx_end*a)+(bcx_0-1)*pi/(bcx_end*2)).*sin(k*pi.*y/(bcy_end*b)+(bcy_0-1)*pi/(bcy_end*2));  %Fourier Series
            R=[R Rjk];
        end
    end
    A=R\z;
case '(COS^2)*(COS^2)'
    a=abs(max(x)-min(x)); b=abs(max(y)-min(y));
    j=1; k=1; R=[];
    for j=1:2:n
        for k=1:2:n
            Rjk=max(abs(z))*((cos(j*pi.*x/(bcx_end*a)+(bcx_0-1)*pi/(bcx_end*2))).^2).*((sin(k*pi.*y/(bcy_end*b)+(bcy_0-1)*pi/(bcy_end*2))).^2);  %Fourier Series
            R=[R Rjk];
        end
    end
    A=R\z;
otherwise
    error('unknown series formulation...');
    return;
end

% ----------
% evaluation of surface elevation and calculation of maximum deviation
% ----------
zfitted=R*A;
Maximum_Deviation=max(abs(zfitted-z));
output.coeff=A;
output.Max_Dev=Maximum_Deviation;
output.zfitted=zfitted;