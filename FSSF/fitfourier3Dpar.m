function [A,A0,dev] = fitfourier3Dpar(data,n)

% fitfourier(data,n,bcx_0,bcx_end,bcy_0,bcy_end)
% Function for fitting a 2-dimensional Fourier series
% to an unstructured 3-dimensional point map.
% Fourier series formulation in the fashion which appears
% in the linear elastic plate bending theory.
% Fitting is performed by means of LSQ
% data :: matrix (of at least n x 3 dimensions) containing the original data
% n :: order of the Fourier series
% function returns the vector of Aij coefficients and plots the resulting surface with deviation data
% additional columns of the input matrix are to be ignored
%(c) Miltiades Salvanos 2003
x = data(:,1);
y = data(:,2);
z = data(:,3);
w = data(:,4);

datasize=size(data);
a=abs(max(x)-min(x)); 
b=abs(max(y)-min(y));
c=abs(max(z)-min(z));
d = max(abs(w));
A0 = [a b c d];

R = zeros(length(x),8*n^3+1,'distributed');
R(:,1) = d*ones(size(x));

i = 2;
for j=1:n
for k=1:n
for l=1:n
    for m=1:length(x)
        R(m,i)=  d*sin(j*pi.*x(m)/(2*a)).*sin(k*pi.*y(m)/(2*b)).*sin(k*pi.*z(m)/(2*c));
        R(m,i+1)=d*cos(j*pi.*x(m)/(2*a)).*sin(k*pi.*y(m)/(2*b)).*sin(k*pi.*z(m)/(2*c));
        R(m,i+2)=d*sin(j*pi.*x(m)/(2*a)).*cos(k*pi.*y(m)/(2*b)).*sin(k*pi.*z(m)/(2*c));
        R(m,i+3)=d*cos(j*pi.*x(m)/(2*a)).*cos(k*pi.*y(m)/(2*b)).*sin(k*pi.*z(m)/(2*c));
        R(m,i+4)=d*sin(j*pi.*x(m)/(2*a)).*sin(k*pi.*y(m)/(2*b)).*cos(k*pi.*z(m)/(2*c));
        R(m,i+5)=d*cos(j*pi.*x(m)/(2*a)).*sin(k*pi.*y(m)/(2*b)).*cos(k*pi.*z(m)/(2*c));
        R(m,i+6)=d*sin(j*pi.*x(m)/(2*a)).*cos(k*pi.*y(m)/(2*b)).*cos(k*pi.*z(m)/(2*c));
        R(m,i+7)=d*cos(j*pi.*x(m)/(2*a)).*cos(k*pi.*y(m)/(2*b)).*cos(k*pi.*z(m)/(2*c));
    end
    i = i + 8;
end
end
end

AA=R\w;

A = gather(AA);
% ----------
% evaluation of surface elevation and calculation of maximum deviation
% ----------
zfitted=R*A;
Maximum_Deviation=max(abs(zfitted-w));
dev = sqrt(mean((zfitted-w).^2));
disp(['Average Deviation: ' num2str(average_dev)])
% ----------
% meshing and plotting of results --colormap Deviation
% ----------
% 		figure;
%         subplot(1,2,1)
%         xline=linspace(min(x),max(x),100);
% 		yline=linspace(min(y),max(y),100);
% 		[X,Y]=meshgrid(xline,yline);
% 		Z=griddata(x,y,w,X,Y,'linear');
% 		Zfitted=griddata(x,y,zfitted,X,Y,'linear');
% 		DevColour=(Z-Zfitted);
% 		surf(X,Y,Zfitted,DevColour,'EdgeColor','white','EdgeAlpha',0,'FaceColor','interp');
% 		axis on; axis tight; grid on; colorbar('horiz'); view(3);
%         % title('Fourier Series fitted Surface-Colour by Deviation');
%         xlabel('a  (x)'); ylabel('b  (y)'); zlabel('w  (z)')
%         subplot(1,2,2)
%         surf(X,Y,Z,DevColour,'EdgeColor','white','EdgeAlpha',0,'FaceColor','interp');
% 		axis on; axis tight; grid on; colorbar('horiz'); view(3);
%         % title('Fourier Series fitted Surface-Colour by Deviation');
%         xlabel('a  (x)'); ylabel('b  (y)'); zlabel('w  (z)')
%         title(['Avg Dev: ' num2str(average_dev)])
% ----------
% end-of-file
% ----------
        