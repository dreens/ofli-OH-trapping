function A=fitfourier3D(data,n)

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
tic
x = data(:,1);
y = data(:,2);
z = data(:,3);
w = data(:,4);

datasize=size(data);
a=abs(max(x)-min(x)); 
b=abs(max(y)-min(y));
c=abs(max(z)-min(z));
R = zeros(length(x),n^3+1);
R(:,1) = max(abs(w))*ones(size(x));
i = 2;
for j=1:n
for k=1:n
for l=1:n
    Rjk1=max(abs(w))*sin(j*pi.*abs(x)/(2*a)).*sin(k*pi.*abs(y)/(2*b)).*sin(k*pi.*z/(2*c));
    Rjk2=max(abs(w))*cos(j*pi.*abs(x)/(2*a)).*cos(k*pi.*abs(y)/(2*b)).*sin(k*pi.*z/(2*c));
    Rjk3=max(abs(w))*sin(j*pi.*abs(x)/(2*a)).*cos(k*pi.*abs(y)/(2*b)).*sin(k*pi.*z/(2*c));
    Rjk4=max(abs(w))*cos(j*pi.*abs(x)/(2*a)).*sin(k*pi.*abs(y)/(2*b)).*sin(k*pi.*z/(2*c));
    Rjk5=max(abs(w))*sin(j*pi.*abs(x)/(2*a)).*sin(k*pi.*abs(y)/(2*b)).*cos(k*pi.*z/(2*c));
    Rjk6=max(abs(w))*cos(j*pi.*abs(x)/(2*a)).*cos(k*pi.*abs(y)/(2*b)).*cos(k*pi.*z/(2*c));
    Rjk7=max(abs(w))*sin(j*pi.*abs(x)/(2*a)).*cos(k*pi.*abs(y)/(2*b)).*cos(k*pi.*z/(2*c));
    Rjk8=max(abs(w))*cos(j*pi.*abs(x)/(2*a)).*sin(k*pi.*abs(y)/(2*b)).*cos(k*pi.*z/(2*c));
    R(:,i:i+7)=[Rjk1 Rjk2 Rjk3 Rjk4 Rjk5 Rjk6 Rjk7 Rjk8];
    i = i + 8;
end
end
end
A=R\w;
% ----------
% evaluation of surface elevation and calculation of maximum deviation
% ----------
zfitted=R*A;
Maximum_Deviation=max(abs(zfitted-w));
average_dev = sqrt(mean((zfitted-w).^2));
disp(['Average Deviation: ' num2str(average_dev)])
toc
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
        