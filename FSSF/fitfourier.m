function A=fitfourier(data,n)

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

datasize=size(data);
a=abs(max(x)-min(x)); b=abs(max(y)-min(y));
R=max(abs(z))*ones(size(x));
for j=1:n
    for k=1:n
        Rjk1=max(abs(z))*sin(j*pi.*abs(x)/(2*a)).*sin(k*pi.*abs(y)/(2*b));
        Rjk2=max(abs(z))*cos(j*pi.*abs(x)/(2*a)).*cos(k*pi.*abs(y)/(2*b));
        Rjk3=max(abs(z))*sin(j*pi.*abs(x)/(2*a)).*cos(k*pi.*abs(y)/(2*b));
        Rjk4=max(abs(z))*cos(j*pi.*abs(x)/(2*a)).*sin(k*pi.*abs(y)/(2*b));
        R=[R Rjk1 Rjk2 Rjk3 Rjk4];
    end
end
A=R\z;
% ----------
% evaluation of surface elevation and calculation of maximum deviation
% ----------
zfitted=R*A;
Maximum_Deviation=max(abs(zfitted-z));
average_dev = sqrt(mean((zfitted-z).^2));
% ----------
% meshing and plotting of results --colormap Deviation
% ----------
		figure;
        subplot(1,2,1)
        xline=linspace(min(x),max(x),100);
		yline=linspace(min(y),max(y),100);
		[X,Y]=meshgrid(xline,yline);
		Z=griddata(x,y,z,X,Y,'linear');
		Zfitted=griddata(x,y,zfitted,X,Y,'linear');
		DevColour=(Z-Zfitted);
		surf(X,Y,Zfitted,DevColour,'EdgeColor','white','EdgeAlpha',0,'FaceColor','interp');
		axis on; axis tight; grid on; colorbar('horiz'); view(3);
        % title('Fourier Series fitted Surface-Colour by Deviation');
        xlabel('a  (x)'); ylabel('b  (y)'); zlabel('w  (z)')
        subplot(1,2,2)
        surf(X,Y,Z,DevColour,'EdgeColor','white','EdgeAlpha',0,'FaceColor','interp');
		axis on; axis tight; grid on; colorbar('horiz'); view(3);
        % title('Fourier Series fitted Surface-Colour by Deviation');
        xlabel('a  (x)'); ylabel('b  (y)'); zlabel('w  (z)')
        title(['Avg Dev: ' num2str(average_dev)])
% ----------
% end-of-file
% ----------
        