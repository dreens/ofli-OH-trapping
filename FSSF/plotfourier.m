function output=plotfourier(string,dataset,bcx_0,bcx_end,bcy_0,bcy_end);
% function to regenerate the fitted surface given the calculated coefficients
% and the original dataset
% string=filename (eg. string='coefficients.dat')
% dataset=the original dataset variable
% bc_... boundary conditions as described in corresponding documentation

coeffraw=csvread(string);
dim=size(coeffraw);
index1=1; index2=1;
coefficients=[];
for index1=1:dim(1)
    for index2=1:dim(2)
        coefficients=[coefficients;coeffraw(index1,index2)];
    end
end
x=dataset(:,1); y=dataset(:,2); z=dataset(:,3);
a=abs(max(x)-min(x)); b=abs(max(y)-min(y));
j=1; k=1; R=[];
    for j=1:2:dim(1)
        for k=1:2:dim(2)
            Rjk=max(abs(z))*sin(j*pi.*x/(bcx_end*a)+(bcx_0-1)*pi/(bcx_end*2)).*sin(k*pi.*y/(bcy_end*b)+(bcy_0-1)*pi/(bcy_end*2));  %Fourier Series
            R=[R Rjk];
        end
    end
    zfitted=R*coefficients;
    xline=linspace(min(x),max(x),50);
    yline=linspace(min(y),max(y),50*floor(abs((max(x)-min(x))/(max(y)-min(y)))));
    [X,Y]=meshgrid(xline,yline);
    Z=griddata(x,y,z,X,Y,'v4');
    Zfitted=griddata(x,y,zfitted,X,Y,'v4');
    DevColour=(Z-Zfitted);
    surface(X,Y,Zfitted,DevColour,'EdgeColor','white','EdgeAlpha',0.5,'FaceColor','interp');
    axis on; axis vis3d; axis equal; axis tight; grid on; colorbar('horiz'); view(3);
    title('Fourier Series fitted Surface-Colour by Deviation');
    xlabel('a  (x)'); ylabel('b  (y)'); zlabel('w  (z)')
    output.X=X;
    output.Y=Y;
    output.Z=Zfitted;
    output.C=DevColour;