%
%
%
clear all; clc

x = -2:0.01:2;
y = -2:0.01:2;
[X,Y]=meshgrid(x,y); r = sqrt(X.^2+Y.^2);

w = (sin(7*r)./(7*r)).^2;
for t=0.01:0.01:3
    sol=w./r.*cos(10*t-4*r+pi);
    
    pcolor(X,Y,sol); shading interp
    caxis([-0.01 0.01]);
    colorbar;
    pause(0.1)
end