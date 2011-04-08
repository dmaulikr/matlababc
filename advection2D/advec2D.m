clear all; clc

cfl = 0.5;

x=-3:0.01:3;
n = length(x);

[X,Y]=meshgrid(x,x);
w = exp(-X.^2-Y.^2);
sol = zeros(n,n+1); sol(:,1:end-1)=w;

for k=1:1000
    sol(:,end) = sol(:,2);
    sol(:,1:end-1) = (1-cfl)*sol(:,2:end) + cfl*sol(:,1:end-1);
    pcolor(X,Y,sol(:,1:end-1)), shading interp
    drawnow
end