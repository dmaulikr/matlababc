%
%
%
%
%

clear all; clc;

cfl=0.5;

% grid
r = 0:0.05:1;  lr = length(r);
z = -1:0.05:1; lz = length(z);
[R,Z] = meshgrid(r,z);

% initialization
w = csc(R).*cos(Z);

% computational domain
sol = zeros(lz+2,lr+2); sol(2:end-1,2:end-1)=w;

for k=1:100
% boundary conditions
% r = 0
sol(:,1) = sol(:,2);
% sol(:,1) = 0;
% r = R
sol(2:end-1,end) = sol(2:end-1,end-1);
% Top
sol(end,:) = sol(2,:);
% Bottom
sol(1,:) = sol(end-1,:);
% advance in time
sol(2:end-1,1:end-1) = sol(2:end-1,1:end-1)+cfl*(sol(2:end-1,2:end)-...
    2*sol(2:end-1,1:end-1)+sol(3:end,1:end-1)); 
surf(sol), pause(0.1)
end