% 1D-wave equation with variation only on the r-direction
%
clear all; clc

% Computational domain size
R = 10;
% Courant-Frederick Levi stability number
cfl = 0.5;
% wave propagation speed
c = 1;
% grid
dr = 0.01; dt = cfl*dr/c;
% time cycles
steps = 15000;
% frequency of reintegration in iterations
rcyc = 200;
%domain size of the reintergration, interations
range = 2*rcyc;
%number of times reintegrations will occur
num = fix(steps/rcyc);
% -------------------------------
% Domain of interior problem
r = 0:dr:R;
% initial condition
w = exp(-(r-5).^2);

%% Auxiliary domain
abcR = 2*R;
aux_r = 0:dr:abcR;

%% Multiplier
a=0.05;
mu = step5(0:dr:R+dr,0.98*R,R,0,1);

%% Build Matrix
% interior problem
e = ones(length(r)+1,1);
val = 2*(1-1/cfl^2);
A = spdiags([e -val*e e], -1:1, length(e), length(e));
A(end,end-1)=0; A(end,end)=1;
A(1,1:2) = 0;

% exterior problem
ee = ones(length(aux_r)+1,1);
Ae = spdiags([ee -val*ee ee], -1:1, length(ee), length(ee)); 
Ae(end,end-1:end)=0; % no flux upper BC
Ae(1,1:2) = 0;

% take one step in time
sol_old = zeros(length(r)+1,1);
sol_old(1:end-1)=w; sol_old(end)=sol_old(end-1);
sol_new = zeros(length(r)+1,1);
sol_new(1:end-1) = sol_old(1:end-1) - cfl*(sol_old(2:end)-sol_old(1:end-1));
sol_new(end) = sol_new(end-1);

% Retarded sources
aux_src = zeros(length(aux_r)+1,1);
cont_old=aux_src; cont_new=aux_src;
aux_new=aux_src;

% things to be removed
inds = zeros(num,1);
new = zeros(length(aux_r)+1,1);
old = new; rem = new; source = new;
src2 = zeros(length(aux_r)+1,steps);

%%
for k=1:steps
    
    sol = -sol_old + cfl^2*A*sol_new;
    src = mu.*(sol+sol_old)-cfl^2*A*(mu.*sol_new);
    src2(1:101,k) = src(end-100:end);
    sol_old = sol_new;
    sol_new = sol;
    
    subplot(4,1,2), plot(src,'m')    

%  Advance auxiliary problem
   cont = - cont_old + cfl^2*Ae*cont_new + src2(:,k);
   cont_old = cont_new;
   cont_new = cont;
   
   n = fix(k/rcyc)+1;
   inds(n) = k;
   
   % do the reintegration
%    if n > 2
%        if (n+1 == fix((k+1)/rcyc)+1)
%            rem(:) = 0; old(:) = 0; new(:) = 0;
%            for z=k-range:k
%                rem = - old + cfl^2*Ae*new + src2(:,z);
%                old = new;
%                new = rem;
%            end
%            subplot(4,1,4), plot(rem,'g')
%            %if n == 10
%            %    pause
%            %end
%            cont_new = rem;
%            cont_old = old;
%            k
%        end
%    end
   
   % apply ABC
   sol_new(end) = cont_new(101);

   % plot solutions
   subplot(4,1,1), plot(sol), axis([0,1001,-5e-3,5e-3])
   subplot(4,1,3), plot(cont,'r')%, axis([0,1501,-1,1])
   drawnow
end