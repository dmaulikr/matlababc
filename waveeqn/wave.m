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
steps = 20000;
% remove after x number of cycles
rcyc = 400; s1 = fix(rcyc/5);
num = fix(steps/rcyc);
% -------------------------------
r = 0:dr:R;
% initial condition
w = exp(-(r-5).^2);

%% Auxiliary domain
abcR = 2*R;
aux_r = 0:dr:abcR;

%% Multiplier
a=0.05;
mu = step5(0:dr:R+dr,0.9*R,R,0,1);

%% Build Matrix
% interior problem
e = ones(length(r)+1,1);
A = spdiags([(1-cfl)*e -cfl*e], 0:1, length(e), length(e));

% exterior problem
ee = ones(length(aux_r)+1,1);
Ae = spdiags([(1-cfl)*ee -cfl*ee], 0:1, length(ee), length(ee)); 

% initialize
sol = zeros(length(r)+1,1); sol_new = sol;
sol(1:end-1)=w; sol(end)=sol(end-1);

% allocate memory for exterior solution
cont = zeros(length(aux_r)+1,1);

% things to be removed
inds = zeros(num,1);
rem = zeros(length(aux_r)+1,1);
src = zeros(length(aux_r)+1,steps);

%%
for k=1:steps
    
    sol_new(1:end-1) = sol(1:end-1) - cfl*(sol(2:end) - sol(1:end-1));
    
    src(1:length(mu)-1,k) = mu(1:end-1).*(sol_new(1:end-1)-sol(1:end-1))/dt + ...
        (mu(2:end).*sol(2:end) - mu(1:end-1).*sol(1:end-1))/dr;
    sol = sol_new;    
    subplot(4,1,2), plot(src(:,k),'m')    

%  Advance auxiliary problem
   cont(1:end-1) = cont(1:end-1) - cfl*(cont(2:end) - cont(1:end-1)) + ...
       src(1:end-1,k);
   cont(end) = cont(end-1);
   
   n = fix(k/rcyc)+1;
   inds(n) = k;
   
   % do the reintegration
   if n > 2
       if (n+1 == fix((k+1)/rcyc)+1)
           rem(:) = 0;
           for z=inds(n-2):k
               v = th(z,inds(n-2)+s1);
               rem(1:end-1) = rem(1:end-1) - cfl*(rem(2:end) - rem(1:end-1)) +...
                   v*src((1:end-1),z);
               rem(end) = rem(end-1);
           end
           subplot(4,1,4), plot(rem,'g')           
           cont = rem;
       end
   end
   
   % apply ABC
   sol(end) = cont(length(sol));

   % plot solutions
   subplot(4,1,1), plot(sol)
   subplot(4,1,3), plot(cont,'r')
   drawnow
end