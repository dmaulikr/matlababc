% 1D-wave equation with variation only on the r-direction
%
clear all; clc

% constants
%f = 0.01;   %frequency
R = 10;

% Courant-Frederick Levi stability number
cfl = 0.75;
% wave propagation speed
c = 1;

% grid
dr = 0.01; dt = cfl*dr/c;
steps = 6000;
r = 0:dr:R;
% initial condition
w = exp(-(r-5).^2);

%% Auxiliary domain
tran  = 100*dr;
sigma = 3/4;
T = 0.1;
delta = 1.5*dr;
T_int = (R+2*T*c)/c;
terms = fix(T_int/(sigma*T))+2;
%abcR = 25;
abcR = R+delta+c/2*T_int;
aux_r = 0:dr:abcR;

%% Multiplier
a=0.05;
mu = step5(0:dr:R+dr,R-tran,R,0,1);

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
sol_new(1:end-1) = sol_old(1:end-1) + cfl*(sol_old(2:end)-sol_old(1:end-1));
sol_new(end) = sol_new(end-1);

% Retarded sources
subs = zeros(length(aux_r)+1,4*terms);
aux_src = zeros(length(aux_r)+1,1);
cont_old=aux_src; cont_new=aux_src; cont_src=aux_src;
aux_new=aux_src;
%mu2 = step5(0:dr:abcR+dr,14,15,0,-1);
%mu2(1:length(r)) = 0;
src2 = aux_src; m=1;

src = zeros(length(mu),steps);

%m = 1:2*terms;
%tstop = T*(sigma*m-1/2)+T_int;

for k=1:steps
    time = k*dt;
    p1 = fix(((time-T_int)/T+0.5)/sigma);
    p2 = fix((time/T+0.5)/sigma)+1;
    sol = -sol_old + cfl^2*A*sol_new;
    src(:,k) = mu.*(sol+sol_old)-cfl^2*A*(mu.*sol_new);
    subplot(3,1,2), plot(r,src(1:end-1,k),'m')    
    sol_old = sol_new;
    sol_new = sol;
%    for j=p1-1:p1+terms
%       if j>=0
%           val = theta_t(sigma,T,time-sigma*T*j);
%           if val > 0
%              aux_src(1:length(src)) = src;
%             aux_old = aux_src;
%              aux_new(1:end-1) = aux_old(1:end-1) + ...
%                  cfl*(aux_old(2:end)-aux_old(1:end-1));
%              t=2*dt; %(j-1+sigma*j)*T
%               while (t <= R/c+T_int)
%                  aux_sol = -aux_old + cfl^2*Ae*aux_new;
%                  aux_old = aux_new;
%                  aux_new = aux_sol;
%                  t = t+dt;
%               end
% if p1 <= 0
%     j = 1;
% else
%     j = p1+1;
% end
% subs(1:length(src),j) = subs(1:length(src),j) + src;
% subplot(3,1,2),plot(subs(:,j),'g')
%           end
%       end
%    end

   if time > 5.11
       %pause
       src2(end-length(mu)+1:end) = src(:,m);
       m = m + 1;
   end

   if k==1
       cont_new(1:length(r)+1) = src(:,k);
       cont = - cont_old + cfl^2*Ae*cont_new;
   else
       cont_src(1:length(r)+1) = src(:,k);
       cont = - cont_old + cfl^2*Ae*cont_new + cont_src - src2;
   end    
    
   cont_old = cont_new;
   cont_new = cont;
   
%       % subtract retarded term
%    clc
%    p1
%    p2
%    if p1 >= 0
%        if (fix(((time-T_int)/T+0.5)/sigma) == fix(((time-T_int-dt)/T+0.5)/sigma)+1)
% %            for m=p1:p2
% %               src2 = src2 + subs(:,m+1);
% %            end
%            aux_old = subs(:,p1+1);
%            aux_new(1:end-1) = aux_old(1:end-1) + ...
%                  cfl*(aux_old(2:end)-aux_old(1:end-1));
%            t=2*dt; %(j-1+sigma*j)*T
%            while (t <= time)
%               aux_sol = -aux_old + cfl^2*Ae*aux_new;
%               aux_old = aux_new;
%               aux_new = aux_sol;
%               t = t+dt;
%            end
%            cont = cont - aux_sol;
%            %cont_src = cont_src - subs(:,p1+1);
%        end
%    end
   
   % apply ABC
   sol_new(end) = cont(length(r)+1);
    
   subplot(3,1,1), plot(r,sol(1:end-1)), axis([0 R -1 1])
   %subplot(3,1,2), plot(r,src(1:end-1),'g')
   subplot(3,1,3), plot(aux_r,cont(1:end-1),'r'), axis([0 abcR -1 1])
   drawnow
end