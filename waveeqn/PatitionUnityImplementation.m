% ABC implementation with:
%  - partition of unity
%  - 1D wave equation initialized with a pulse

clear all; clc

% constants
R = 10;

% Courant-Frederick Levi stability number
cfl = 0.5;
% wave propagation speed
c = 1;

% grid
dr = 0.01; dt = cfl*dr/c;
steps = 5000;
r = 0:dr:R;
% initial condition
w = exp(-(r-5).^2);

%% Auxiliary domain
tran  = 50;
sigma = 3/4;
T = 1.1;
delta = 1.5*T;
T_int = (R+2*T*c)/c;
terms = fix(T_int/(sigma*T))+2;
%abcR = R+delta+c/2*T_int;
abcR = delta+c/2*T_int;
aux_r = 0:dr:abcR;

%% Multiplier
a=0.05;
mu = step5(0:dr:R+dr,R-tran*dr,R,0,1);

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

% allocate memory for n-1 and n solutions
sol_old = zeros(length(r)+1,1);
sol_new = zeros(length(r)+1,1);
% take one step in time using advection equation
% to account for n-1 solution
sol_old(1:end-1)=w; sol_old(end)=sol_old(end-1); %right boundary condition
sol_new(1:end-1) = sol_old(1:end-1) - cfl*(sol_old(2:end)-sol_old(1:end-1));
%right boundary condition
sol_new(end) = sol_new(end-1);

% Exterior domain memory allocations
subs = zeros(length(aux_r)+1,4*terms);
subs_new = subs; subs_old = subs;
solution = zeros(length(aux_r)+1,1);
aux_src = zeros(length(aux_r)+1,1);
cont_old=aux_src; cont_new=aux_src; cont_src=aux_src;

% loop over timesteps
for k=1:steps
    % current time
    time = (k+1)*dt;
    
    % lower bondary of the solution summation
    p1 = fix(((time-T_int)/T+0.5)/sigma);
    % upper boundary of the solution summation
    p2 = fix((time/T+0.5)/sigma)+1;
    
    %interior solution update
    sol = -sol_old + cfl^2*A*sol_new;
    %transition region source calculation
    src = mu.*(sol+sol_old)-cfl^2*A*(mu.*sol_new);
    subplot(4,1,2), plot(r,src(1:end-1),'g'), title('Source')
    sol_old = sol_new;
    sol_new = sol;
    
    %exterior solution update
    if k==1
       cont_new(1:tran+1) = src(end-tran:end);
       cont = - cont_old + cfl^2*Ae*cont_new;
    else
       cont_src(1:tran+1) = src(end-tran:end);
       cont = - cont_old + cfl^2*Ae*cont_new + cont_src;
    end    
    cont_old = cont_new;
    cont_new = cont;
    
    %i=fix(time)-fix(terms/2+1):fix(time)+fix(terms/2+1);
    %aux_src(1:tran+1) = src(end-tran:end);
    for m=p1:p2+2
        if m > 0
            n = m;
        else 
            n = 1;
        end
        theta = theta_t(sigma,T,time-sigma*T*n);
        subs(:,n) = - subs_old(:,n) + cfl^2*Ae*subs_new(:,n) + theta*cont_src;
        subs_old(:,n) = subs_new(:,n);
        subs_new(:,n) = subs(:,n);
    end
    if p1 > 2
       if (fix(((time-T_int)/T+0.5)/sigma) == fix(((time-T_int-dt)/T+0.5)/sigma)+1)
           solution = 0;
           for sss=p1:p2
               solution = solution + subs(:,sss);
           end
           cont_old = solution;
           cont_new(1:end-1) = solution(1:end-1) - ...
               cfl*(solution(2:end)-solution(1:end-1));
           sol_new(end) = sol_new(end-1);
           cont_old(1) = solution(2);
           
           subplot(4,1,4), plot(solution,'m'), title('Remaining After Removing Source')
           k
       end
   end
    % apply ABC
    sol_new(end) = cont_new(tran+1);
    
    subplot(4,1,1), plot(r,sol(1:end-1),'b')%, title('Interior Solution'), axis([0 10 -1 1])
    subplot(4,1,3), plot(aux_r,cont(1:end-1),'r'), title('Exterior Solution')
    drawnow
end