% -------------------------------------------------------------------------
% solver for 1D Euler system
% uses lacunae-based approach as ABC to truncate domain
% dq/dt + dF/dx = 0
% -------------------------------------------------------------------------
clc; clear all; close all;
global q0 dt theta h nqty npi npa npt mu Q St trhs St0 St1 gamma
cput0=cputime;
trhs=0;

% -------------------------------------------------------------------------
% main user inputs are initialized here
% -------------------------------------------------------------------------
theta = .6;

p0 = 1;     % background pressure
rho0 = 1;   % background density
eps = 1;   % perturbation on background pressure
             % if eps is large, the solution will exh5ibit nonlinearity
            
n = 200;                % number of interior cells
nt = 20;                % number of transition cells
na = 3*nt;              % number of auxiliary cells
npi = n + 1;            % number of interior points
npa = na + 1;           % number of auxiliary points
npt = nt + 1;           % number of transition points

nr = 20;    % number of steps to take between reintegration of ABC 
            % source term
nk_fac = 1;

cfl = 2;
T = 1;

% -------------------------------------------------------------------------
% define interior domain
% length, L
% number of divisions, n
% spacing, h
% -------------------------------------------------------------------------
L = 1;
h = L/n;

% -------------------------------------------------------------------------
% time stepping parameters
% courant-friedrichs-lewy limit, cfl
% time step size, dt
% period, T
% number of steps, nstep
% -------------------------------------------------------------------------
gamma = 5/3;
cs = sqrt(gamma*(p0 + eps)/rho0);
dt = cfl*h/cs;
nstep = round(T/dt);

% -------------------------------------------------------------------------
% define an auxiliary region that will contain the interior domain.
% within the interior, a transition region is defined along with an
% associated transition function, mu.
% transition region size, Lt
% number of divisions in the auxiliary region, na
% number of divisions in the transition region, nt
% time for source effects to clear the transition region, tt
% -------------------------------------------------------------------------
Lt = h*nt;
tt = h*nt/cs;
nkeep = round(tt/dt);
nkeep = round(nk_fac*nkeep + nr);
nstep
nkeep
nr
% -------------------------------------------------------------------------
% initialize position arrays and initial condition
% the first npi points are the interior solution.  the next npa points
% represent the auxiliary solution, ua.  npt of these points are actually
% overlapping.  the BC for the interior problem is that the npi point must
% equal the auxiliary solution at the (npt)th point of ua, that is, the
% (npi+npt)th point of the full solution vector, u.
% -------------------------------------------------------------------------
xi = zeros(npi,1);
xa = zeros(npa,1);
xt = zeros(npt,1);

mu = xt;

nqty=3;
Q = zeros(nqty,1);
Q(1) = rho0;
Q(3) = p0/(gamma-1);

St_tmp = zeros(nqty,npt,nstep+1); % storage for St at every step.

q = zeros(nqty,npi+npa,nstep+1);
for iqty = 1:nqty, q(iqty,:,1) = Q(iqty); end

for i = 1:npi
    xi(i) = h*(i - 1);
%     if xi(i) >= L/4 && xi(i) < L/2
%         bump = exp(1 - 1/(1 - ((xi(i) - 3*L/8)/(L/8))^2));
    if xi(i) < L/6
        bump = exp(1 - 1/(1 - (xi(i)/(L/6))^2));
        q(3,i,1) = q(3,i,1) + eps*bump;
    end
end

for i = 1:npa
    xa(i) = h*(i - 1) + L - Lt;
end

for i = 1:npt
    xt(i) = h*(i - 1) + L - Lt;
    % define mu following the STEP5 fortran routine
    mu(i) = ((-L + Lt + xt(i))^3*(6*L^2 + Lt^2 + 3*L*(Lt - 4*xt(i)) ...
        - 3*Lt*xt(i) + 6*xt(i)^2))/Lt^5;
end

% -------------------------------------------------------------------------
% solve using Newton-Krylov solver
% -------------------------------------------------------------------------
tol = [1e-6,1e-6];
ir = nr;
for is = 1:nstep
    if is < nkeep || ir < nr
        q0 = reshape(q(:,:,is),nqty*(npi+npa),1);
        [q1,it_hist,ierr] = NKS(q0,@resid,tol); chkerr(ierr);
        q(:,:,is+1) = reshape(q1,nqty,npi+npa);
        St_tmp(:,:,is+1) = St;
        if is >= nkeep
            ir = ir + 1;
        end
    else
        is
        qa = zeros(nqty,npa,nkeep-nr);
        for iqty = 1:nqty, qa(iqty,:,1) = Q(iqty); end
        for ik = 1:nkeep-nr-1
            imain = is - (nkeep - nr) + ik;
            q0 = reshape(qa(:,:,ik),nqty*npa,1);
            St0 = St_tmp(:,:,imain);
            St1 = St_tmp(:,:,imain+1);
            [q1,it_hist,ierr] = NKS(q0,@resida,tol); chkerr(ierr);
            qa(:,:,ik+1) = reshape(q1,nqty,npa);
        end
        q(:,npi+1:npi+npa,is) = qa(:,:,nkeep-nr);
        q0 = reshape(q(:,:,is),nqty*(npi+npa),1);
        [q1,it_hist,ierr] = NKS(q0,@resid,tol); chkerr(ierr);
        q(:,:,is+1) = reshape(q1,nqty,npi+npa);
        St_tmp(:,:,is+1) = St;
        ir = 1;
    end
end

cput1 = cputime;
computation_time = cput1-cput0 %#ok<NOPTS>
rhs_time = trhs %#ok<NOPTS>