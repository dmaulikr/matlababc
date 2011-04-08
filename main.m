% Application of Tsynkov Lacunae based Artificial boundary conditions
% 
% Applications described in: Artificial Boundary Conditions for the
% Numerical Simulation of Unsteady Electromagnetic Waves by S.V Tsynkov
% ---- unpublished report

clear all; clc;

%% Constants
c  = 1;      % speed of light
cfl= 0.5;    % Courant number

%% Construction of the Computational Grid
% nodes
Nr = 64;
Nz = 2*Nr;
% radius
R  = pi;
dr = R/Nr;
ratio = R/(R+dr);
% length
Z  = 2*pi;
dz = Z/Nz;
% timesteping
dt = min (cfl/c*[dr dz]);

%% Interior domain
% electric field grid
r = linspace(0,R,Nr);
z = linspace(-Z/2,Z/2,Nz);
[x,y] = meshgrid(r,z);
% r magnetic field grid
rr = linspace(0,R,Nr);
zr = linspace(-Z/2+dz/2,Z/2-dz/2,Nz-1);
[xr,yr] = meshgrid(rr,zr);
% z magnetic field grid
rz = linspace(0+dr/2,R-dr/2,Nr-1);
zz = linspace(-Z/2,Z/2,Nz);
[xz,yz] = meshgrid(rz,zz);

%clear r z rr zr rz zz
% -------------------------------------------------------------------------
%% -------- Initialization --------------------
Et = cos(pi/R*(x-(R/2))).*sin(1*pi/Z*y); Et=Et/max(max(Et));
Hr = -sin(2*pi/R*xr).*cos(0.5*pi/Z*yr); Hr=Hr/max(max(Hr));
Hz = -2*Z/pi*sin(pi/(2*Z)*yz).*(cos(2*pi/R*xz)+...
    sin(2*pi/R*xz)./xz);Hz=Hz/max(max(Hz));

% -------------------------------------------------------------------------
% allocate memory for interior solution and initialize it
SolEt = zeros(Nz+2,Nr+2); SolEt(2:end-1,2:end-1) = Et;
SolHr = zeros(Nz+1,Nr+2); SolHr(2:end-1,2:end-1) = Hr;
SolHz = zeros(Nz+2,Nr+1); SolHz(2:end-1,2:end-1) = Hz;

% -------------------------------------------------------------------------
%% Create transition region
% size of transition region
tran = 5*dr;
% parameters needed for lacunae integration
delta=3/2*dr; sigma=3/4; T=dt;%9/16*pi;
% interval that the waves will remain in the domain
T_int = (R+2*T*c)/c;
% terms in lacuna integration
terms = fix(T_int/(sigma*T))+2;
%% Allocate memory for auxiliary domain
abcR = R+delta+c/2*T_int;
abcZ = Z+2*delta+c*T_int;
st = length(-abcZ/2:dz:abcZ/2);
st = st/2-Nr;
% theta electric field grid
[abcEx,abcEy] = meshgrid(0:dr:abcR,-abcZ/2:dz:abcZ/2);
% r magnetic field grid
[abcHr_x,abcHr_y] = meshgrid(0:dr:abcR,-abcZ/2+dz/2:dz:abcZ/2-dz/2);
% z magnetic field grid
[abcHz_x,abcHz_y] = meshgrid(0+dr/2:dr:abcR-dr/2,-abcZ/2:dz:abcZ/2);

%% transition multiplier (mu)
a = 0.05;
m1 = 1./(1+exp(-(r-(R-tran))/a));
m2 = 1./(1+exp(-(rz-(R-tran))/a));
muE = ones(size(Et));
muHr= ones(size(Hr));
muHz= ones(size(Hz));
for l=1:length(z)-1
    muE(l,:) = m1;
    muHr(l,:) = m1;
end
muE(length(z),:) = m1;
for n=1:length(zz)
    muHz(n,:) = m2;
end

%% alocate space for source terms
ecurt = zeros(size(abcEx));
mcurr = zeros(size(abcHr_x));
mcurz = zeros(size(abcHz_x));

%% allocate memory for exterior solution values to be subtracted
SolEtb = zeros(length(-abcZ/2:dz:abcZ/2),length(0:dr:abcR),terms);
SolHrb = zeros(length(-abcZ/2+dz/2:dz:abcZ/2-dz/2),length(0:dr:abcR),terms);
SolHzb = zeros(length(-abcZ/2:dz:abcZ/2),length(0+dr/2:dr:abcR-dr/2),terms);

%% allocate memory for the exterior solution explicit integration
auxEt = zeros(length(-abcZ/2:dz:abcZ/2),length(0:dr:abcR));
auxHr = zeros(length(-abcZ/2+dz/2:dz:abcZ/2-dz/2),length(0:dr:abcR));
auxHz = zeros(length(-abcZ/2:dz:abcZ/2),length(0+dr/2:dr:abcR-dr/2));

p1 = 0:1:10;
p1_time = T*(sigma*p1-1/2)+T_int;
 
%% Step the problem forward
for k=1:150
    % current time
    time = k*dt;
    % Apply Boundary Conditions
    [SolEt,SolHr,SolHz] = periodicINz(SolEt,SolHr,SolHz);
    [SolEt,SolHr,SolHz] = r_ZERObc(SolEt,SolHr,SolHz,cfl);
    [SolEt,SolHr,SolHz] = wguide(SolEt,SolHr,SolHz);
%    [SolEt,SolHr,SolHz] = OUTbc(SolEt,SolHr,SolHz,ratio);
    vis_solution(x,y,SolEt(2:end-1,2:end-1),xr,yr,...
        SolHr(2:end-1,2:end-1),xz,yz,SolHz(2:end-1,2:end-1))
        
    % Take a Step
    SolEtn = electricSTEP(SolEt,SolHr,SolHz,dt,cfl,0);
    [SolHrn,SolHzn] = magneticSTEP(SolEt,SolHr,SolHz,x,xz,dt,cfl,0,0);
    
%% Calculate the sources
    % source: electric current
    ecurt(st+1:st+Nz-1,1:Nr-1) = -(c/(4*pi))*(1/(c*dt)*...
     muE(1:end-1,1:end-1).*(SolEtn(2:end-2,2:end-2)-SolEt(2:end-2,2:end-2))- ...
     (1/dz*muHr(:,1:end-1).*(SolHrn(2:end-1,2:end-2)-SolHrn(1:end-2,2:end-2))- ...
     1/dr*muHz(1:end-1,:).*(SolHzn(2:end-2,2:end-1)-SolHzn(2:end-2,1:end-2))));

    % source: magnetic currents
    mcurr(st+1:st+Nz-1,1:Nr-1) = -(c/(4*pi))*(1/(c*dt)*...
     muHr(:,1:end-1).*(SolHrn(2:end-1,2:end-2)-SolHr(2:end-1,2:end-2))- ...
     (1/dz*muE(1:end-1,1:end-1).*(SolEtn(3:end-1,2:end-2)-SolEtn(2:end-2,2:end-2))));

    mcurz(st+1:st+Nz-1,1:Nr-1) = -(c/(4*pi))*(1/(c*dt)*...
     muHz(1:end-1,:).*(SolHzn(2:end-2,2:end-1)-SolHz(2:end-2,2:end-1))- ...
     (1/dr*muE(1:end-1,1:end-1)./xz(1:end-1,:).*(x(1:end-1,2:end).*...
     SolEtn(2:end-2,3:end-1)-x(1:end-1,1:end-1).*SolEtn(2:end-2,2:end-2))));

 % Advance the exterior solution expliticly
[auxEt,auxHr,auxHz] = periodicINz(auxEt,auxHr,auxHz);
[auxEt,auxHr,auxHz] = r_ZERObc(auxEt,auxHr,auxHz,cfl);
[auxEt,auxHr,auxHz] = OUTbc(auxEt,auxHr,auxHz,ratio);
auxEt = electricSTEP(auxEt,auxHr,auxHz,dt,cfl,ecurt(2:end-1,3:end-1));
[auxHr,auxHz] = magneticSTEP(auxEt,auxHr,auxHz,...
    abcEx(2:end-1,2:end-1),abcHz_x(2:end-1,2:end-1),dt,cfl,...
    mcurr(2:end,2:end-1),mcurz(2:end-1,1:end-3));
 
i=fix(time)-3:fix(time)+3;
for m=1:length(i)
    if (i(m) >=0)
        val = theta_t(sigma,T,time-sigma*T*i(m));
        if (val ~= 0)
            ect = ecurt*val;
            mcr = mcurr*val;
            mcz = mcurz*val;
            % Advance source to substraction time
            [Etb,Hrb,Hzb] = explicit_abc(ect,mcr,mcz,time,...
                p1_time(i(m)+1),c,dr,abcEx,abcHz_x);
            % Add contributions that are to be subtract at the same time
            % this is possible due to superposition in the exterior domain
            SolEtb(:,:,i(m)+1) = SolEtb(:,:,i(m)+1) + Etb(2:end-1,2:end-1);
            SolHrb(:,:,i(m)+1) = SolHrb(:,:,i(m)+1) + Hrb(2:end-1,2:end-1);
            SolHzb(:,:,i(m)+1) = SolHzb(:,:,i(m)+1) + Hzb(2:end-1,2:end-1);
        end
    end
end
%% Switch condition which subtracts the retarded sources
% lower boundary of summation
   p1 = fix(((time-T_int)/T+0.5)/sigma);
if p1 >= 0
    if (fix(((time+dt-T_int)/T+0.5)/sigma) == fix(((time-T_int)/T+0.5)/sigma)+1)
        auxEt = auxEt-SolEtb(:,:,p1+1);
        auxHr = auxHr-SolHrb(:,:,p1+1);
        auxHz = auxHz-SolHzb(:,:,p1+1);
    end
end
SolEtn(2:end-1,end)=auxEt(st+1:st+2*Nr,Nr+1);
SolHrn(2:end,end)=auxHr(st+1:st+2*Nr,Nr+1);
SolHzn(2:end-1,end)=auxHz(st+1:st+2*Nr,Nr+1);

% copy
SolEt = SolEtn; SolHr=SolHrn; SolHz=SolHzn;
clc
k
end