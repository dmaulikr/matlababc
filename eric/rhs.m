% -------------------------------------------------------------------------
% function rhs
% define right-hand side for 1D Euler system
% -------------------------------------------------------------------------
function rhsx_vec = rhs(q_vec)
global h nqty npi npa npt mu St trhs Q
cput0=cputime;
% -------------------------------------------------------------------------
% reshape solution vector
% -------------------------------------------------------------------------
q(:,:) = reshape(q_vec,nqty,npi+npa);
qt(:,:) = q(:,npi-npt+1:npi);
muqt = zeros(nqty,npt);
for iqty = 1:nqty
    muqt(iqty,:) = mu'.*(qt(iqty,:) - Q(iqty)) + Q(iqty);
end

% -------------------------------------------------------------------------
% prepare quantities for matrix formation
% -------------------------------------------------------------------------
d1i = ones(npi,1)/(2*h);
d1a = ones(npa,1)/(2*h);
d1t = ones(npt,1)/(2*h);
Za = sparse(npi,npa);
Zi = sparse(npa,npi);

% -------------------------------------------------------------------------
% form operator matrix for centered first derivative (2nd order accurate)
% -------------------------------------------------------------------------
D1i = spdiags([-d1i d1i],[-1 1],npi,npi);
D1a = spdiags([-d1a d1a],[-1 1],npa,npa);
D1t = spdiags([-d1t d1t],[-1 1],npt,npt);

D1 = cat(1,cat(2,D1i,Za),cat(2,Zi,D1a));

% -------------------------------------------------------------------------
% compute flux
% -------------------------------------------------------------------------
F = zeros(nqty,npi+npa);
Ft = zeros(nqty,npt);
Fta = zeros(nqty,npt);
for iqty = 1:nqty
    F(iqty,:) = getF(q,iqty);
    Ft(iqty,:) = getF(qt,iqty);
    Fta(iqty,:) = getF(muqt,iqty);
end

% -------------------------------------------------------------------------
% compute source
% note that mu, mux, muxx, and thus St are zero at St(1) and St(npt)
% -------------------------------------------------------------------------
St = zeros(nqty,npt);
for iqty = 1:nqty
     St(iqty,:) = D1t*Fta(iqty,:)' - mu.*(D1t*Ft(iqty,:)');
end
St(:,1) = 0;
St(:,npt) = 0;

S = zeros(nqty,npi+npa);
S(:,npi+1:npi+npt) = St;

% -------------------------------------------------------------------------
% compute rhs
% -------------------------------------------------------------------------
rhsx = zeros(nqty,npi+npa);
for iqty = 1:nqty
    rhsx(iqty,:) = -D1*F(iqty,:)' + S(iqty,:)';
end
rhsx_vec = reshape(rhsx,nqty*(npi+npa),1);

cput1 = cputime;
trhs = trhs + cput1 - cput0;