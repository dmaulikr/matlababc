% -------------------------------------------------------------------------
% check ierr for NKS
% -------------------------------------------------------------------------
function chkerr(ierr)

if ierr == 1
    NKS_stop_maxit
elseif ierr == 2
    NKS_stop_LS_fail
end