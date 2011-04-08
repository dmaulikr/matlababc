function value = theta_t(s,T,t)

m = -1.0/(T*(1-s));

if abs(t) > 0.5*T
    value = 0.0;
elseif ((abs(t) >= 0) && (abs(t) <=(s-0.5)*T))
    value = 1.0;
else
    value = m*(abs(t)-0.5*T);
end