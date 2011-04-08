function s = th(k,s1)

d = k/s1;

if k > s1
    s = 1;
else
    s = d^3*(10-15*d+6*d^2);
end

