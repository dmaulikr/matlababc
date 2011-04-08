function mu = step5(x,x0,x1,h0,h1)

mu = zeros(length(x),1);
a = h1-h0;

for i=1:length(x)
    d = (x(i)-x0)/(x1-x0);
    if x(i) <= x0
        mu(i) = h0;
    elseif x(i) >= x1
        mu(i) = h1;
    else
        mu(i) = h0 + a*d^3*(10-15*d+6*d^2);
    end
end

