function [Hr,Hz] = magneticSTEP(Et,Hr,Hz,x,xz,dt,cfl,jmr,jmz)


Hr(2:end,2:end-1) = Hr(2:end,2:end-1) + cfl*(Et(3:end,2:end-1)-...
    Et(2:end-1,2:end-1))-4*pi*dt*jmr;

Hz(2:end-1,3:end-1) = Hz(2:end-1,3:end-1) - cfl./xz(:,2:end).*...
    (x(:,3:end).*Et(2:end-1,4:end-1)-x(:,2:end-1).*Et(2:end-1,3:end-2))-4*pi*dt*jmz;