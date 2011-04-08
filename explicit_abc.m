function [Etb,Hrb,Hzb] = explicit_abc(ecurt,mcurr,mcurz,t_start,t_end,c,dr,re,rh)

ratio=re(end,end)/(re(end,end)+dr);
dt = (t_end-t_start)/20;
cfl= c*dt/dr;

if (cfl > 0.5)
    cfl = 0.5;
    dt = cfl*dr/c;
end

e=size(ecurt); Etb=zeros(e(1)+2,e(2)+2);
r=size(mcurr); Hrb=zeros(r(1)+2,r(2)+2);
z=size(mcurz); Hzb=zeros(z(1)+2,z(2)+2);

Etb(2:end-1,2:end-1) = -4*pi*dt*ecurt;
Hrb(2:end-1,2:end-1) = -4*pi*dt*mcurr;
Hzb(2:end-1,2:end-1) = -4*pi*dt*mcurz;

time = t_start+dt;
while(time <= t_end)
    [Etb,Hrb,Hzb] = periodicINz(Etb,Hrb,Hzb);
    [Etb,Hrb,Hzb] = r_ZERObc(Etb,Hrb,Hzb,cfl);
    [Etb,Hrb,Hzb] = OUTbc(Etb,Hrb,Hzb,ratio);
    Etb = electricSTEP(Etb,Hrb,Hzb,dt,cfl,0);
    [Hrb,Hzb] = magneticSTEP(Etb,Hrb,Hzb,re,rh,dt,cfl,0,0);
    time = time+dt;
end