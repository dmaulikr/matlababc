function bc = abc1d(u,src,c,R)

sigma = 3/4;
T = 1.1;
delta = 0.5;
T_int = (R+2*T*c)/c;

abcR = R+delta+c/2*T_int;