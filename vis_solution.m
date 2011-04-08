function vis_solution(x,y,Et,xr,yr,Hr,xz,yz,Hz)

subplot(1,3,1), pcolor(x ,y ,Et), shading interp, axis equal
subplot(1,3,2), pcolor(xr,yr,Hr), shading interp, axis equal
subplot(1,3,3), pcolor(xz,yz,Hz), shading interp, axis equal

pause(0.1)