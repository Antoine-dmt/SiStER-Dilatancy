clear all
run SiStER_Input_File_oceanic_core_complex.m
for k = 10:10:800
    load(num2str(k))
    fastscatter(xm(im>1),ym(im>1),ep(im>1))
    %fastscatter(xm(ep>MAT(2).ecrit),ym(ep>MAT(2).ecrit),ep(ep>MAT(2).ecrit))
    %caxis([-200 200])
    axis ij
    axis equal
    
    colorbar
    %xlim([7e4 9e4])
    %ylim([1.2e4 2e4])
    %grid on
    %caxis([0 1])  
    pause(0.3)
end
%hold on
%quiver(X,Y,vx,vy)