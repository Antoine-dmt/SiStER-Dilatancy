%Quick script to produce quality ep images with isoT visible on it
% A. Demont 5/2021

clear all
close all
run SiStER_Input_File_oceanic_core_complex.m
snapshot = input('Enter snapshot number');
fluids_depth = 5; % in km
load(num2str(snapshot))
fastscatter(xm(im>1)/1e3,ym(im>1)/1e3-10,xim(im>1));axis ij; axis equal
hold on
axis ([40 110 -5 20])
plot(topo_x/1e3,topo_y/1e3-10,'r.','MarkerSize',17)% surface markers
plot(xm(abs(Tm-400)<2)/1e3,ym(abs(Tm-400)<2)/1e3-10,'g.')% Fluid T limit, and serpentinization approximative temparature limit
plot(xm(abs(Tm-650)<2)/1e3,ym(abs(Tm-650)<2)/1e3-10,'y.')%BDT
plot(topo_x/1e3,topo_y/1e3-10+fluids_depth,'c.','MarkerSize',8)% depth limit for fluid percolation
caxis([0 1])
colorbar
%legend('ep/ecrit','Topography','400°C','600°C','Fluids limit')
set(gca,'FontSize',20)
xlabel('Distance to left side (km)')
ylabel('Depth from seafloor (km)')
