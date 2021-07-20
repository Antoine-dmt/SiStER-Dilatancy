%Quick script for topo visualization

load 50
plot(topo_x/1e3, topo_y/1e3-10,'g.','LineWidth',8)

axis ij
axis equal
hold on
load 90
plot(topo_x/1e3, topo_y/1e3-10,'k.','LineWidth',8)
load 140
plot(topo_x/1e3, topo_y/1e3-10,'r.','LineWidth',8)
set(gca,'FontSize',20)
xlabel('Distance from left side (km)','FontSize',18)
ylabel('Depth (km)','FontSize',18)
axis ([40 110 -5 10])
legend('A','B','C')
