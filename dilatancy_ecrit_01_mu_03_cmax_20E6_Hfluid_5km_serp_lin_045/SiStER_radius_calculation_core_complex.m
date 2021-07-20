%Script for elevation extraction from GeoMap App profiles and hand picking
%for curve fitting

%Profile loading
clear
load 140
%plot(topo_x,topo_y)

%plot(dist_13N20,elev_13N20)
begin_fit_point_13N20 = 66914; % manually picked on curve
end_fit_point_13N20 = 75926;%here end of first core complex 
fit_dist_13N20 = topo_x(find(topo_x>begin_fit_point_13N20 & topo_x<end_fit_point_13N20)); %fitting topography
fit_elev_13N20 = topo_y(topo_x>begin_fit_point_13N20 & topo_x<end_fit_point_13N20);
%plot(fit_dist_13N20,fit_elev_13N20)

%fitting now curvature radius
fit_deg = 2; %degreee for polynomial fitting
polyfit_curve_13N20 = polyfit(fit_dist_13N20,fit_elev_13N20, fit_deg);
x1 = fit_dist_13N20(1):1:fit_dist_13N20(end);
y1 = polyval(polyfit_curve_13N20,x1);
plot(fit_dist_13N20,fit_elev_13N20,'o',x1,y1,'--')

grad_first_13N20 = gradient(y1(:)) ./ gradient(x1(:));
curvature_array_13N20 = gradient(grad_first_13N20(:)) ./ gradient(x1(:));
radius_array_13N20 = 1./ (abs(curvature_array_13N20));
radius_average_13N20 = mean(radius_array_13N20);
