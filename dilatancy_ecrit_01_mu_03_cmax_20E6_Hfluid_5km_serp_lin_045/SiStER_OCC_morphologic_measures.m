
clear all
run SiStER_Input_File_oceanic_core_complex.m

max_snapshot = 140 ;
fault_zone_threshold = MAT(2).ecrit/1.2 ; % can be tuned to restrain fault zone
radius_array = zeros((max_snapshot/10 -1),1); % array for radius calculations
dip_array = zeros((max_snapshot/10 -1),1); % array to store fault dip through time
time_array = zeros((max_snapshot/10 -1),1);
dip_shallow_array = zeros((max_snapshot/10 -1),1);
dip_deep_array = zeros((max_snapshot/10 -1),1);
fault_length_array = zeros((max_snapshot/10 -1),1);
fault_vertical_array = zeros((max_snapshot/10 -1),1);
ep_cell_array = cell(2,max_snapshot/10 -1);
xim_cell_array = cell(2,max_snapshot/10 -1);
depth_ep_check = 1.8*1e4;
depth_xim_check = 1.5*1e4;
fit_deg = 2;%degree for radius fitting
nbins = 30;% number of bins for radius fitting

for k = 10:10:max_snapshot  
    
    load(num2str(k))
    Xfault = xm(ep>fault_zone_threshold);
    Yfault = ym(ep>fault_zone_threshold);
    tol = 0;
    while isempty(Xfault)
        tol = tol +0.1;
        Xfault = xm(ep>(MAT(2).ecrit/(1.2 +tol)));
        Yfault = ym(ep>(MAT(2).ecrit/(1.2 +tol)));
    end
    %
    
    vect_bins= linspace(min(Xfault),max(Xfault),nbins);
    Y_bins = zeros(1,nbins);

    for i = 2:1:nbins-1
        Y_bins(i) = mean(Yfault(Xfault> vect_bins(i-1) & Xfault< vect_bins(i+1)));
    end

    vect_bins = vect_bins(Y_bins>0);
    Y_bins = Y_bins(Y_bins>0);
    
    
    %polyfit_curve_fault_2 = polyfit(Xfault,Yfault, 1);
    polyfit_curve_fault_1 = polyfit(vect_bins,Y_bins, fit_deg);
    polyfit_curve_dip = polyfit(Xfault,Yfault,2);
    polyfit_curve_dip_first_order = polyfit(Xfault,Yfault,1);
    %x1 = Xfault(1):1:Xfault(end);
    %y1 = polyval(polyfit_curve_fault_1,x1);
    x1 = vect_bins(1):1:vect_bins(end);
    y1 = polyval(polyfit_curve_fault_1,x1);
    %
    fault_length_array(k/10) = sqrt((max(Xfault) - min(Xfault))^2 + (max(Yfault) - min(Yfault))^2);
    %
    fault_vertical_array(k/10) = max(Yfault)- min(Yfault);
    %
    x_shallow = min(Xfault);
    x_deep = max(Xfault);
    dip_shallow_array(k/10) = atand(x_shallow*polyfit_curve_dip(1)*2 +polyfit_curve_dip(2));
    dip_deep_array(k/10) = atand(x_deep*polyfit_curve_dip(1)*2 +polyfit_curve_dip(2));
    %y2 = polyval(
    %plot(Xfault,Yfault,'o',x1,y1,'--')
    dip_array(k/10) = atand(polyfit_curve_dip_first_order(1));
    grad_first = gradient(y1(:)) ./ gradient(x1(:));
    curvature_array = gradient(grad_first(:)) ./ gradient(x1(:));
    radius_array_time = 1./ (abs(curvature_array));
    radius_array(k/10) = mean(radius_array_time);
    time_array(k/10) = time;
    %%%% ep map through time
    ep_cell_array{1,k/10} = ep(ep> 0.1*MAT(2).ecrit & ep<MAT(2).ecrit/1.2 & ym < (depth_ep_check + 100) &ym > (depth_ep_check - 100));
    ep_cell_array{2,k/10} = xm(ep> 0.1*MAT(2).ecrit & ep<MAT(2).ecrit/1.2 & ym < (depth_ep_check + 100) &ym > (depth_ep_check - 100));
    %%%% xim map though time
    xim_cell_array{1,k/10} = xim(xim>0.1 & ym < (depth_xim_check + 100) &ym > (depth_xim_check - 100));
    xim_cell_array{2,k/10} = xm(xim>0.1 & ym < (depth_xim_check + 100) &ym > (depth_xim_check - 100));
    
end

figure
plot(time_array./(365*24*3600), radius_array,'*k')
xlabel('Time ( Myr)')
ylabel('Fault curvature radius (m)')
title('Time evolution of fault curvature radius')

figure
%plot(time_array./(365*24*3600), dip_array,'*r') average dip, not
%meaningfull
hold on
plot(time_array./(365*24*3600), dip_shallow_array,'*g')
plot(time_array./(365*24*3600), dip_deep_array,'*k')
xlabel('Time ( Myr)')
ylabel('Fault dip (°)')
title('Time evolution of fault dip')

figure
plot(time_array./(365*24*3600), fault_length_array,'*r')
xlabel('Time ( Myr)')
ylabel('Fault length (m)')
title('Time evolution of fault length')

figure
plot(time_array./(365*24*3600), fault_vertical_array,'*b')
xlabel('Time ( Yr)')
ylabel('Fault vertical extension (m)')
title('Time evolution of fault vertical extension, boosted diffusion below 400°C & 10km')

figure
hold on
%pointsize = 10;
for i= 1:(max_snapshot/10 -1)
    time_vect = time_array(i)/(365*24*3600*1e6)*ones(length(ep_cell_array{2,i}),1);
    scatter(ep_cell_array{2,i},time_vect,5,ep_cell_array{1,i}/MAT(2).ecrit)
end
xlabel('Distance from left border (m)')
ylabel('Time elapsed (Myr)')
title('Time localisation of ep without brittle failure @8 km under seafloor')


figure
hold on
%pointsize = 10;
for i= 1:(max_snapshot/10 -1)
    time_vect = time_array(i)/(365*24*3600*1e6)*ones(length(xim_cell_array{2,i}),1);
    scatter(xim_cell_array{2,i},time_vect,5,xim_cell_array{1,i})
end
xlabel('Distance from left border (m)')
ylabel('Time elapsed (Myr)')
title('Time localisation of serpentinisation @5 km under seafloor')

total_divergence = time * BC.right(3)*2;


%Script for elevation extraction from GeoMap App profiles and hand picking
%for curve fitting

%Profile loading
load 140
%plot(topo_x,topo_y)

%Break-away height calculation
top_break_away = 8274;
bottom_break_away = 9896 ;
h_break_away = bottom_break_away - top_break_away;


%plot(dist_13N20,elev_13N20)
begin_fit_point_13N20 = 66914; % manually picked on curve
end_fit_point_13N20 = 75926;%here end of first core complex 

OCC_length = end_fit_point_13N20 - begin_fit_point_13N20;

fit_dist_13N20 = topo_x(find(topo_x>begin_fit_point_13N20 & topo_x<end_fit_point_13N20)); %fitting topography
fit_elev_13N20 = topo_y(topo_x>begin_fit_point_13N20 & topo_x<end_fit_point_13N20);
%plot(fit_dist_13N20,fit_elev_13N20)

%fitting now curvature radius
fit_deg = 2; %degreee for polynomial fitting
polyfit_curve_13N20 = polyfit(fit_dist_13N20,fit_elev_13N20, fit_deg);
x1 = fit_dist_13N20(1):1:fit_dist_13N20(end);
y1 = polyval(polyfit_curve_13N20,x1);
figure
plot(fit_dist_13N20,fit_elev_13N20,'o',x1,y1,'--')

grad_first_13N20 = gradient(y1(:)) ./ gradient(x1(:));
curvature_array_13N20 = gradient(grad_first_13N20(:)) ./ gradient(x1(:));
radius_array_13N20 = 1./ (abs(curvature_array_13N20));
radius_average_13N20 = mean(radius_array_13N20);


% saving for systematic analysis

save('dilatancy_ecrit_01_mu_03_cmax_20E6_Hfluid_5km_serp_linear_045_data','radius_array','dip_shallow_array','dip_deep_array','time_array','fault_length_array','fault_vertical_array','radius_average_13N20','OCC_length','h_break_away','total_divergence') 
