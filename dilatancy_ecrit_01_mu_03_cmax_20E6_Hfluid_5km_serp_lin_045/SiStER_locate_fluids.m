function [fcm] = SiStER_locate_fluids(xm,ym,Tm,topo_x,topo_y,PARAMS,xsize,im)

ytopo_marker = interp1(topo_x,topo_y,xm);
depth_marker = ym-ytopo_marker;
fcm = zeros(size(xm));
fcm(Tm<=PARAMS.Tboundary_bot & abs(xm-(xsize/2)) <= PARAMS.fwidth & im>1 & depth_marker<=PARAMS.hfluids) = 1;