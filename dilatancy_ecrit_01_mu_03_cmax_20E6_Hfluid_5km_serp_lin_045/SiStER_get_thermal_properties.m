function [km, cpm]=SiStER_get_thermal_properties(im,MAT,PARAMS,Tm,xm,xsize,topo_x,topo_y,ym)
% [km, cpm]=SiStER_get_thermal_properties(im,MAT)
% obtain thermal conductivity and heat capacity

km = zeros(size(im));
cpm = km;

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    km(imInd) = MAT(types(i)).k;
    cpm(imInd) = MAT(types(i)).cp;
end

%ytopo_marker = interp1(topo_x,topo_y,xm);
%depth_marker = ym-ytopo_marker;
[fcm] = SiStER_locate_fluids(xm,ym,Tm,topo_x,topo_y,PARAMS,xsize,im);
km(fcm ==1) = PARAMS.kboost_bot*km(fcm == 1);
%km(Tm<PARAMS.Tboundary_up & abs(xm-(xsize/2)) < PARAMS.fwidth & im>1 & depth_marker<PARAMS.hfluids) = PARAMS.kboost_up* km(Tm<PARAMS.Tboundary_up & abs(xm-(xsize/2))< PARAMS.fwidth & im>1 & depth_marker<PARAMS.hfluids);
%km(Tm<PARAMS.Tboundary_bot & abs(xm-(xsize/2)) < PARAMS.fwidth & im>1 & depth_marker<PARAMS.hfluids) = PARAMS.kboost_bot* km(Tm<PARAMS.Tboundary_bot & abs(xm-(xsize/2))< PARAMS.fwidth & im>1 & depth_marker<PARAMS.hfluids);
        
  