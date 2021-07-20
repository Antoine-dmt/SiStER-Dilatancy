function [xim]=SiStER_serp_rate(im,ep,MAT,Rm,dx,dy,vx,vy,PARAMS,xim,Tm)

% compute serpentinization rate, based on cinetic reaction and tectonic
% dilatance
% J.-A. Olive 3/2021


%xim = zeros(size(im));

%types = unique(im);
dt_m=SiStER_set_timestep(dx,dy,vx,vy,PARAMS);

xim_cinetic_new = (PARAMS.k*exp(-PARAMS.Ea./(8.314.*(Tm+273.15))).*((1-exp((((Tm+273.15)-PARAMS.Tref)*PARAMS.alpha)./(8.314.*(Tm+273.15))))).^PARAMS.lambda).*(1-xim)*dt_m +xim;
xim_cinetic_new(Tm + 273.15>PARAMS.Tref) = 0;
xim_tecto_new = (1 -(Rm/((1+ PARAMS.beta_sr)*PARAMS.Xseq).*dt_m)).*xim+(Rm/((1+ PARAMS.beta_sr)*PARAMS.Xseq).*dt_m);
%xim_cin_new = [];
xim(Rm >0 & xim <1) = min(xim_tecto_new(Rm >0& xim<1),xim_cinetic_new(Rm >0& xim<1));
xim(xim>=1) =1;

return

