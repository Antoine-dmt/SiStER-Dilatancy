 function [T] = SiStER_hscm_function(BCM,x,y)
 T = BCM.tsurface +(BCM.tmagma-BCM.tsurface)*erf((y-BCM.hocean)./(2*sqrt(BCM.kappa_hscm *(abs(x-BCM.xsize/2)/BCM.u))));
 
end
