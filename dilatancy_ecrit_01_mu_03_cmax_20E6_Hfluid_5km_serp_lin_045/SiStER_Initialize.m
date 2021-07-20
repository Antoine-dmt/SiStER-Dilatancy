% SiStER Initialize

PARAMS.Nphase = Nphase; % for convenience

% construct staggered grids
[X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid(xsize,ysize,GRID);

% initialize marker arrays and positions
[xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad);

% locate markers with respect to grid
[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);

% assign marker phases
[im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym);
% initialize marker plastic strain (to zero) and strain rate (to one)
ep=zeros(size(xm));
epNH=ep;
epsIIm=ones(size(xm));
xim = zeros(size(xm));%serpentinization reaction 
fcm = zeros(size(xm));%fluid circulation zone
% initialize marker stresses
sxxm=zeros(size(xm));
sxym=sxxm;

% initialize marker index (a unique number to identify and track each marker)
idm=1:length(xm);

% initialize temperature structure on nodes
%T=PARAMS.a0+PARAMS.a1*Y+PARAMS.a2*Y.^2+PARAMS.a3*Y.^3;
%T=T+PARAMS.amp*sin(2*pi*X/PARAMS.lam);
n = size(X,1);
p = size(X,2);


Xfake= zeros([n p]);
for i = 1:n
    for j = 1:p
        if abs(X(i,j)-xsize/2)> BCM.width
            Xfake(i,j) = xsize/2 + BCM.width;
        else
            Xfake(i,j) = X(i,j);
        end
    end
end

%New temperature profile solving heat equation, and introducing erf()
%function,model based on following paper : Bickert et al 2020
T = BCM.tsurface +(BCM.tmagma-BCM.tsurface)*erf((Y-BCM.hocean)./(BCM.d_ax+(abs(Xfake-xsize/2)/BCM.width)*(BCM.d_far-BCM.d_ax)));
%T = BCM.tsurface + BCM.flat_gradient*(Y-BCM.hocean);% flat T profile


if PARAMS.ynTreset==1 % reset T=T0 in top layer
    T(T<PARAMS.T0)=PARAMS.T0;
end
% pass initial nodal T to markers
[Tm]=SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);
Tm0=Tm;

% initialize nodal strain rate and other useful arrays
EXX=zeros(size(X));
EXY=zeros(size(X));
vx=zeros(size(X));
vy=zeros(size(X));
v=vx;
p=1e12*ones(size(EXX));  %initialize to be high so plasticity doesnt activate at t=1, pit=1;
etan_new=zeros(Ny,Nx);
%-------------------------------------------------------------------------
% initialize dt_m small to keep things elastic & no plasticity at t=1, G.Ito
%-------------------------------------------------------------------------
if (exist('dt_m','var')==0)
    dt_m=1e2;
end

% initialize marker chain to track base of layer 1 (sticky layer)
Ntopo=PARAMS.Ntopo_markers;
topo_x=linspace(0,xsize,Ntopo);
topo_y=GEOM(1).bot*ones(size(topo_x));
topo_marker_spacing=mean(diff(topo_x)); % initial mean spacing of topography markers


%initiating rm % WIP
Rm=zeros(size(xm));

%creation of a initial fault-like weakness zone,with an initialisation at
%epcrit, with a  gradient(down) specified in input file, fault will be
%located at the simulation's center


gradient = tand(PARAMS.init_angle);
origin_abs = xsize/2 -(PARAMS.center_depth/gradient);

%curve eq 

w = PARAMS.band_width;

iF = zeros(size(xm));

vect_dir = [xsize/2 - origin_abs , PARAMS.center_depth];

coeff_param_c = -vect_dir(2)*xsize/2 +vect_dir(1)*PARAMS.center_depth;

for k = 1:size(xm,2)% le but est de parcourir les marqueurs et de ne modifier un propriété que pour ceux qui sont à une distance W/2 de la droi
    if abs(vect_dir(2)*xm(k)-vect_dir(1)*ym(k)+coeff_param_c)/sqrt(vect_dir(2)^2 +vect_dir(1)^2) < w/2
        iF(k) = 1;
    end
end
ep(iF==1)=MAT(2).ecrit;