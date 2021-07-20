% SiStER_Input_File


% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
Nt=800; % max number of time iterations
dt_out=10; % output files every "dt_out" iterations


% DOMAIN SIZE AND GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsize=150e3;
ysize=50e3;
% gridding- from 0 to GRID.x(1), grid size is GRID.dx(1)
% from GRID.x(1) to GRID.x(2), grid size is GRID.dx(1) etc...
% same for y
GRID.dx(1)=2000;
GRID.x(1)=60e3;
GRID.dx(2)=400;
GRID.x(2)=90e3;
GRID.dx(3)=2000;
GRID.dy(1)=2000;
GRID.y(1)=8.5e3;
GRID.dy(2)=400;
GRID.y(2)=22e3;
GRID.dy(3)=2000;


% LAGRANGIAN MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=8; % number of markers in the smallest quadrant
Mquad_crit=4; % minimum number of markers allowed in smallest quadrant (for reseeding)

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nphase=2; % number of phases

% phase 1
GEOM(1).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(1).top=0;
GEOM(1).bot=10e3;

% phase 2
GEOM(2).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(2).top=10e3;
GEOM(2).bot=50e3;



% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creep laws of the form: pre^(-1/n)*epsII^((1-n)/n)*exp(E/(nRT))
% harmonically averaging diffusion creep, dislocation creep 
% (and plastic creep to simulate brittle failure)

% phase 1
MAT(1).phase=1;
% density parameters
MAT(1).rho0=1000;
MAT(1).alpha=0;
% thermal parameters
MAT(1).k=3*50;
MAT(1).cp=1000;
% elasticity 
MAT(1).G=1e18;
% diffusion creep parameters
MAT(1).pre_diff=.5/1e18;
MAT(1).Ediff=0;
MAT(1).ndiff=1;
% dislocation creep parameters
MAT(1).pre_disc=.5/1e18;
MAT(1).Edisc=0;
MAT(1).ndisc=1;
% plasticity
MAT(1).mu=0.6;
MAT(1).mu_serp = 0.45;
MAT(1).mumin=0.3;
MAT(1).psi = 15;%dilantacy
MAT(1).Cmax=40e6;
MAT(1).Cmin=0.01e6;
MAT(1).ecrit=0.1;


% phase 2
MAT(2).phase=2;
% density parameters
MAT(2).rho0=3300;
MAT(2).alpha=0;
% thermal parameters
MAT(2).k=3;
MAT(2).cp=1000;
% elasticity 
MAT(2).G=30e9;
% diffusion creep parameters
MAT(2).pre_diff=.5/1e40;
MAT(2).Ediff=0;
MAT(2).ndiff=1;
% dislocation creep parameters
MAT(2).pre_disc=(2^3.5)*3.77e-14;%hirth and kohlstedt 2003
MAT(2).Edisc=520e3;
MAT(2).ndisc=3.5;
% plasticity
MAT(2).mu=0.6;
MAT(2).mu_serp = 0.45;
MAT(2).mumin=0.3;
MAT(2).psi = 15;%dilatancy
MAT(2).Cmax=20e6;% cohesion corrected, initially from Lavier
MAT(2).Cmin=0.01e6;
MAT(2).ecrit=0.1;




% ADDITIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.YNElast=1; % elasticity on (1) or off (0)
PARAMS.YNPlas=1; % plasticity on (1) or off (0)
PARAMS.tau_heal=1e12; % healing time for plasticity (s)
PARAMS.gx=0; % gravity along x
PARAMS.gy=9.8; % gravity along y
PARAMS.fracCFL=0.5; % distance by which a marker is allowed to move over a time step, expressed as a fraction of the smallest cell size
PARAMS.R=8.314; % gas constant
PARAMS.etamax=1e25; % maximum viscosity
PARAMS.etamin=1e18; % minimum viscosity
PARAMS.Tsolve=1; % yes (1) or no (0) solve for temperature
PARAMS.init_angle = 60; %initializing fault zone by weakening the model with specified angle, in degrees , down side, (0) for no weakening
PARAMS.center_depth = GEOM(1).bot + 5e3 ; %implanting fault point where it crosses simulation center
PARAMS.band_width = 2e3; % Weakness band width simulating fault

%SERPENTINIZATION CREATION PARAMETERS
% J.A Olive & B. Malvoisin & A. Demont, 3/2021, from Malvoisin et al. 2012
% & Malvoisin et al. 2021
%Tectonic parameters
PARAMS.Xseq = 0.12; %equilibrium state for serpentinization
PARAMS.beta_sr = 1.7; % volume created in tectonic expansion eq.
%Cinetic parameters
PARAMS.alpha = 95.457769814481940;
PARAMS.lambda =2;
PARAMS.Ea = 3.259163498830051e+04;%activation energy
PARAMS.k = 1.455589572232530e-05;%pre exponential factor
PARAMS.Tref = 6.410462984999999e+02;% IN KELVIN FOR CALCULATIONS

%FLUID CIRCULATION
PARAMS.fwidth = xsize;% region from center where fluid percolate, usually whole bos is concerned, xsize/2 minimum for whole box size
PARAMS.hfluids = 5e3;% depth for fluid action
PARAMS.kboost_up = 5;% boost for upper part of T in fluid circulation
PARAMS.kboost_bot = 5;% boost for lower part of fluid circulation
PARAMS.Tboundary_up = 300;%top boost boundary, make it the same as bot for linear profile
PARAMS.Tboundary_bot = 400;%T fluid limit boundary
PARAMS.lambda_f = 0.3; %pore fluid ratio(hydrostatic)

% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)
PARAMS.a0=0;
PARAMS.a1=0;
PARAMS.a2=0;
PARAMS.a3=0;
PARAMS.amp=0; % amplitude of sinusoidal perturbation
PARAMS.lam=1; % wavelength of sinusoidal perturbation
PARAMS.ynTreset=1; % if ==1, reset T=T0 where im==1 (sticky layer)
PARAMS.T0=0;
% reference values for the constant diffusivity thermal solver
% (kappa = kref / (rhoref*cpref))
PARAMS.rhoref=MAT(2).rho0; 
PARAMS.kref=3;
PARAMS.cpref=1000;

% TOPOGRAPHY EVOLUTION (interface between rock and sticky air/water layer)
PARAMS.Ntopo_markers=1000; % number of markers in marker chain tracking topography
PARAMS.YNSurfaceProcesses=0; % surface processes (diffusion of topography) on or off
PARAMS.topo_kappa=1e-8; % diffusivity of topography (m^2/s)


% Solver iterations
PARAMS.Npicard_min=10; % minimum number of Picard iterations per time step
PARAMS.Npicard_max=100; % maximum number of Picard iterations per time step
PARAMS.conv_crit_ResL2=1e-9;
PARAMS.pitswitch=0; % number of Picard iterations at which the solver switches to quasi-Newton



% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.pressure_seafloor = 20e6;


% pressure
PARAMS.p0cell=PARAMS.pressure_seafloor-(MAT(1).rho0*PARAMS.gy*GEOM(1).bot); % pressure in the top-left corner of the domain (anchor point)
%PARAMS.p0cell = 0;

% flow

% boundary conditions
% entries in BC correspond to
% 1/ rollers? (1=yes, 0=no)
% 2/ type of velocity normal to boundary (0=constant)
% 3/ value of normal velocity 

BC.top=[1 3 1.0563e-10];
BC.bot=[1 0 -1.0563e-10];
BC.left=[1 0 -3.1688e-10/2];
BC.right=[1 0 3.1688e-10/2];%1 cm/yr divergence rate
PARAMS.BalanceStickyLayer=1; % if set to 1, the code will reset the inflow 
% / outflow BCs to balance the inflow / outflow of sticky layer material,
% and rock separately, based on the position of the sticky layer / air
% interface


% thermal 

% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 2=Neumann)
% 2/ value

BCtherm.top=[1 0];
BCtherm.bot=[1 1350];
BCtherm.left=[2 0];
BCtherm.right=[2 0];


%Temperature parameters setting a half space cooling model, Turcotte&
%Oxburgh 19XX
BCM.hscm_on = 0;% define if we use half space cooling model on bottom, 1 if yes, 0 for no
BCM.tsurface = BCtherm.top(2);
BCM.tmagma = BCtherm.bot(2);
BCM.hocean = GEOM(1).bot;
BCM.kappa_hscm = .6e-6;% kappa parameter
BCM.xsize = xsize;
BCM.u = BC.right(3);
BCM.size = 50e3;


%PARAMETERS FOR TRIANGLE T INIT
BCM.d_ax = 20e3;%controls the mid ridge position of isoT BCTherm.bot
BCM.d_far = 30e3;% controls borders of the model T settings
BCM.width = 70e3;

%PARAMETER FOR FLAT INIT
BCM.flat_gradient = 0.06;
