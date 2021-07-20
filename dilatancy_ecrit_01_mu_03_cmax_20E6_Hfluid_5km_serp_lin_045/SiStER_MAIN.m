% SiStER_MAIN.m
%
% Simple Stokes solver with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
% jaolive <at> ldeo.columbia.edu
% March 2011 - April 2017

close all

% INITIALIZATION

% Input File: loads parameter values, model geometry, boundary conditions
if exist('running_from_SiStER_RUN','var')==0
    clear 
    %InpFil = input('Input file ? ','s');
    InpFil = 'SiStER_Input_File_oceanic_core_complex';
end
run(InpFil)

% construct grid and initialize marker / node arrays
SiStER_Initialize

% BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0;

for t=1:Nt % time loop
    
    disp(['STARTING ITERATION: ' num2str(t) ' out of ' num2str(Nt)])
    
    % update time
    time=time+dt_m;
    %if t <= 10
    %    time = time+dt_m/10;
    %else
    %    time=time+dt_m;
    %end
    
    % Here we prepare nodal arrays to feed the Stokes solver 
    SiStER_material_props_on_nodes

    %%% SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE 
    SiStER_flow_solve
    
    % GET STRAIN RATE FROM CURRENT SOLUTION
    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
    
    % USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
    SiStER_update_marker_stresses;
    
    % BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
    if (PARAMS.YNPlas==1) 
        SiStER_update_ep;
    end
    
    %DILATANCY
    Rm = SiStER_interp_normal_nodes_to_markers(R_n,xc,yc,xm,ym,icn,jcn);
    
    %ALTERATION UPDATE IN SERPENTINIZED AREAS
    xim = SiStER_serp_rate(im,ep,MAT,Rm,dx,dy,vx,vy,PARAMS,xim,Tm);
    
    % OUTPUT VARIABLES OF INTEREST (prior to rotation & advection)
    if (mod(t,dt_out)==0 && dt_out>0) || t==1 || t==Nt % SAVING SELECTED OUTPUT
        disp('SAVING SELECTED VARIABLES TO OUTPUT FILE') 
        filename=num2str(t);
        [etam]=SiStER_interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn); % to visualize viscosity on markers
        %[Rm] = SiStER_interp_normal_nodes_to_markers(R_n,xc,yc,xm,ym,icn,jcn);
        %[xim] = SiStER_get_serp_rate(im,ep,MAT,Rm,dx,dy,vx,vy,PARAMS,xim);
        save(filename,'X','Y','vx','vy','p','time','xm','ym','etam','rhom','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','icn','jcn','qd','topo_x','topo_y','Rm','xim','fcm')
    end
    
    % SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS);

    % ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
    if (PARAMS.YNElast==1) 
        SiStER_rotate_stresses;
    end
    
    % EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
    if PARAMS.Tsolve==1
        SiStER_thermal_update;
    end

    % MARKER ADVECTION, REMOVAL, AND ADDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SiStER_move_remove_and_reseed_markers;
    % advect markers in current flow field
    % remove markers if necessary
    % add markers if necessary
    SiStER_update_topography_markers
    % here we do the same for the marker chain that keeps track of topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('---------------')
    disp(['END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs.'])
    disp('--------------------------------')
    disp('--------------------------------')
    

end

disp('FIN')

    