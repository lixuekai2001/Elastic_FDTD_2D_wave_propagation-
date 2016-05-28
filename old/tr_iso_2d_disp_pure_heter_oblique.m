%                             Transverse isotropic case - Pure with
%                             heterogeneous operators from Zahradnik, 1995
% O.Ovcharenko, R.Martin, N.Fuji

close all;  %close all extra windows
clc;  %clear console
clear all; %clear all variables

% total number of grid points in each direction of the grid
 mg=1;
 NX =100*mg;  %X
 NY =100*mg;  %Y
 
  time=0.02d0;
%  time=0.1d0;
  
 % time step in seconds
 DELTAT = 0.5d-4;	%[sec] 
%  DELTAT = 0.2d-3;	%[sec] 
 DELTAT = DELTAT/mg;
 % total number of time steps
 %  NSTEP = 2000;
 NSTEP = round(time/DELTAT);
 time_vec = [1:NSTEP]'*DELTAT;
 
YMAX=50.d0; %[m]
XMAX=50.d0; %[m]

% YMAX=800.d0; %[m]
% XMAX=800.d0; %[m]

XMIN=0.d0;
YMIN=0.d0;

% DELTAT = 1.d-4;	%[sec]
% DELTAT = 1.d-4;	%[sec]
% DELTAT = 50.d-8;	%[sec]

DELTAX=(XMAX-XMIN)/NX; %[m]
DELTAY=(YMAX-YMIN)/NY; %[m]
%--------------------------------------------------------------------------
%---------------------- FLAGS ---------------------------------------------
% total number of time steps
% NSTEP = 2000;

% display information on the screen from time to time
IT_DISPLAY = 25;

%Take instant snapshotfalse
SNAPSHOT=false;
snapshot_time=400:50:NSTEP; %on what steps

%Use explosive source or gaussian?
EXPLOSIVE_SOURCE=false;

%Show source position on the graph?
SHOW_SOURCE_POSITION=true;

%Pause a little bit each iteration
PAUSE_ON=false;
pause_time=0.1; %[sec]

% To show or don't show wavefield 
SAVE_VX_JPG =false; %doesn't work, because I didn't pay attention to it yet
SAVE_VY_JPG =true;

%Record video - corresponding SAVE_VX or VY must be turned on
%because video is being created by capturing of current frame
%Matlab 2012 + required, saves video to a current folder
MAKE_MOVIE_VX=false;
MAKE_MOVIE_VY=false;
% tagv='ex0';
% tagv=['phetf' num2str(NX)];
tagv='aph1000';
% flags to add PML layers to the edges of the grid
USE_PML_XMIN = false;
USE_PML_XMAX = false;
USE_PML_YMIN = false;
USE_PML_YMAX = false;

DISP_NORM=true;
VEL_NORM=false;

DATA_TO_BINARY_FILE=false;
tag='mz_';

RED_BLUE=false;
COLORBAR_ON=true;
FE_BOUNDARY=true;

SAVE_SEISMOGRAMS=false;
% seis_tag='mz2';
% seis_tag='ex0';
seis_tag=['het2Dgauss' num2str(NX)];
LEFT_R = true;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
eps=0.00000000001d0;

if FE_BOUNDARY
     fprintf('ON. Lithological boundary\n');
    %  cp_above_eb=1800.d0;
    %  cp_below_eb=3300.d0;
    %  rho_above_eb=2400.d0;
    %  rho_below_eb=3200.d0;
     cp_below_eb=1800.d0;
     cp_above_eb=3300.d0;
%      cp_below_eb=3300.d0;
%      cp_above_eb=1800.d0;
rho_below_eb=2400.d0;
%      rho_above_eb=3200.d0;
     rho_above_eb=2400.d0;
     fprintf('  cp_above=%.2f rho_above=%.2f\n', cp_above_eb, rho_above_eb);
     fprintf('  cp_below=%.2f rho_below=%.2f\n', cp_below_eb, rho_below_eb);
else
     fprintf('OFF. Lithological boundary\n');
    %  cp_above_eb=3300.d0;
    %  cp_below_eb=3300.d0;
    %  rho_above_eb=1800.d0;
    %  rho_below_eb=1800.d0;
     cp_above_eb=1800.d0;
     cp_below_eb=1800.d0;
     rho_above_eb=2400.d0;
     rho_below_eb=2400.d0;
     fprintf('  cp_above=%.2f rho_above=%.2f\n', cp_above_eb, rho_above_eb);
     fprintf('  cp_below=%.2f rho_below=%.2f\n', cp_below_eb, rho_below_eb);
end

% % zinc, from Komatitsch et al. (2000)
% c11 = 16.5d10;
% c13 = 5.d10;
% c33 = 6.2d10;
% c44 = 3.96d10;
% rho = 7100.d0;
% f0 = 170.d3;
% 
% % apatite, from Komatitsch et al. (2000)
% c11 = 16.7d10;
% c13 = 6.6d10;
% c33 = 14.d10;
% c44 = 6.63d10;
% rho = 3200.d0;
% f0 = 300.d3;

% % isotropic material a bit similar to apatite
% c11 = 16.7d10;
% c13 = c11/3.d0;
% c33 = c11;
% c44 = (c11-c13)/2.d0;  % = c11/3.d0
% density = 3200.d0;
% f0 = 300.d3;

% % model I from Becache, Fauqueux and Joly, which is stable
% scale_aniso = 1.d10;
% c11 = 4.d0 * scale_aniso;
% c13 = 3.8d0 * scale_aniso;
% c33 = 20.d0 * scale_aniso;
% c44 = 2.d0 * scale_aniso;
% density = 4000.d0;  % used to be 1.
% f0 = 200.d3;

% model II from Becache, Fauqueux and Joly, which is stable
% scale_aniso = 1.d10;
% c11 = 20.d0 * scale_aniso;
% c13 = 3.8d0 * scale_aniso;
% c33 = c11;
% c44 = 2.d0 * scale_aniso;
% density = 4000.d0;  % used to be 1.
% f0 = 200.d3;

% cp = max(sqrt(c33/density),sqrt(c11/density));
  
% True isotropic
density= rho_above_eb;
cp = cp_above_eb;	%[km/s]
cs = cp / 1.732d0;	%[km/s]
lambda =density*(cp*cp - 2.d0*cs*cs);
mu = density*cs*cs;

c11 = (lambda + 2.d0*mu);
c13 = lambda;
c33 = c11;
c44 = mu;

 f0 = 170.d0;%0.d3;
% f0=100;


%From Roland Martin code
% from Becache et al., INRIA report, equation 7 page 5 http://hal.inria.fr/docs/00/07/22/83/PDF/RR-4304.pdf
  if(c11*c33 - c13*c13 <= 0.d0)
      disp('problem in definition of orthotropic material');
      break
  end

% check intrinsic mathematical stability of PML model for an anisotropic material
% from E. B\'ecache, S. Fauqueux and P. Joly, Stability of Perfectly Matched Layers, group
% velocities and anisotropic waves, Journal of Computational Physics, 188(2), p. 399-433 (2003)
  aniso_stability_criterion = ((c13+c44)^2 - c11*(c33-c44)) * ((c13+c44)^2 + c44*(c33-c44));
  fprintf('PML anisotropy stability criterion from Becache et al. 2003 = %e\n', aniso_stability_criterion);
  if(aniso_stability_criterion > 0.d0 && (USE_PML_XMIN  ||  USE_PML_XMAX  ||  USE_PML_YMIN  ||  USE_PML_YMAX))
     fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 1');
     break
  end

  aniso2 = (c13 + 2*c44)^2 - c11*c33;
  fprintf('PML aniso2 stability criterion from Becache et al. 2003 = %e\n',aniso2);
  if(aniso2 > 0.d0 && (USE_PML_XMIN  ||  USE_PML_XMAX  ||  USE_PML_YMIN  ||  USE_PML_YMAX))
     fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 2');
     break
  end

  aniso3 = (c13 + c44)^2 - c11*c33 - c44^2;
  fprintf('PML aniso3 stability criterion from Becache et al. 2003 = %e\n',aniso3);
  if(aniso3 > 0.d0 && (USE_PML_XMIN  ||  USE_PML_XMAX  ||  USE_PML_YMIN  ||  USE_PML_YMAX))
     fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 3');
     break
  end
 
%------------------------------------------------------------------


%Check if it is possible to save video
if MAKE_MOVIE_VY && ~SAVE_VY_JPG
    disp('Error. It is necesary to have SAVE_VY_JPG=true.');
    MAKE_MOVIE_VY=false;
end
if MAKE_MOVIE_VX && ~SAVE_VX_JPG
    disp('Error. It is necesary to have SAVE_VX_JPGquasi_cp_max=true.');
    MAKE_MOVIE_VX=false;
end

%vectors for visualisation using imagesec
nx_vec=[0:NX]*DELTAX;	%[m]
ny_vec=[0:NY]*DELTAY;


% thickness of the PML layer in grid points
NPOINTS_PML = 10;

% P-velocity, S-velocity and density
% cp = 3300.d0;	%[km/s]
% cs = cp / 1.732d0;	%[km/s]
% density = 2800.d0;	%[kg/m3]

% parameters for the source
% f0 = 10.d0;
t0 = 1.20d0 / f0;
% factor = 1.d7;
factor = 1.d6;
% factor=100.d0;

% ISOURCE = NX - round(NX/3);
% JSOURCE = round(NY / 3) + 1;
% ISOURCE = round(NX / 2);
% JSOURCE = round(NY / 2);

% xsource = (ISOURCE - 1) * DELTAX;
% ysource = (JSOURCE - 1) * DELTAY;
% angle of source force clockwise with respect to vertical (Y) axis
% ANGLE_FORCE = 75.d0;
if LEFT_R
    ANGLE_FORCE = 135.d0;
else
    ANGLE_FORCE = 45.d0;
end
% value of PI
PI = 3.141592653589793238462643d0;

% conversion from degrees to radians
DEGREES_TO_RADIANS = PI / 180.d0;

% zero
ZERO = 0.d0;

% large value for maximum
HUGEVAL = 1.d+30;

% velocity threshold above which we consider that the code became unstable
STABILITY_THRESHOLD = 1.d+25;


% source
% xsource = 33.d0;
% ysource = 16.5d0;
% xsource=round(5*NX/8)*DELTAX;
% % ysource=round(NY/4)*DELTAY;
% xsource=round(24*NX/40)*DELTAX;
% ysource=round(22*NY/40)*DELTAY;
if LEFT_R
    %left \
    xsource=round(18*NX/40)*DELTAX;
    ysource=round(16*NY/40)*DELTAY;
else
    %right /
    xsource=round(22*NX/40)*DELTAX;
    ysource=round(16*NY/40)*DELTAY;
end
% xsource=500.d0;
% ysource=500.d0;
% xsource=round(NX/2)*DELTAX;
% ysource=round(NY/2)*DELTAY;
dist = HUGEVAL;
for j = 2:NY
for i = 2:NX
  distval = sqrt((DELTAX*double(i-1) - xsource)^2 + (DELTAY*double(j-1) - ysource)^2);
  if(distval < dist)
    dist = distval;
    ISOURCE = i;
    JSOURCE = j;
  end
end
end
ISOURCE = ISOURCE - 1;
JSOURCE = JSOURCE - 1;
xsource = ISOURCE*DELTAX;
ysource = JSOURCE*DELTAY;

% main arrays
%displacements over X and Y
ux=zeros(3,NX+1,NY+1);
uy=zeros(3,NX+1,NY+1);
% variables that save wavefield at two previous time steps

velx=zeros(NX+1,NY+1);
vely=zeros(NX+1,NY+1);

%elastic parameters
rho=zeros(NX+1,NY+1);

%   total_energy_kinetic=zeros(NSTEP);
%   total_energy_potential=zeros(NSTEP);

% power to compute d0 profile
NPOWER = 2.d0;

K_MAX_PML = 1.d0; % from Gedney page 8.11
ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0); % from Festa and Vilotte

% arrays for the memory variables
% could declare these arrays in PML only to save a lot of memory, but proof of concept only here
memory_dux_dxx=zeros(NX+1,NY+1);
memory_duy_dyy=zeros(NX+1,NY+1);
memory_dux_dxy=zeros(NX+1,NY+1);
memory_duy_dxy=zeros(NX+1,NY+1);

% 1D arrays for the damping profiles
d_x=zeros(NX+1,1);
K_x=zeros(NX+1,1);
alpha_x=zeros(NX+1,1);
a_x=zeros(NX+1,1);
b_x=zeros(NX+1,1);

d_y=zeros(NY+1,1);
K_y=zeros(NY+1,1);
alpha_y=zeros(NY+1,1);
a_y=zeros(NY+1,1);
b_y=zeros(NY+1,1);

 
%Initiate video object for vx
if MAKE_MOVIE_VX
    movie_name_vx=[tagv '_vx_' num2str(NX) '_' num2str(NY) '_' num2str(DELTAX) '_' num2str(f0) '.avi'];
    vidObj_vx=VideoWriter(movie_name_vx);
    open(vidObj_vx);
end

%Initiate video object for vy
if MAKE_MOVIE_VY
    movie_name_vy=[tagv '_vy_' num2str(NX) '_' num2str(NY) '_' num2str(DELTAX) '_' num2str(f0) '.avi'];
    vidObj_vy=VideoWriter(movie_name_vy);
    open(vidObj_vy);
end

 
 %----------------------------------------
 %--- program starts here ----------------
 %----------------------------------------

fprintf('2D elastic finite-difference code in displacement formulation with C-PML\n\n');
fprintf('NX = %d\n',NX);
fprintf('NY = %d\n\n',NY);
fprintf('size of the model along X = %.2f\n',NX*DELTAX);
fprintf('size of the model along Y = %.2f\n\n',NY*DELTAY);
fprintf('Total number of grid points = %.2f\n\n',NX * NY);

%--- define profile of absorption in PML region ---
% thickness of the PML layer in meters
  thickness_PML_x = NPOINTS_PML * DELTAX;
  thickness_PML_y = NPOINTS_PML * DELTAY;

% reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0;

% check that NPOWER is okaymarkers=zeros(nx+1,ny+1);
  if(NPOWER < 1)       
      disp('NPOWER must be greater than 1');
      break;
  end

% R. Courant et K. O. Friedrichs et H. Lewy (1928)
  Courant_number = cp * DELTAT * sqrt(1.d0/DELTAX^2.d0 + 1.d0/DELTAY^2.d0);
  fprintf('Courant number = %.4f\n',Courant_number); 
  if Courant_number > 1.d0 
      disp('Error. Time step is too large, simulation will be unstable.');
      break;
  end
  
  fprintf(' dx=%.2f\n dy=%.2f\n dt=%e\n\n',DELTAX,DELTAY,DELTAT);

% compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1.d0) * cp * log(Rcoef) / (2.d0 * thickness_PML_x);
  d0_y = - (NPOWER + 1.d0) * cp * log(Rcoef) / (2.d0 * thickness_PML_y);
 
  fprintf('d0_x = %.2f\n',d0_x);
  fprintf('d0_y = %.2f\n\n',d0_y);

  d_x(:) = ZERO;
  K_x(:) = 1.d0;
  alpha_x(:) = ZERO;
  a_x(:) = ZERO;

  d_y(:) = ZERO;
  K_y(:) = 1.d0;
  alpha_y(:) = ZERO;
  a_y(:) = ZERO;
  %break;
%--------------------------------------------------------------------------
% damping in the X direction

% origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x;
  xoriginright = (NX-1)*DELTAX - thickness_PML_x;

for i = 1:NX+1
    % abscissa of current grid point along the damping profile
    xval = DELTAX * double(i-1);
    %---------- left edge
    if(USE_PML_XMIN)
        % define damping profile at the grid points
        abscissa_in_PML = xoriginleft - xval;
        if(abscissa_in_PML >= ZERO)
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            d_x(i) = d0_x * abscissa_normalized^NPOWER;
            % this taken from Gedney page 8.2
            K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
            alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
        end
    end

%---------- right edge
   if(USE_PML_XMAX)
        % define damping profile at the grid points
        abscissa_in_PML = xval - xoriginright;
        if(abscissa_in_PML >= ZERO)
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            d_x(i) = d0_x * abscissa_normalized^NPOWER;
            % this taken from Gedney page 8.2
            K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
            alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
        end
   end

    % just in case, for -5 at the end
    if(alpha_x(i) < ZERO) 
        alpha_x(i) = ZERO;
    end

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT);
  
    % this to avoid division by zero outside the PML
    if(abs(d_x(i)) > 1.d-6) 
        a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)));
    end    
  end

%--------------------------------------------------------------------------
% damping in the Y direction

% origin of the PML layer (position of right edge minus thickness, in meters)
  yoriginbottom = thickness_PML_y;
  yorigintop = (NY-1)*DELTAY - thickness_PML_y;

  for j = 1:NY+1
    % abscissa of current grid point along the damping profile
    yval = DELTAY * double(j-1);
    %---------- bottom edge
    if(USE_PML_YMIN)
      % define damping profile at the grid points
      abscissa_in_PML = yoriginbottom - yval;
      if(abscissa_in_PML >= ZERO)
        abscissa_normalized = abscissa_in_PML / thickness_PML_y;
        d_y(j) = d0_y * abscissa_normalized^NPOWER;
        % this taken from Gedney page 8.2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
      end
    end

%---------- top edge
    if(USE_PML_YMAX)
      % define damping profile at the grid points
      abscissa_in_PML = yval - yorigintop;
      if(abscissa_in_PML >= ZERO)
        abscissa_normalized = abscissa_in_PML / thickness_PML_y;
        d_y(j) = d0_y * abscissa_normalized^NPOWER;
        % this taken from Gedney page 8.2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
      end
    end

    b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DELTAT);
    
    % this to avoid division by zero outside the PML
    if(abs(d_y(j)) > 1.d-6)
        a_y(j) = d_y(j) * (b_y(j) - 1.d0) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)));
    end  
  end

if SAVE_SEISMOGRAMS
    fprintf('ON. Save seismograms\n');
    NREC=41;
    fprintf('Set %d recievers:\n',NREC);
%     xdeb=25.d0;
%     xfin=25.d0;
%     xdeb=xsource;
%     xfin=xsource;
    xdeb=0.2d0*XMAX;
    xfin=0.8d0*XMAX;
    ydeb=0.8d0*YMAX;
    yfin=0.2d0*YMAX;
    fprintf('  x0=%.2f  x1=%.2f\n  y0=%.2f  y1=%.2f\n', xdeb,xfin,ydeb,yfin);
    % way=0.0:0.50:100.0;
    % way=15.0:0.50:50.0;     %Reciever line

    % for receivers
    ix_rec=zeros(NREC,1);
    iy_rec=zeros(NREC,1);
    xrec=zeros(NREC,1);
    yrec=zeros(NREC,1); 

    % for seismograms
    seisux=zeros(NSTEP,NREC);
    seisuy=zeros(NSTEP,NREC);

    % [sharedVals,idxsIntoA] = intersect(xrec,way);
    xspacerec = (xfin-xdeb) / double(NREC-1);
    yspacerec = (yfin-ydeb) / double(NREC-1);
    for irec=1:NREC
         xrec(irec) = xdeb + double(irec-1)*xspacerec;
         yrec(irec) = ydeb + double(irec-1)*yspacerec;
    end
    NREC=length(xrec);
    % find closest grid point for each receiver
    for irec=1:NREC
       dist = HUGEVAL;
       for j = 2:NY
        for i = 2:NX
          distval = sqrt((DELTAX*double(i-1) - xrec(irec))^2 + (DELTAY*double(j-1) - yrec(irec))^2);
          if(distval < dist)
            dist = distval;
            ix_rec(irec) = i;
            iy_rec(irec) = j;
          end
        end
       end
       fprintf('Reciever %d at x= %.2f y= %.2f\n',irec,xrec(irec), yrec(irec));
    %    fprintf('closest grid point found at distance %.2f in i = %d\n',dist,ix_rec(irec));
    end
    fprintf('Source position:  '); % print position of the source
    fprintf('x = %.2f  ',xsource);
    fprintf('y = %.2f\n',ysource);

    fprintf('%d files will be saved to %s',2*NREC,pwd)
    fprintf('\n ...OK\n');
end


%--------------------------------------------------------------------------
% initialize arrays
  ux(:,:,:) = ZERO;
  uy(:,:,:) = ZERO;
  
  velx(:,:)=ZERO;
  vely(:,:)=ZERO;

% PML
  memory_dux_dxx(:,:) = ZERO;
  memory_duy_dyy(:,:) = ZERO;
  memory_duy_dxy(:,:) = ZERO;
  memory_dux_dxy(:,:) = ZERO;

% initialize seismograms
  seisux(:,:) = ZERO;
  seisuy(:,:) = ZERO;

% % initialize total energy
%   total_energy_kinetic(:) = ZERO;
%   total_energy_potential(:) = ZERO;

    
  %Set red-blue colormap for images
  CMAP=zeros(256,3);
  c1=[0 0 1]; %blue
  c2=[1 1 1]; %white
  c3=[1 0 0]; %red
  for nc=1:128
	  f=(nc-1)/128;
	  c=(1-sqrt(f))*c1+sqrt(f)*c2;
	  CMAP(nc,:)=c;
	  c=(1-f^2)*c2+f^2*c3;
	  CMAP(128+nc,:)=c;
	end
  if RED_BLUE
      colormap(CMAP);
  end
  
%grid point coordinates in physical domain
gr_x=zeros(NX+1,NY+1);
gr_y=zeros(NX+1,NY+1);
%calculate cartesian grid points
for i=1:NX+1
    for j=1:NY+1
        gr_x(i,j)=(i-1)*DELTAX;
        gr_y(i,j)=(j-1)*DELTAY;    
    end    
end
  
  
fprintf('Find closest grid nodes')  
curvature=0.2;
% curvature=0.0000001;
% xdscr=[0:NX]*DELTAX;
% ymid=((YMAX+YMIN)/2.d0)/3.d0;
% ydscr=-curvature*YMAX*sin(1.25*PI*xdscr/max(xdscr)+0.25*PI);

if LEFT_R
    ymid=YMAX;
    xdscr=linspace(XMIN,XMAX,(40*NX)+1);
    ydscr=-YMAX*xdscr/max(xdscr);
    ydscr=ymid+ydscr;
else
    ymid=0;
    xdscr=linspace(XMIN,XMAX,(40*NX)+1);
    ydscr=YMAX*xdscr/max(xdscr);
    ydscr=ymid+ydscr;
end

% ydscr=abs(min(ydscr))+ydscr+YMAX/4;
%calculate involved grid points, descritized coordinates of curve,
%normal vectors, coordinates of middles of the descritized samples.
%All the output variables are vectors
% [markers, xt_dis, yt_dis, nvecx, nvecy, xmn, ymn] = func_p_find_closest_grid_nodes(NX,NY,1,gr_x,gr_y ,xdscr, ydscr);
% clearvars xt_dis nvecx nvecy yt_dis xmn ymn;  
fprintf('...OK\n');

%------------------------------------------------------------------------
C=zeros(NX+1,NY+1,4);
        
% compute the Lame parameters and density  
% Create Cijkl matrix

nice_matrix=zeros(NX+1,NY+1);
densitya = rho_above_eb;
cpa = cp_above_eb;	%[km/s]
csa = cpa / 1.732d0;	%[km/s]
lambdaa =densitya*(cpa*cpa - 2.d0*csa*csa);
mua = densitya*csa*csa;
densityb = rho_below_eb;
cpb = cp_below_eb;	%[km/s]
csb = cpb / 1.732d0;	%[km/s]
lambdab =densityb*(cpb*cpb - 2.d0*csb*csb);
mub = densityb*csb*csb;
topo_szx=length(xdscr);
tgrx=round(topo_szx/NX);

f0=600.d0;
c11a = (lambdaa + 2.d0*mua);
c13a = lambdaa;
c33a = c11a;
c44a = mua;

c11b = (lambdab + 2.d0*mub);
c13b = lambdab;
c33b = c11b;
c44b = mub;

%  scale_aniso = 1.d10;
% c11b = 4.d0 * scale_aniso;
% c13b = 3.8d0 * scale_aniso;
% c33b = 20.d0 * scale_aniso;
% c44b = 2.d0 * scale_aniso;
% rho_below_eb = 4000.d0;  % used to be 1.
% % f0 = 200.d3;
% % f0 = 450.d0;
% f0 = 10.d0;
% % scale_aniso = 1.d10;
% % c11b = 4.d0 * scale_aniso;
% % c13b = 3.8d0 * scale_aniso;
% % c33b = 20.d0 * scale_aniso;
% % c44b = 2.d0 * scale_aniso;
% % rho_below_eb = 4000.d0;  % used to be 1.
% % % f0 = 200.d3;
% % f0 = 170.d0;
% 
% % % model II from Becache, Fauqueux and Joly, which is stable
%  scale_aniso = 1.d10;
% %  c11a = 20.d0 * scale_aniso;
% %  c13a = 3.8d0 * scale_aniso;
% %  c33a = c11a;
% %  c44a = 2.d0 * scale_aniso;
% %  rho_above_eb = 4000.d0;  % used to be 1.
% 
% %  scale_aniso = 1.d10;
% c11a = 4.d0 * scale_aniso;
% c13a = 3.8d0 * scale_aniso;
% c33a = 20.d0 * scale_aniso;
% c44a = 2.d0 * scale_aniso;
% rho_above_eb = 4000.d0;  % used to be 1.
% % 
% f0 = 40.d0;
% % % model II from Becache, Fauqueux and Joly, which is stable
% %  scale_aniso = 1.d10;
%  c11b = 20.d0 * scale_aniso;
%  c13b = 3.8d0 * scale_aniso;
%  c33b = c11a;
%  c44b = 2.d0 * scale_aniso;
%  rho_below_eb = 4000.d0;  % used to be 1.


fprintf('\nCreate C 6D %d elements\n',NX*NY*4);
for i = 1:NX
    x_trial=(1+(i-1)*tgrx):(i*tgrx);
    for j = 1:NY
        y_trial=ny_vec(j);
        if y_trial>=ydscr(x_trial)
            C(i,j,:)=[c11a c13a c33a c44a];
            rho(i,j) = rho_above_eb;
            nice_matrix(i,j)=1.d0;
            if i==NX
                C(i+1,j,:)=[c11a c13a c33a c44a];
                rho(i+1,j) = rho_above_eb;
                nice_matrix(i+1,j)=1.d0;
            end
            if j==NY
                C(i,j+1,:)=[c11a c13a c33a c44a];
                rho(i,j+1) = rho_above_eb;
                nice_matrix(i,j+1)=1.d0;
            end
        else
            C(i,j,:)=[c11b c13b c33b c44b];
            rho(i,j) = rho_below_eb;
            nice_matrix(i,j)=0.d0;
            if i==NX
                C(i+1,j,:)=[c11b c13b c33b c44b];
                rho(i+1,j) = rho_below_eb;
                nice_matrix(i+1,j)=1.d0;
            end
            if j==NY
                C(i,j+1,:)=[c11b c13b c33b c44b];
                rho(i,j+1) = rho_below_eb;
                nice_matrix(i,j+1)=1.d0;
            end
        end
    end
end

% dlmwrite('cijtr', C);
% fprintf('C(i,j,4) saved to %s\n', pwd);
clearvars densitya cpa csa lambdaa mua densityb cpb csb lambdab mub;
clearvars x_trial y_trial topo_szx tgrx;
clearvars c11a c13a c33a c44a c11b c13b c33b c44b;
fprintf('C(i,j,4) of size: %s  ...OK\n',num2str(size(C)));

fprintf('\n');


clearvars densitya cpa csa lambdaa mua densityb cpb csb lambdab mub;
clearvars x_trial y_trial topo_szx tgrx;
clearvars c11a c13a c33a c44a c11b c13b c33b c44b;

dx = DELTAX; 
dy = DELTAY;
dx2 = DELTAX^2.d0;
dy2 = DELTAY^2.d0;
ddx = 2.d0*DELTAX;
ddy = 2.d0*DELTAY;
dxdy4 = 4.d0*DELTAX*DELTAY;

fprintf('Used memory: %.2f mb\n', monitor_memory_whos);
input('\nPress Enter to start time loop ...');
%---------------------------------
%---  beginning of the time loop -----
%---------------------------------
for it = 1:NSTEP
    tic;
    ux(3,:,:)=ZERO;
    uy(3,:,:)=ZERO;
    for i = 2:NX
        for j = 2:NY
            rhov=rho(i,j);
            
            ae1 = (C(i+1,j,1) + C(i,j,1)) / 2.d0; %c11v * value_dux_dxx
            aw1 = (C(i-1,j,1) + C(i,j,1)) / 2.d0; %c11v * value_dux_dxx
            
            as3 = (C(i,j-1,3) + C(i,j,3)) / 2.d0; %c33v * value_duy_dyy
            an3 = (C(i,j+1,3) + C(i,j,3)) / 2.d0; %c33v * value_duy_dyy
            
            as4 = (C(i,j-1,4) + C(i,j,4)) / 2.d0; %c44v * value_dux_dyy
            an4 = (C(i,j+1,4) + C(i,j,4)) / 2.d0; %c44v * value_dux_dyy
            
            ae4 = (C(i+1,j,4) + C(i,j,4)) / 2.d0; %c44v * value_duy_dxx
            aw4 = (C(i-1,j,4) + C(i,j,4)) / 2.d0; %c44v * value_duy_dxx
            
            ase1 = (0.5d0*(C(i+1,j,1) + C(i,j,1)) + 0.5d0*(C(i+1,j-1,1) + C(i,j-1,1))) / 2.d0;
            asw1 = (0.5d0*(C(i-1,j,1) + C(i,j,1)) + 0.5d0*(C(i-1,j-1,1) + C(i,j-1,1))) / 2.d0;
            anw1 = (0.5d0*(C(i-1,j+1,1) + C(i,j+1,1)) + 0.5d0*(C(i-1,j,1) + C(i,j,1))) / 2.d0;
            ane1 = (0.5d0*(C(i+1,j+1,1) + C(i,j+1,1)) + 0.5d0*(C(i+1,j,1) + C(i,j,1))) / 2.d0;
            
            ase2 = (0.5d0*(C(i+1,j,2) + C(i,j,2)) + 0.5d0*(C(i+1,j-1,2) + C(i,j-1,2))) / 2.d0;
            asw2 = (0.5d0*(C(i-1,j,2) + C(i,j,2)) + 0.5d0*(C(i-1,j-1,2) + C(i,j-1,2))) / 2.d0;
            anw2 = (0.5d0*(C(i-1,j+1,2) + C(i,j+1,2)) + 0.5d0*(C(i-1,j,2) + C(i,j,2))) / 2.d0;
            ane2 = (0.5d0*(C(i+1,j+1,2) + C(i,j+1,2)) + 0.5d0*(C(i+1,j,2) + C(i,j,2))) / 2.d0;
            
            ase3 = (0.5d0*(C(i+1,j,3) + C(i,j,3)) + 0.5d0*(C(i+1,j-1,3) + C(i,j-1,3))) / 2.d0;
            asw3 = (0.5d0*(C(i-1,j,3) + C(i,j,3)) + 0.5d0*(C(i-1,j-1,3) + C(i,j-1,3))) / 2.d0;
            anw3 = (0.5d0*(C(i-1,j+1,3) + C(i,j+1,3)) + 0.5d0*(C(i-1,j,3) + C(i,j,3))) / 2.d0;
            ane3 = (0.5d0*(C(i+1,j+1,3) + C(i,j+1,3)) + 0.5d0*(C(i+1,j,3) + C(i,j,3))) / 2.d0;
            
            ase4 = (0.5d0*(C(i+1,j,4) + C(i,j,4)) + 0.5d0*(C(i+1,j-1,4) + C(i,j-1,4))) / 2.d0;
            asw4 = (0.5d0*(C(i-1,j,4) + C(i,j,4)) + 0.5d0*(C(i-1,j-1,4) + C(i,j-1,4))) / 2.d0;
            anw4 = (0.5d0*(C(i-1,j+1,4) + C(i,j+1,4)) + 0.5d0*(C(i-1,j,4) + C(i,j,4))) / 2.d0;
            ane4 = (0.5d0*(C(i+1,j+1,4) + C(i,j+1,4)) + 0.5d0*(C(i+1,j,4) + C(i,j,4))) / 2.d0;
            
            % c11v * value_dux_dxx
            value_dux_dxx = (ae1*(ux(2,i+1,j)-ux(2,i,j)) - aw1*(ux(2,i,j)-ux(2,i-1,j))) / dx2;
            % c44v * value_duy_dxx
            value_duy_dxx = (ae4*(uy(2,i+1,j)-uy(2,i,j)) - aw4*(uy(2,i,j)-uy(2,i-1,j))) / dx2;
            % c33v * value_duy_dyy
            value_duy_dyy = (as3*(uy(2,i,j-1)-uy(2,i,j)) - an3*(uy(2,i,j)-uy(2,i,j+1))) / dy2;
            % c44v * value_dux_dyy
            value_dux_dyy = (as4*(ux(2,i,j-1)-ux(2,i,j)) - an4*(ux(2,i,j)-ux(2,i,j+1))) / dy2;
            
            % c13v * value_dux_dxy
            value_dux_dxy = (-ase2*(ux(2,i+1,j) + ux(2,i+1,j-1) - ux(2,i,j) - ux(2,i,j-1)) - ...
                              asw2*(ux(2,i,j) + ux(2,i,j-1) - ux(2,i-1,j) - ux(2,i-1,j-1)) + ...
                              ane2*(ux(2,i+1,j+1) + ux(2,i+1,j) - ux(2,i,j+1) - ux(2,i,j)) + ...
                              anw2*(ux(2,i,j+1) + ux(2,i,j) - ux(2,i-1,j+1) - ux(2,i-1,j))) / dxdy4;
            % c11v * value_duy_dxy            
            value_duy_dxy = (-ase4*(uy(2,i+1,j) + uy(2,i+1,j-1) - uy(2,i,j) - uy(2,i,j-1)) - ...
                              asw4*(uy(2,i,j) + uy(2,i,j-1) - uy(2,i-1,j) - uy(2,i-1,j-1)) + ...
                              ane4*(uy(2,i+1,j+1) + uy(2,i+1,j) - uy(2,i,j+1) - uy(2,i,j)) + ...
                              anw4*(uy(2,i,j+1) + uy(2,i,j) - uy(2,i-1,j+1) - uy(2,i-1,j))) / dxdy4;
            % c44v * value_dux_dyx             
            value_dux_dyx = (ase4*(ux(2,i,j) + ux(2,i+1,j) - ux(2,i,j-1) - ux(2,i+1,j-1)) + ...
                             ane4*(ux(2,i,j+1) + ux(2,i+1,j+1) - ux(2,i,j) - ux(2,i+1,j)) - ...
                             anw4*(ux(2,i-1,j+1) + ux(2,i,j+1) - ux(2,i-1,j) - ux(2,i,j)) - ...
                             asw4*(ux(2,i-1,j) + ux(2,i,j) - ux(2,i-1,j-1) - ux(2,i,j-1))) / dxdy4;
            % c13v * value_duy_dyx            
            value_duy_dyx = (ase2*(uy(2,i,j) + uy(2,i+1,j) - uy(2,i,j-1) - uy(2,i+1,j-1)) + ...
                             ane2*(uy(2,i,j+1) + uy(2,i+1,j+1) - uy(2,i,j) - uy(2,i+1,j)) - ...
                             anw2*(uy(2,i-1,j+1) + uy(2,i,j+1) - uy(2,i-1,j) - uy(2,i,j)) - ...
                             asw2*(uy(2,i-1,j) + uy(2,i,j) - uy(2,i-1,j-1) - uy(2,i,j-1))) / dxdy4;
            

            %               memory_dux_dxx(i,j) = b_x(i) * memory_dux_dxx(i,j) + a_x(i) * value_dux_dxx;
            %               memory_duy_dyy(i,j) = b_y(j) * memory_duy_dyy(i,j) + a_y(j) * value_duy_dyy;

            %               value_dux_dxx = value_dux_dxx / K_x(i) + memory_dux_dxx(i,j);
            %               value_duy_dyy = value_duy_dyy / K_y(j) + memory_duy_dyy(i,j);
            %--------------------------------------------------------------------------------------------------------------------
            
                dt2rho=(DELTAT^2.d0)/rhov;

                sigmas_ux= value_dux_dxx + value_duy_dyx + value_dux_dyy + value_duy_dxy;
                sigmas_uy= value_dux_dyx + value_duy_dxx + value_dux_dxy + value_duy_dyy;

                ux(3,i,j) = 2.d0 * ux(2,i,j) - ux(1,i,j) + sigmas_ux * dt2rho;
                uy(3,i,j) = 2.d0 * uy(2,i,j) - uy(1,i,j) + sigmas_uy * dt2rho;
                                                               
%             if VEL_NORM
%                 velx(i,j)=(ux(3,i,j)-ux(1,i,j))/(2.d0*DELTAT);
%                 vely(i,j)=(uy(3,i,j)-uy(1,i,j))/(2.d0*DELTAT);
%             end
        end
    end
         
    % add the source (force vector located at a given grid point)
    a = pi*pi*f0*f0;
    t = double(it-1)*DELTAT;
    % Gaussian
     source_term = factor * exp(-a*(t-t0)^2);
     %source_term = factor * (t-t0);  
     % first derivative of a Gaussian
%       source_term =  -factor*2.d0*a*(t-t0)*exp(-a*(t-t0)^2);
    % Ricker source time function (second derivative of a Gaussian)
%      source_term = factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);

    force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
    force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
%     force_x =  source_term;
%     force_y =  source_term;
%       force_x=factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
%       force_y=factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);

    % define location of the source
    i = ISOURCE;
    j = JSOURCE;

    ux(3,i,j) = ux(3,i,j) + force_x * DELTAT / rhov;
    uy(3,i,j) = uy(3,i,j) + force_y * DELTAT / rhov;
    
%     if source_term>eps
%         ux(3,i,j) = force_x * DELTAT / rhov;
%         uy(3,i,j) = force_y * DELTAT / rhov;
%     end

    % Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
    ux(3,1,:) = ZERO;
    ux(3,NX+1,:) = ZERO;

    ux(3,:,1) = ZERO;
    ux(3,:,NY+1) = ZERO;

    uy(3,1,:) = ZERO;
    uy(3,NX+1,:) = ZERO;

    uy(3,:,1) = ZERO;
    uy(3,:,NY+1) = ZERO;

    % store seismograms
    if SAVE_SEISMOGRAMS
        for irec = 1:NREC
                seisux(it,irec) = ux(3,ix_rec(irec),iy_rec(irec));
                seisuy(it,irec) = uy(3,ix_rec(irec),iy_rec(irec));
        end   
    end
    
    %Set previous timesteps
    ux(1,:,:)=ux(2,:,:);
    ux(2,:,:)=ux(3,:,:);
    
    uy(1,:,:)=uy(2,:,:);
    uy(2,:,:)=uy(3,:,:);


    % compute total energy in the medium (without the PML layers)

    % compute kinetic energy first, defined as 1/2 rho ||v||^2
    % in principle we should use rho_half_x_half_y instead of rho for vy
    % in order to interpolate density at the right location in the staggered grid cell
    % but in a homogeneous medium we can safely ignore it
    
    %total_energy_kinetic(it) = 0.5d0 .*sum(rho((NPOINTS_PML+1):(NX-NPOINTS_PML),(NPOINTS_PML+1):(NY-NPOINTS_PML))*( ...
    %    vx((NPOINTS_PML+1):(NX-NPOINTS_PML),(NPOINTS_PML+1):(NY-NPOINTS_PML)).^2 +  ...
    %    vy((NPOINTS_PML+1):(NX-NPOINTS_PML),(NPOINTS_PML+1):(NY-NPOINTS_PML)).^2));

    % add potential energy, defined as 1/2 epsilon_ij sigma_ij
    % in principle we should interpolate the medium parameters at the right lo thencation
    % in the staggered grid cell but in a homogeneous medium we can safely ignore it
%     total_energy_potential(it) = ZERO;
%     for j = NPOINTS_PML+1: NY-NPOINTS_PML
%         for i = NPOINTS_PML+1: NX-NPOINTS_PML
%             epsilon_xx = ((lambda(i,j) + 2.d0*mu(i,j)) * sigmaxx(i,j) - lambda(i,j) * ...
%                 sigmayy(i,j)) / (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)));
%             epsilon_yy = ((lambda(i,j) + 2.d0*mu(i,j)) * sigmayy(i,j) - lambda(i,j) * ...
%                 sigmaxx(i,j)) / (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)));
%             epsilon_xy = sigmaxy(i,j) / (2.d0 * mu(i,j));
%             total_energy_potential(it) = total_energy_potential(it) + ...
%                 0.5d0 * (epsilon_xx * sigmaxx(i,j) .+ epsilon_yy * sigmayy(i,j) + 2.d0 * epsilon_xy * sigmaxy(i,j));
%         end
%     end
        % check stability of the code, exit if unstable
                % print maximum of norm of velocity
        velocnorm = max(sqrt(ux(3,:,:).^2 + uy(3,:,:).^2));
        if(velocnorm > STABILITY_THRESHOLD)
            break 
            disp('code became unstable and blew up');
        end   
    % output information
    if(mod(it,IT_DISPLAY) == 0 || it == 5)
        fprintf('Time step: %d\n',it)
        fprintf('Time: %.2f seconds\n',single((it-1)*DELTAT));
        toc;

        %fprintf('Max norm velocity vector V (m/s) = %.2f\n',velocnorm);
        %     fprintf('total energy = ',total_energy_kinetic(it) + total_energy_potential(it)
        %     print * 
    
        if(SAVE_VX_JPG || SAVE_VY_JPG)
            clf;	%clear current frame
            if DISP_NORM
                u=sqrt(ux(3,:,:).^2+uy(3,:,:).^2);
            elseif VEL_NORM
                u=sqrt((velx/max(max(velx))).^2+(vely/max(max(vely))).^2);
            else
                if SAVE_VX_JPG 
                    u=ux(3,:,:); 
                elseif SAVE_VY_JPG
                    u=uy(3,:,:);
                end
            end
            %velnorm(ISOURCE-1:ISOURCE+1,JSOURCE-1:JSOURCE+1)=ZERO;
            imagesc(nx_vec,ny_vec,squeeze(u(1,:,:))'); hold on;
            title(['Step = ',num2str(it),' Time: ',num2str(single((it-1)*DELTAT)),' sec']); 
            xlabel('m');
            ylabel('m');
            set(gca,'YDir','normal');
            if FE_BOUNDARY
                plot(xdscr,ydscr,'m');
            end
            if COLORBAR_ON
                colorbar();
            end
            drawnow;  hold on;
            if SHOW_SOURCE_POSITION
                scatter(xsource, ysource,'g','filled'); drawnow;
            end
           
            if SNAPSHOT
                if  nnz(snapshot_time==it)>0
                    snapshat = getframe(gcf);
                    imgg = frame2im(snapshat);
                    scrsht_name=['im' num2str(it) '.png'];
                    imwrite(imgg,scrsht_name);
                    fprintf('Screenshot %s saved to %s\n', scrsht_name, pwd);
                    clearvars scrsht_name imgg snapshat
                end  
            end
            
            if MAKE_MOVIE_VY
                F_y=getframe(gcf);  %-  capture figure or use gcf to get current figure. Or get current
                writeVideo(vidObj_vy,F_y);  %- add frame to the movie
                fprintf('Frame for %s captured\n',movie_name_vy);
            end
            
            if DATA_TO_BINARY_FILE
                filename=[tag 'u_' 'disp_t_' num2str(it) '.txt'];
                dlmwrite(filename, u);
                fprintf('Data file %s saved to %s\n',filename, pwd);
            end
        end
        fprintf('\n'); 
        if it==3150
            input('Next?');
        end
    end
    if PAUSE_ON
        pause(pause_time);
    end
end
  % end of time loop
  
  current_folder=pwd;	%current path
  if MAKE_MOVIE_VX
	  close(vidObj_vx);     %- close video file
      printf('Video %s saved in %s\n',movie_name_vx,current_folder);
  end
  
  if MAKE_MOVIE_VY
	  close(vidObj_vy);     %- close video file
      fprintf('Video %s saved in %s\n',movie_name_vy, current_folder);
  end
  
 if SAVE_SEISMOGRAMS
      for i=1:NREC
          filename=[seis_tag 'ux' '4x' num2str((ix_rec(i)-1)*dx,'%.2f') 'y' num2str((iy_rec(i)-1)*dy,'%.2f') '_' num2str(i) '.txt'];
          dlmwrite(filename, [time_vec, seisux(:,i)]);
          fprintf('ux. Seismogram for rec at %.2f %.2f saved as %s to %s\n', (ix_rec(i)-1)*dx, (iy_rec(i)-1)*dy, filename, pwd);
          filename=[seis_tag 'uy' '4y' num2str((iy_rec(i)-1)*dy,'%.2f') 'x' num2str((ix_rec(i)-1)*dx,'%.2f') '_' num2str(i) '.txt'];
          dlmwrite(filename, [time_vec, seisuy(:,i)]);
          fprintf('uy. Seismogram for rec at %.2f %.2f saved as %s to %s\n', (ix_rec(i)-1)*dx, (iy_rec(i)-1)*dy, filename, pwd);
      end
 end
 
  disp('End');
 
