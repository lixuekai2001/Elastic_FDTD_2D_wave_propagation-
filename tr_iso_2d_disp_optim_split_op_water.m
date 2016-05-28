%                             Transverse isotropic case
%           Operators splitted for dx2 dy2 dxdy

close all;  %close all extra windows
clc;        %clear console
clear all;  %clear all variables

% total number of grid points in each direction of the grid
 mg= 1;
 NX =100*mg;  %X
 NY =100*mg;  %Y
 
%  time=0.02d0;
  time=0.16d0;
 
 % time step in seconds
%  DELTAT = 0.5d-4;	%[sec] 
 DELTAT = 0.2d-3;	%[sec] 
 DELTAT = DELTAT/mg;
 % total number of time steps
 %  NSTEP = 2000;
 NSTEP = round(time/DELTAT);
 time_vec = [1:NSTEP]'*DELTAT;
 
% YMAX=50.d0; %[m]
% XMAX=50.d0; %[m]

YMAX=200.d0; %[m]
XMAX=200.d0; %[m]

XMIN=0.d0;
YMIN=0.d0;

% DELTAT = 1.d-4;	%[sec]
% DELTAT = 1.d-4;	%[sec]
% DELTAT = 50.d-9;	%[sec]

DELTAX=(XMAX-XMIN)/NX; %[m]
DELTAY=(YMAX-YMIN)/NY; %[m]
%--------------------------------------------------------------------------
%---------------------- FLAGS ---------------------------------------------
% total number of time steps
% NSTEP = 350;
% NSTEP=2000;
% display information on the screen from time to time
IT_DISPLAY = 25;

%Take instant snapshot
SNAPSHOT=false;
snapshot_time=50:25:NSTEP; %on what steps

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
% tagv='mz2triso';
tagv='mzm100';

% flags to add PML layers to the edges of the grid
USE_PML_XMIN = false;
USE_PML_XMAX = false;
USE_PML_YMIN = false;
USE_PML_YMAX = false;

DISP_NORM=true;    %show normal displacement
VEL_NORM=false;     %show normal velocity

DATA_TO_BINARY_FILE=false;  %save data to .txt files
tag='mz_';

SAVE_SEISMOGRAMS=false;
seis_tag=['mzcurvetriso2D' num2str(NX)];

RED_BLUE=false;      %use custom red-blue only colormap
COLORBAR_ON=true;   %show colorbar
FE_BOUNDARY=true;   %homogeneous or heterogeneous media
WATER = false;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
eps=0.00000000001d0;

if FE_BOUNDARY
     fprintf('ON. Lithological boundary\n');     
     % Water
     cp_above_eb=1500.d0;
     rho_above_eb=1000.d0;
     % Seabed
     cp_below_eb=3000.d0;
     rho_below_eb=2400.d0;

     fprintf('  cp_above=%.2f rho_above=%.2f\n', cp_above_eb, rho_above_eb);
     fprintf('  cp_below=%.2f rho_below=%.2f\n', cp_below_eb, rho_below_eb);
else
     fprintf('OFF. Lithological boundary\n');
%      cp_above_eb=1800.d0;
%      cp_below_eb=1800.d0;
%      rho_above_eb=2400.d0;
%      rho_below_eb=2400.d0;
     cp_above_eb=1500.d0;
     rho_above_eb=1000.d0;
     % Seabed
     cp_below_eb=1500.d0;
     rho_below_eb=1000.d0;
     fprintf('  cp_above=%.2f rho_above=%.2f\n', cp_above_eb, rho_above_eb);
     fprintf('  cp_below=%.2f rho_below=%.2f\n', cp_below_eb, rho_below_eb);
end

check_CFL(cp_below_eb, DELTAT, DELTAX, DELTAY);
% % zinc, from Komatitsch et al. (2000)
% c11 = 16.5d10;
% c13 = 5.d10;
% c33 = 6.2d10;
% c44 = 3.96d10;
% rho = 7100.d0;
% % f0 = 170.d3;
% f0 = 100;
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
% f0 = 450.d3;

% model II from Becache, Fauqueux and Joly, which is stable
% scale_aniso = 1.d10;
% c11 = 20.d0 * scale_aniso;
% c13 = 3.8d0 * scale_aniso;
% c33 = c11;
% c44 = 2.d0 * scale_aniso;
% density = 4000.d0;  % used to be 1.
% f0 = 200.d3;
% f0=170.d0;
% density= rho;
% cp = max(sqrt(c33/density),sqrt(c11/density));
  
% True isotropic
% density= rho_above_eb;
% cp = cp_above_eb;	%[km/s]
% cs = cp / 1.732d0;	%[km/s]
% lambda =density*(cp*cp - 2.d0*cs*cs);
% mu = density*cs*cs;

% c11 = (lambda + 2.d0*mu);
% c13 = lambda;
% c33 = c11;
% c44 = mu;
 
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


% P-velocity, S-velocity and density
% cp = 3300.d0;	%[km/s]
% cs = cp / 1.732d0;	%[km/s]
% density = 2800.d0;	%[kg/m3]

% parameters for the source
f0 = 50.d0;
t0 = 1.20d0 / f0;
factor = 1.d6;

% source
ANGLE_FORCE = 90.d0;

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


% xsource=round(NX/2)*DELTAX;
% ysource=round(NY/4)*DELTAY;

xsource=round(NX/4)*DELTAX;
ysource=round(0.95*NY)*DELTAY;

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
ISOURCE = ISOURCE -1;
JSOURCE = JSOURCE -1;
xsource=ISOURCE*DELTAX;
ysource=JSOURCE*DELTAY;

% main arrays
ux=zeros(3,NX+1,NY+1);
uy=zeros(3,NX+1,NY+1);

%elastic parameters
rho=zeros(NX+1,NY+1);

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
fprintf('NX = %d  ',NX);
fprintf('NY = %d  ',NY);
fprintf('%d in total\n',NX*NY);
fprintf(' dx = %f  dy=%f  dt=%e\n',DELTAX,DELTAY,DELTAT);
fprintf('Size of the model: %.2f m x ',NX*DELTAX);
fprintf('%.2f\n',NY*DELTAY);
fprintf('\n');

if SAVE_SEISMOGRAMS
    fprintf('ON. Save seismograms\n');
    NREC=41;
    fprintf('Set %d recievers:\n',NREC);
    % xdeb=XMIN;
    % xfin=XMAX;
%     xdeb=25.d0;
%     xfin=25.d0;
%     xdeb=xsource;
%     xfin=xsource;
    ydeb=0.2d0*YMAX;
    yfin=0.8d0*YMAX;
    xdeb=xsource;
    xfin=xsource;
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


%Reflection and transition coefficients
Refl_coef=(rho_below_eb*cp_below_eb-rho_above_eb*cp_above_eb)/(rho_below_eb*cp_below_eb+rho_above_eb*cp_above_eb);
Trans_coef=2.d0*rho_below_eb*cp_below_eb/(rho_below_eb*cp_below_eb+rho_above_eb*cp_above_eb);
if Refl_coef<ZERO
  tmps=', inverse polarity';
else
  tmps='';
end
fprintf('Below --> Above:\n');
fprintf('  R= %.2f - reflection%s\n  T= %.2f - transmition\n', Refl_coef, tmps ,Trans_coef);
clearvars  Refl_coef Trans_coef tmps;

% initialize arrays
  ux(:,:,:) = ZERO;
  uy(:,:,:) = ZERO;
  
  velx(:,:)=ZERO;
  vely(:,:)=ZERO;

% % PML
%   memory_dux_dxx(:,:) = ZERO;
%   memory_duy_dyy(:,:) = ZERO;
%   memory_duy_dxy(:,:) = ZERO;
%   memory_dux_dxy(:,:) = ZERO;

% initialize seismograms
  seisux(:,:) = ZERO;
  seisuy(:,:) = ZERO;

% % initialize total energy
%   total_energy_kinetic(:) = ZERO;
%   total_energy_potential(:) = ZERO;

fprintf('Set custom colormap');
 if RED_BLUE 
      CMAP = make_red_blue_colormap();
      colormap(CMAP);
     fprintf('...OK\n');   
  end
  
fprintf('Cartesian grid generation');
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
fprintf('...OK\n');
  
fprintf('Find closest grid nodes')  
% curvature=0.2;
curvature=0.0000001;
% xdscr=[0:NX]*DELTAX;
ymid=DELTAY/3.d0+3*(YMAX+YMIN)/5.d0;
xdscr=linspace(0,XMAX+0.1*DELTAX,(40*NX)+1);
ydscr=-curvature*YMAX*sin(1.25*PI*xdscr/max(xdscr)+0.25*PI);
% ydscr=YMAX*ydscr/4;
ydscr=ymid+ydscr;
% ydscr=abs(min(ydscr))+ydscr+YMAX/4;
%calculate involved grid points, descritized coordinates of curve,
%normal vectors, coordinates of middles of the descritized samples.
%All the output variables are vectors
[markers, xt_dis, yt_dis, nvecx, nvecy, xmn, ymn] = func_p_find_closest_grid_nodes(NX,NY,1,gr_x,gr_y ,xdscr, ydscr);
 clearvars xt_dis nvecx nvecy yt_dis xmn ymn;  
fprintf('...OK\n')

%------------------------------------------------------------------------
C=zeros(NX+1,NY+1,4);
        
% compute the Lame parameters and density  
% Create Cijkl matrix

nice_matrix=zeros(NX+1,NY+1);

densitya = rho_above_eb;
cpa = cp_above_eb;	%[km/s]
csa = cpa / 1.732d0;	%[km/s]

if WATER
     csa = 0.d0;
end

lambdaa =densitya*(cpa*cpa - 2.d0*csa*csa);
mua = densitya*csa*csa;

densityb = rho_below_eb;
cpb = cp_below_eb;	%[km/s]
csb = cpb / 1.732d0;	%[km/s]
lambdab =densityb*(cpb*cpb - 2.d0*csb*csb);
mub = densityb*csb*csb;

topo_szx=length(xdscr);
tgrx=round(topo_szx/NX);

c11a = (lambdaa + 2.d0*mua);
c13a = lambdaa;
c33a = c11a;
c44a = mua;

% c11a = lambdaa;
% c13a =lambdaa;
% c33a = lambdaa;
% c44a = lambdaa;


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
% f0 = 200.d3;
% f0 = 450.d0;
% f0 = 250.d0;
% scale_aniso = 1.d10;
% c11b = 4.d0 * scale_aniso;
% c13b = 3.8d0 * scale_aniso;
% c33b = 20.d0 * scale_aniso;
% c44b = 2.d0 * scale_aniso;
% rho_below_eb = 4000.d0;  % used to be 1.
% % f0 = 200.d3;
% f0 = 170.d0;

% % model II from Becache, Fauqueux and Joly, which is stable
%  scale_aniso = 1.d10;
% %  c11a = 20.d0 * scale_aniso;
% %  c13a = 3.8d0 * scale_aniso;
% %  c33a = c11a;
% %  c44a = 2.d0 * scale_aniso;
% %  rho_above_eb = 4000.d0;  % used to be 1.

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
                nice_matrix(i+1,j)=0.d0;
            end
            if j==NY
                C(i,j+1,:)=[c11b c13b c33b c44b];
                rho(i,j+1) = rho_below_eb;
                nice_matrix(i,j+1)=0.d0;
            end
        end
    end
end
clearvars densitya cpa csa lambdaa mua densityb cpb csb lambdab mub;
clearvars x_trial y_trial topo_szx tgrx;
clearvars c11a c13a c33a c44a c11b c13b c33b c44b;
fprintf('C(i,j,4) of size: %s  ...OK\n',num2str(size(C)));

fprintf('\n');

check_material_stability(C);

fprintf('Constructing coeff{i,j}');
% fprintf('Keep calm. It can take couple of minutes.\n');

arr_eta0x=zeros(NX,NY,9);
arr_eta1x=zeros(NX,NY,9);
arr_eta0y=zeros(NX,NY,9);
arr_eta1y=zeros(NX,NY,9);

coeffux=cell(NX,NY);
coeffuy=cell(NX,NY);
% coeffuxm=zeros(NX,NY,3,4);
% coeffuym=zeros(NX,NY,3,4);

dx = DELTAX; 
dy = DELTAY;
dx2 = DELTAX^2.d0;
dy2 = DELTAY^2.d0;
ddx = 2.d0*DELTAX;
ddy = 2.d0*DELTAY;
dxdy4 = 4.d0*DELTAX*DELTAY;

% one_over_2dx2 = 1.d0/(2.d0*DELTAX^2.d0);
% one_over_2dy2 = 1.d0/(2.d0*DELTAY^2.d0);
% one_over_2dxdy4 = 1.d0/(2.d0*dxdy4);

% tmp_coeff=create_coeff(dx,dy);
%ux(2,i+1,j+1); ux(2,i+1,j); ux(2,i+1,j-1); ux(2,i,j+1); ux(2,i,j); ux(2,i,j-1); ux(2,i-1,j+1); ux(2,i-1,j); ux(2,i-1,j-1)
% tmp_coeff=[0 0 0 0 1.d0 0 0 0 0;...         %u
%            0 1.d0 0 0 0 0 0 -1.d0 0;...     %ux
%            0 0 0 1.d0 0 -1.d0 0 0 0;...     %uy
%            0 1.d0 0 0 -2.d0 0 0 1.d0 0;...  %uxx
%            0 0 0 1.d0 -2.d0 1.d0 0 0 0;...  %uyy
%            1.d0 0 -1.d0 0 0 0 -1.d0 0 1.d0];%uxyassembly
%        
% denom_for_tmp_coeff=ones(size(tmp_coeff));  %u
% denom_for_tmp_coeff(2,:)=1.d0/ddx;          %ux
% denom_for_tmp_coeff(3,:)=1.d0/ddx;          %uy
% denom_for_tmp_coeff(4,:)=1.d0/dx2;          %uxx
% denom_for_tmp_coeff(5,:)=1.d0/dy2;          %uyy
% denom_for_tmp_coeff(6,:)=1.d0/dxdy4;        %uxy
% tmp_coeff=tmp_coeff.*denom_for_tmp_coeff;

tmp_dx2=[1.d0 -2.d0 1.d0]/dx2;
tmp_dy2=[1.d0 -2.d0 1.d0]/dy2;
tmp_dxdy=[1.d0 -1.d0 -1.d0 1.d0]/dxdy4;

tic;
for i=2:NX %over OX
    for j=2:NY %over OY

       % Construct eta0 and eta1 arrays for each marked point
        if markers(i,j)>0
%             pt0x=gr_x(i,j);
%             pt0y=gr_y(i,j);
%             ctr=0;
%             for ik=1:-1:-1 %over columns of points from right to left
%                 for jk=1:-1:-1  %over rows of points from top to bottom
%                         pt1x=gr_x(i+ik,j);
%                         pt1y=gr_y(i,j+jk);
%                         ctr=ctr+1;
%                         x_trial=linspace(pt0x,pt1x,20);
%                         y_trial=linspace(pt0y,pt1y,20);
%                         [xi,yi]=curveintersect(x_trial,y_trial,xdscr, ydscr);
%                         if ~isempty([xi,yi]) % check if there is an intersection
%                             if size(xi,1)*size(xi,2)>1  %get rid of multiple intersections
%                                 xi=xi(1);
%                                 yi=yi(1);
%                             end
%                             delta_x=abs(pt1x-pt0x);
%                             delta_y=abs(pt1y-pt0y);
%                             if delta_x<eps || delta_y<eps
%                                 if delta_x<eps %vertical line
%                                     arr_eta0y(i,j,ctr)=abs(yi-pt0y)/delta_y;
%                                     arr_eta1y(i,j,ctr)=1.d0-arr_eta0y(i,j,ctr);
%                                     arr_eta0x(i,j,ctr)=ZERO;
%                                     arr_eta1x(i,j,ctr)=ZERO;                                    
%                                 end
%                                 if delta_y<eps %hoizontal line
%                                     arr_eta0x(i,j,ctr)=abs(xi-pt0x)/delta_x;
%                                     arr_eta1x(i,j,ctr)=1.d0-arr_eta0x(i,j,ctr);
%                                     arr_eta0y(i,j,ctr)=ZERO;
%                                     arr_eta1y(i,j,ctr)=ZERO;
%                                 end
%                                 if delta_x<eps && delta_y<eps %if single point
%                                     arr_eta0y(i,j,ctr)=ZERO;
%                                     arr_eta1y(i,j,ctr)=ZERO;
%                                     arr_eta0x(i,j,ctr)=ZERO;
%                                     arr_eta1x(i,j,ctr)=ZERO;
%                                 end
%                             else %if evrything ok with deltas
%                                     arr_eta0y(i,j,ctr)=abs((yi-pt0y)/delta_y);
%                                     arr_eta1y(i,j,ctr)=1.d0-arr_eta0y(i,j,ctr);
%                                     arr_eta0x(i,j,ctr)=abs((xi-pt0x)/delta_x);
%                                     arr_eta1x(i,j,ctr)=1.d0-arr_eta0x(i,j,ctr);
%                             end
%                             
%                             %Define normal in point
%                             tmp=abs(xdscr-xi);
%                             [cvalue,idx]=min(tmp);
%                             if idx==1 
%                                 idx=2;
%                             end
%                             p1x=xdscr(idx-1); p2x=xi; p3x=xdscr(idx+1);
%                             p1y=ydscr(idx-1); p2y=yi; p3y=ydscr(idx+1);
%                             s12 = sqrt((p2x-p1x)^2+(p2y-p1y)^2);
%                             s23 = sqrt((p3x-p2x)^2+(p3y-p2y)^2);
%                             dxds = (s23^2*(p2x-p1x)+s12^2*(p3x-p2x))/(s12*s23*(s12+s23));
%                             dyds = (s23^2*(p2y-p1y)+s12^2*(p3y-p2y))/(s12*s23*(s12+s23));
%                             tvx=dxds;
%                             tvy=dyds;
%                             nvx=-dyds;
%                             nvy=dxds;
%                         else  % if there is no intersection
%                                 arr_eta1x(i,j,ctr)=1.d0*abs(ik);
%                                 arr_eta1y(i,j,ctr)=1.d0*abs(jk);
%                                 arr_eta0x(i,j,ctr)=abs(1.d0*ik)*abs((1.d0-arr_eta1x(i,j,ctr)));
%                                 arr_eta0y(i,j,ctr)=abs(1.d0*jk)*abs((1.d0-arr_eta1y(i,j,ctr)));
%                                 nvx=0.d0;
%                                 nvy=0.d0;
%                         end   % checking if there are intersections
%                         
%                         %Apply Mizutani operators
%                         eta0x = arr_eta0x(i,j,ctr);
%                         eta1x = arr_eta1x(i,j,ctr);
%                         eta0y = arr_eta0y(i,j,ctr);
%                         eta1y = arr_eta1y(i,j,ctr);                    
% 
%                         [A0,B0,A1,B1]=A0B0A1B1triso2(i,j,ik,jk,0,0,nvx,nvy,C,rho,dx,dy, ik*eta0x,-ik*eta1x, jk*eta0y, -jk*eta1y); % + 
%                         CJI=svdinv(B0*A0)*(B1*A1);
%                         coeffAux(ctr,:)=CJI(1,1:6);
%                         coeffAuy(ctr,:)=CJI(7,7:12);
%                 end      %end of jk loop
%             end  %%end of ik loop
%             
%             pre_dx2=svdinv([coeffAux(2,[1,2,4]); coeffAux(5,[1,2,4]); coeffAux(8,[1,2,4])]);
%             pre_dy2=svdinv([coeffAux(4,[1,3,5]); coeffAux(5,[1,3,5]); coeffAux(6,[1,3,5])]);
%             pre_dxdy=svdinv([coeffAux(1,:); coeffAux(3,:); coeffAux(7,:); coeffAux(9,:)]);
%             coeffux_dx2 = C(i,j,1)*pre_dx2(3,:);
%             coeffux_dy2 = C(i,j,4)*pre_dy2(3,:);
%             coeffux_dxdy= C(i,j,2)*pre_dxdy(6,:);
%       
%             pre_dx2=svdinv([coeffAuy(2,[1,2,4]); coeffAuy(5,[1,2,4]); coeffAuy(8,[1,2,4])]);
%             pre_dy2=svdinv([coeffAuy(4,[1,3,5]); coeffAuy(5,[1,3,5]); coeffAuy(6,[1,3,5])]);
%             pre_dxdy=svdinv([coeffAuy(1,:); coeffAuy(3,:); coeffAuy(7,:); coeffAuy(9,:)]);
%             coeffuy_dx2 = C(i,j,4)*pre_dx2(3,:);
%             coeffuy_dy2 = C(i,j,3)*pre_dy2(3,:);
%             coeffuy_dxdy= C(i,j,4)*pre_dxdy(6,:);
            [Aux, Auy] = construct_interface_operators(i,j, gr_x, gr_y, xdscr, ydscr, C, rho);
            
            Aux(1,:) = C(i,j,1)*Aux(1,:);
            Aux(2,:) = C(i,j,4)*Aux(2,:);
            Aux(3,:) = C(i,j,2)*Aux(3,:);
            Aux(4,:) = C(i,j,4)*Aux(4,:);
            coeffux{i,j}=Aux;
            
            Auy(1,:) = C(i,j,4)*Auy(1,:);
            Auy(2,:) = C(i,j,3)*Auy(2,:);
            Auy(3,:) = C(i,j,4)*Auy(3,:);
            Auy(4,:) = C(i,j,2)*Auy(4,:);
            coeffuy{i,j}=Auy;
            
%             coeffux{i,j}=[[coeffux_dx2 0]; [coeffux_dy2 0]; coeffux_dxdy];
%             coeffuy{i,j}=[[coeffuy_dx2 0]; [coeffuy_dy2 0]; coeffuy_dxdy];
%             coeffuxm(i,j,:,:) = [[coeffux_dx2 0]; [coeffux_dy2 0]; coeffux_dxdy];
%             coeffuym(i,j,:,:) = [[coeffuy_dx2 0]; [coeffuy_dy2 0]; coeffuy_dxdy];
            
        end  %end of if markers(i,j)
        
        %if any conditions were used - use conventional heterogeneous operator
        if isempty(coeffux{i,j}) && isempty(coeffuy{i,j})
            [Aux, Auy] = construct_Zahradnik_operators(i,j,C, DELTAX, DELTAY);
            coeffux{i,j} = Aux;
            coeffuy{i,j} = Auy;
        end    
    end % end of j loop
end  %end of i loop
fprintf('...OK\n')


fprintf('Check coeff{i,j} for explosions');

cux=[1.d0; 1.d0; 1.d0; C(i,j,1); C(i,j,4); C(i,j,2)];
cuy=[1.d0; 1.d0; 1.d0; C(i,j,4); C(i,j,3); C(i,j,1)];
mmAB=0;
for i=2:size(coeffux,1)-2
    for j=2:size(coeffux,2)-2
        A=max(max(coeffux{i,j}));
        B=max(max(coeffuy{i,j}));
        mAB=max(A,B);
        if mAB>mmAB
            mmAB=mAB;
        end
    end
end

if mmAB>100*max(cux) || mmAB>100*max(cuy)
    fprintf('...FAILED\n');
else
    fprintf('...OK\n');
    clearvars i j A B mAB mmAB;
end

toc;

%Clean up memory from temporary variables
clearvars coeffAux coeffAuy cux cuy tmp_coeff eta0x eta0y eta1x eta1y;
clearvars A0 B0 A1 B1 CJI ctr pt0x pt0y pt1x pt1y;
clearvars ik jk ii jj i j denom_for_tmp_coeff nvx nvy x_trial y_trial;
clearvars arr_eta0x arr_eta1x arr_eta0y arr_eta1y delta_x delta_y;
clearvars idx cvalue p1x p1y p2x p2y p3x p3y s12 s23 tmp tvx tvy
clearvars xi yi xdeb xfin ydeb yfin dxds dyds;
clearvars xoriginleft xoriginright nc;
clearvars cp cs cp_above_eb cp_below_eb density rho_below_eb rho_above_eb;
% clearvars tmp_dx2 tmp_dy2 tmp_dxdy coeffux_dx2 coeffux_dy2 coeffux_dxdy;
% clearvars coeffuy_dx2 coeffuy_dy2 coeffuy_dxdy;

fprintf('Used memory: %.2f mb\n', monitor_memory_whos);
input('\nPress Enter to start time loop ...');

%---------------------------------
%---  beginning of the time loop -----
%---------------------------------
for it = 1:NSTEP
    
    % calculate next step
    tic;
%     [ux, uy] = solver_mx_VTI_elastic(ux, uy, DELTAT, coeffux, coeffuy, rho, C);
%     [ux, uy] = solver_mx_acoustic(ux, uy, DELTAT, coeffux, coeffuy, rho);
    
%             if WATER
%                 sigmas_ux = value_dux_dxx + value_dux_dyy*C(i,j,1)/C(i,j,4);
%                 sigmas_uy = value_duy_dxx*C(i,j,1)/C(i,j,4) + value_duy_dyy;
%             end

    ux(3,:,:)=0.d0;
    uy(3,:,:)=0.d0;
    
    
    for i = 2:NX
        for j = 2:NY          
            rhov=rho(i,j);
            A_ux=coeffux{i,j};
            value_dux_dxx=A_ux(1,1:3)*[ux(2,i-1,j); ux(2,i,j); ux(2,i+1,j)];
            value_dux_dyy=A_ux(2,1:3)*[ux(2,i,j-1); ux(2,i,j); ux(2,i,j+1)];
            value_dux_dxy=A_ux(3,:)*[ux(2,i+1,j+1); ux(2,i+1,j-1); ux(2,i-1,j+1); ux(2,i-1,j-1)];
            value_dux_dyx=A_ux(4,:)*[ux(2,i+1,j+1); ux(2,i+1,j-1); ux(2,i-1,j+1); ux(2,i-1,j-1)];
            
            A_uy=coeffuy{i,j};
            value_duy_dxx=A_uy(1,1:3)*[uy(2,i-1,j); uy(2,i,j); uy(2,i+1,j)];
            value_duy_dyy=A_uy(2,1:3)*[uy(2,i,j-1); uy(2,i,j); uy(2,i,j+1)];
            value_duy_dxy=A_uy(3,:)*[uy(2,i+1,j+1); uy(2,i+1,j-1); uy(2,i-1,j+1); uy(2,i-1,j-1)];
            value_duy_dyx=A_uy(4,:)*[uy(2,i+1,j+1); uy(2,i+1,j-1); uy(2,i-1,j+1); uy(2,i-1,j-1)];
%             value_dux_dyx=value_dux_dxy*C(i,j,4)/C(i,j,2);
%             value_duy_dyx=value_duy_dxy*C(i,j,2)/C(i,j,4);

            %--------------------------------------------------------------------------------------------------------------------
            
            dt2rho=(DELTAT^2.d0)/rhov;
% 
%           sigmas_ux= c11v * value_dux_dxx + c13v * value_duy_dyx + c44v * value_dux_dyy + c44v * value_duy_dxy;
%           sigmas_uy= c44v * value_dux_dyx + c44v * value_duy_dxx + c13v * value_dux_dxy + c33v * value_duy_dyy;
            
            if WATER && nice_matrix(i,j)
                sigmas_ux= value_dux_dxx + value_dux_dyy;
                sigmas_uy= value_duy_dxx + value_duy_dyy;
            else
                sigmas_ux= value_dux_dxx + value_duy_dyx + value_dux_dyy + value_duy_dxy;
                sigmas_uy= value_dux_dyx + value_duy_dxx + value_dux_dxy + value_duy_dyy;
            end

            ux(3,i,j) = 2.d0 * ux(2,i,j) - ux(1,i,j) + sigmas_ux * dt2rho;
            uy(3,i,j) = 2.d0 * uy(2,i,j) - uy(1,i,j) + sigmas_uy * dt2rho;

        end
    end

    % add source
    t = double(it-1)*DELTAT;
    [force_x, force_y] = source_function(f0, t0, factor, ANGLE_FORCE, t);
    i = ISOURCE;
    j = JSOURCE;
    rhov = rho(i,j);
    ux(3,i,j) = ux(3,i,j) + force_x * DELTAT^2.d0 / rhov;
    uy(3,i,j) = uy(3,i,j) + force_y * DELTAT^2.d0 / rhov;
    
    
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


  
%     velocnorm = max(sqrt(ux(3,:,:).^2 + uy(3,:,:).^2));
%     if(velocnorm > STABILITY_THRESHOLD)
%         %break 
%         disp('code became unstable and blew up');
%     end   
    
    % output information
    if(mod(it,IT_DISPLAY) == 0 || it == 5)
        fprintf('Time step: %d\n',it)
        fprintf('Time: %.4f sec\n',single((it-1)*DELTAT));
        toc;
    
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
            title(['Step = ',num2str(it),' Time: ',sprintf('%.4f',single((it-1)*DELTAT)),' sec']); 
            xlabel('m');
            ylabel('m');
            set(gca,'YDir','normal');
            if FE_BOUNDARY
                plot(xdscr,ydscr,'m'); 
            end
            if COLORBAR_ON
                colorbar();
%                 set(gca, 'CLim', [minbar, maxbar]);
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
        if it==3100
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
 