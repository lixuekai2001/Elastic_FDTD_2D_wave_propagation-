%                             Transverse isotropic case - Pure
% conventional operators are given explicitly. I need it to compare with Roland's solution and understand
%if this approach works correctly

close all;  %close all extra windows
clc;  %clear console
clear all; %clear all variables

 mg=8;
 NX =100*mg + 1;  %X
 NY =100*mg + 1;  %Y

%  NX = 201;
%  NY = 201;
 
 %source parameters
 f0 = 16;
%  t0 = 1.20d0 / f0;
 factor = 1.d10;
 ANGLE_FORCE = 0.d0;
 
 t0 = 0.15d0;
%  DELTAT = 0.8d-3;
%  DELTAT = DELTAT/mg;
 DELTAT = 0.4d-3;
 time = 0.6d0;
%  NSTEP = 1500;
 NSTEP = round(time/DELTAT);
 time_vec = -t0 + [1:NSTEP]'*DELTAT;
 
YMAX=2000.d0; %[m]
XMAX=2000.d0; %[m]

XMIN=0.d0;
YMIN=0.d0;

DELTAX=(XMAX-XMIN)/(NX-1); %[m]
DELTAY=(YMAX-YMIN)/(NY-1); %[m]


density = 2700.d0;
cp = 3000.d0;
cs = 1732.05d0;

Courant_number = cp * DELTAT * sqrt(1.d0/DELTAX^2.d0 + 1.d0/DELTAY^2.d0);
fprintf(' Courant number = %.4f\n',Courant_number); 
if Courant_number > 1.d0 
  disp('Error. Time step is too large, simulation will be unstable.');
  %break;
end

wavelength = cp/f0;
nodes_per_wavelength_x = wavelength/DELTAX;
nodes_per_wavelength_y = wavelength/DELTAY;
fprintf(' Shortest wavelength = %.2f m\n Nodes per wavelength:\n \t %.2f OX\n \t %.2f OY\n', wavelength, nodes_per_wavelength_x, nodes_per_wavelength_y);

if nodes_per_wavelength_x < 10 || nodes_per_wavelength_y < 10
    disp('Too few nodes per wavelength. Decrease f0 or increase NX and NY');
end

%--------------------------------------------------------------------------
%---------------------- FLAGS ---------------------------------------------

% display information on the screen from time to time
IT_DISPLAY = 100;

%Take instant snapshotfalse
SNAPSHOT=false;
snapshot_time=1200:100:1600; %on what steps

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

DISP_NORM=true;

DATA_TO_BINARY_FILE=false;
tag='mz_';

RED_BLUE=false;
COLORBAR_ON=true;
FE_BOUNDARY=false;
SHOW_REC_POSITIONS=true;

SAVE_SEISMOGRAMS=true;

seis_tag=['BB'];
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
eps=0.00000000001d0;

nx_vec=[0:NX]*DELTAX;	%[m]
ny_vec=[0:NY]*DELTAY;

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

xsource=round(NX/2)*DELTAX;
ysource=round(NY/2)*DELTAY;

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
% xsource=ISOURCE*DELTAX;
% ysource=JSOURCE*DELTAY;
xsource = 1000;
ysource = 1000;

ux=zeros(3,NX+1,NY+1);
uy=zeros(3,NX+1,NY+1);

 %----------------------------------------
 %--- program starts here ----------------
 %----------------------------------------

fprintf('2D elastic finite-difference code in displacement formulation with C-PML\n\n');
fprintf('NX = %d\n',NX);
fprintf('NY = %d\n\n',NY);
fprintf('size of the model along X = %.2f\n',NX*DELTAX);
fprintf('size of the model along Y = %.2f\n\n',NY*DELTAY);
fprintf('Total number of grid points = %.2f\n\n',NX * NY);

if SAVE_SEISMOGRAMS
    fprintf('ON. Save seismograms\n');
    NREC=2;
    fprintf('Set %d recievers:\n',NREC);

    xdeb=500;
    xfin=1000;
    ydeb=1000;
    yfin=1000;
    fprintf('  x0=%.2f  x1=%.2f\n  y0=%.2f  y1=%.2f\n', xdeb,xfin,ydeb,yfin);


    % for receivers
    ix_rec=zeros(NREC,1);
    iy_rec=zeros(NREC,1);
    xrec=zeros(NREC,1);
    yrec=zeros(NREC,1); 

    % for seismograms
    seisux=zeros(NSTEP,NREC);
    seisuy=zeros(NSTEP,NREC);

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

rho=zeros(NX+1,NY+1);
C=zeros(NX+1,NY+1,4);

% lambda = 20.d0;
% mu = 10.d0;


lambda =density*(cp*cp - 2.d0*cs*cs);
mu = density*cs*cs;

c11 = (lambda + 2.d0*mu);
c13 = lambda;
c33 = c11;
c44 = mu;    
 
fprintf('\nCreate C 6D %d elements\n',NX*NY*4);
for i = 1:NX
    for j = 1:NY
            C(i,j,:)=[c11 c13 c33 c44];
            rho(i,j) = density;
    end
end

% clearvars densitya cpa csa lambdaa mua densityb cpb csb lambdab mub;
% clearvars x_trial y_trial topo_szx tgrx;
% clearvars c11a c13a c33a c44a c11b c13b c33b c44b;
fprintf('C(i,j,4) of size: %s  ...OK\n',num2str(size(C)));

fprintf('\n');


% clearvars densitya cpa csa lambdaa mua densityb cpb csb lambdab mub;
% clearvars x_trial y_trial topo_szx tgrx;
% clearvars c11a c13a c33a c44a c11b c13b c33b c44b;

dx=DELTAX; 
dy=DELTAY;
dx2=DELTAX^2.d0;
dy2=DELTAY^2.d0;
ddx=2.d0*DELTAX;
ddy=2.d0*DELTAY;
dxdy4=4.d0*DELTAX*DELTAY;

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
            
            value_dux_dxx = (ux(2,i-1,j) - 2.d0*ux(2,i,j) + ux(2,i+1,j)) / dx2;
            value_duy_dxx = (uy(2,i-1,j) - 2.d0*uy(2,i,j) + uy(2,i+1,j)) / dx2;
            
            value_dux_dyy = (ux(2,i,j-1) - 2.d0*ux(2,i,j) + ux(2,i,j+1)) / dy2;
            value_duy_dyy = (uy(2,i,j-1) - 2.d0*uy(2,i,j) + uy(2,i,j+1)) / dy2;

            value_dux_dxy = (ux(2,i+1,j+1) - ux(2,i+1,j-1) - ux(2,i-1,j+1) + ux(2,i-1,j-1)) / dxdy4;
            value_duy_dxy = (uy(2,i+1,j+1) - uy(2,i+1,j-1) - uy(2,i-1,j+1) + uy(2,i-1,j-1)) / dxdy4;            

            value_dux_dyx = value_dux_dxy;
            value_duy_dyx = value_duy_dxy;
            
            c11v=C(i,j,1);
            c13v=C(i,j,2);
            c33v=C(i,j,3);
            c44v=C(i,j,4);

            dt2rho=(DELTAT^2.d0)/rhov;

            sigmas_ux = c11v * value_dux_dxx + c13v * value_duy_dyx + c44v * value_dux_dyy + c44v * value_duy_dxy;
            sigmas_uy = c44v * value_dux_dyx + c44v * value_duy_dxx + c13v * value_dux_dxy + c33v * value_duy_dyy;

            ux(3,i,j) = 2.d0 * ux(2,i,j) - ux(1,i,j) + sigmas_ux * dt2rho;
            uy(3,i,j) = 2.d0 * uy(2,i,j) - uy(1,i,j) + sigmas_uy * dt2rho;
   
        end
    end

    t = double(it-1)*DELTAT;
%     if t<=t0             
        % add the source (force vector located at a given grid point)
        a = pi*pi*f0*f0;
        % Gaussian
        %source_term = factor * exp(-a*(t-t0)^2);
        %source_term = factor * (t-t0);  
        % first derivative of a Gaussian
        %       source_term =  -factor*2.d0*a*(t-t0)*exp(-a*(t-t0)^2);
        % Ricker source time function (second derivative of a Gaussian)
        source_term = -factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
        %     source_term=factor * exp(-a*(t-t0)^2);
        force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
        force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
        %       force_x=factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
        %       force_y=factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
        % define location of the source
        i = ISOURCE;
        j = JSOURCE;

        rhov=rho(i,j);
%         ux(3,i,j) = ux(3,i,j) + (force_x * DELTAT)/rhov;
%         uy(3,i,j) = uy(3,i,j) + (force_y * DELTAT)/rhov;
        ux(3,i,j) = (force_x * DELTAT^2.0)/rhov;
        uy(3,i,j) = (force_y * DELTAT^2.0)/rhov;
%     fprintf('%e \n', force_y * DELTAT / rhov)
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
    
    if DISP_NORM
        unorm = sqrt(ux(3,:,:).^2+uy(3,:,:).^2);
    end
    
    % store seismograms
    if SAVE_SEISMOGRAMS
        for irec = 1:NREC
                seisux(it,irec) = ux(3,ix_rec(irec),iy_rec(irec));
                seisuy(it,irec) = uy(3,ix_rec(irec),iy_rec(irec));
%                 if DISP_NORM
%                     seisuy(it,irec) = unorm(1,ix_rec(irec),iy_rec(irec));
%                 end
        end   
    end
    
    %Set previous timesteps
    ux(1,:,:)=ux(2,:,:);
    ux(2,:,:)=ux(3,:,:);
    
    uy(1,:,:)=uy(2,:,:);
    uy(2,:,:)=uy(3,:,:);

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
                u=unorm;
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
            
            if SHOW_REC_POSITIONS && SAVE_SEISMOGRAMS
                for i=1:NREC
                    scatter(xrec(i), yrec(i),'filled','r');hold on;
                end
                drawnow; 
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
          
            if DATA_TO_BINARY_FILE
                filename=[tag 'u_' 'disp_t_' num2str(it) '.txt'];
                dlmwrite(filename, u);
                fprintf('Data file %s saved to %s\n',filename, pwd);
            end
        end
        fprintf('\n'); 
    end
    if PAUSE_ON
        pause(pause_time);
    end
end
  % end of time loop
  
  current_folder=pwd;	%current path

if SAVE_SEISMOGRAMS
      for i=1:NREC
%           filename=[seis_tag '_XX' '_' num2str(i) '.txt'];
          filename=[seis_tag '.S' sprintf('%04d',i) '.BXX.semd'];
          filename = ['../OUTPUT_FILES/' filename];
%           dlmwrite(filename, [time_vec, seisux(:,i)],'Delimiter','\t');
          dlmwrite(filename, [time_vec, seisux(:,i)],'Delimiter','\t', 'precision', '%E');
          fprintf('ux. Seismogram for rec at %.2f %.2f saved as %s to %s\n', (ix_rec(i)-1)*dx, (iy_rec(i)-1)*dy, filename, pwd);
%           filename=[seis_tag '_XZ' '_' num2str(i) '.txt'];
          filename = [seis_tag '.S' sprintf('%04d',i) '.BXZ.semd'];
          filename = ['../OUTPUT_FILES/' filename];
          dlmwrite(filename, [time_vec, seisuy(:,i)],'Delimiter','\t', 'precision', '%E');
          fprintf('uy. Seismogram for rec at %.2f %.2f saved as %s to %s\n', (ix_rec(i)-1)*dx, (iy_rec(i)-1)*dy, filename, pwd);
      end
 end
 
  disp('End');
 