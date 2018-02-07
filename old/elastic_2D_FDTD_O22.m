% Simple FDTD seismic wave propagation in 2D elastic isotropic medium 
% surrounded by simple sponge boundaries with exponential decay 
% (Cerjan, 1985).
%
% We solve wave equation in time domain and displacement formulation
% getting wavefield in terms of displacement vector [ux, uz].
%
% Elastic medium is parametrized by Lame parameters and density, we show
% Courant condition and number of points per wavelength prior running %
% loop over time steps.
%
% Conventional FD star-stencils deliver accuracy O(2,2)
% [1 -2 1]/dx2 and [1 -1 -1 1]/4dxdz
%
% --------------------------------------------------------------
% The code is intentionally writen in a single file
% to simplify start up. MATLAB is designed to operate matrices
% but it is naturally easier to understand loops.
%
% The program does not save any files, add such option manually.
%
% The goal is to provide a simple example code of wave propagation
% in elastic medium.
% --------------------------------------------------------------
% Oleg Ovcharenko, 2017
% oleg.ovcharenko@kaust.edu.sa
%
% King Abdullah University of Science and Technology
% Thuwal, Saudi Arabia
% --------------------------------------------------------------

close all;
clear all;

% Output periodicity in time steps
IT_DISPLAY = 20;

%% MODEL

% Model dimensions, [m]
zmax = 2000.0;
xmax = 2000.0;
xmin = 0.0;
zmin = 0.0;

% Grid dimensions
nx = 201;
nz = 201;

% Elastic parameters
vp = 3300.0 * ones(nz, nx);         % velocity of compressional waves, [m/s]
vs = vp / 1.732;                    % velocity of shear waves, [m/s]
rho = 2800.0 * ones(size(vp));      % density, [kg/m3]

% Uncomment below to load Marmousi II benchmark model
% load('./data/marmousi_ii_all_1_25m.mat','vp','vs','rho', 'xmax', 'zmax')
% vp = imresize(vp, [nz nx]);
% vs = imresize(vs, [nz nx]);
% rho = imresize(rho, [nz nx]) * 1000;
% xmax = xmax * 1000;
% zmax = zmax * 1000;

% Grid step, [m]
dx = (xmax-xmin)/(nx-1);
dz = (zmax-zmin)/(nz-1);

% Lame parameters
lambda = rho.*(vp.^2 - 2*vs.^2);    % first Lame parameter
muu = rho.*vs.^2;                   % shear modulus, [N/m2]

%% TIME STEPPING

t_total = 0.6;                      % [sec] recording duration
% dt = 2e-3;                          % [sec] time step
dt = 1/(max(vp(:)) * sqrt(1.0/dx^2 + 1.0/dz^2));
nt = round(t_total/dt);             % number of time steps

%% SOURCE

f0 = 10.0;                          % dominant frequency of the wavelet
t0 = 1.20 / f0;                     % excitation time
factor = 1e10;                      % amplitude coefficient
angle_force = 90.0;                 % spatial orientation

isrc = round(nx/2);                 % source location along OX
jsrc = round(nz/2);                 % source location along OZ

%% ABSORBING BOUNDARY (ABS)

abs_thick = floor(0.15*nx);         % thicknes of the layer
abs_rate = 0.015*15/abs_thick;      % decay rate

%% SUMMARY

fprintf('2D elastic finite-difference code in displacement formulation with Cerjan boundary conditions\n\n');
fprintf('Model %dx%d, dx = %.2f  dz=%.2f  dt=%.2e \n',nx, nz, dx,dz,dt);
fprintf('%.2f x %.2f m\n',nx*dx, nz*dz);
fprintf('Allocated memory: %.2f mb\n\n', monitor_memory_whos);
% Check Courant stability condition
fprintf('CFL number = %.2f\n', max(vp(:))*dt * sqrt(1.0/dx^2 + 1.0/dz^2));
% Check number of nodes per wavelength, ~ 10 is recommended
fprintf('Shortest wavelength = %.2f m\nNodes per wavelength:\n \t %d OX\n \t %d OY\n', min(vp(:))/f0, floor(min(vp(:))/f0/dx), floor(min(vp(:))/f0/dz));

%% MODELLING LOOP

ux = zeros(3,nx+2,nz+2);            % Wavefields for time steps
uz = zeros(3,nx+2,nz+2);            % t-2 (1), t-1 (2) and t (3)
nx_vec = [0:nx]*dx;                 % vectors with axis for plottig
nz_vec = [0:nz]*dz;

dxdy4 = 4.0 * dx * dz;              % denominators for stencils
dx2 = dx^2.0;                       % assigned here to slightly reduce
dz2 = dz^2.0;                       % computations

for it = 1:nt                       % Loop over TIME
    ux(3,:,:)=0.0;                  % init wavefield at t
    uz(3,:,:)=0.0;
    for i = 2:nx                    % Loop over OX
        for j = 2:nz                % Loop over OZ
            lam = lambda(i,j);
            mu = muu(i,j);
            lam_2mu = lam + 2 * mu;
            
            dux_dxx = [1 -2 1]/dx2 * [ux(2,i-1,j); ux(2,i,j); ux(2,i+1,j)];
            dux_dzz = [1 -2 1]/dz2 * [ux(2,i,j-1); ux(2,i,j); ux(2,i,j+1)];
            dux_dxz = [1 -1 -1 1]/dxdy4 * [ux(2,i+1,j+1); ux(2,i+1,j-1); ux(2,i-1,j+1); ux(2,i-1,j-1)];
            dux_dzx = [1 -1 -1 1]/dxdy4 * [ux(2,i+1,j+1); ux(2,i+1,j-1); ux(2,i-1,j+1); ux(2,i-1,j-1)];
            
            duz_dxx = [1 -2 1]/dx2 * [uz(2,i-1,j); uz(2,i,j); uz(2,i+1,j)];
            duz_dzz = [1 -2 1]/dz2 * [uz(2,i,j-1); uz(2,i,j); uz(2,i,j+1)];
            duz_dxz = [1 -1 -1 1]/dxdy4 * [uz(2,i+1,j+1); uz(2,i+1,j-1); uz(2,i-1,j+1); uz(2,i-1,j-1)];
            duz_dzx = [1 -1 -1 1]/dxdy4 * [uz(2,i+1,j+1); uz(2,i+1,j-1); uz(2,i-1,j+1); uz(2,i-1,j-1)];
            
            sigmas_ux = lam_2mu * dux_dxx + lam * duz_dzx + mu * (dux_dzz + duz_dxz);
            sigmas_uz = mu * (dux_dzx + duz_dxx) + lam * dux_dxz + lam_2mu * duz_dzz;
            
            dt2rho=(dt^2.d0)/rho(i,j);
            ux(3,i,j) = 2.d0 * ux(2,i,j) - ux(1,i,j) + sigmas_ux * dt2rho;
            uz(3,i,j) = 2.d0 * uz(2,i,j) - uz(1,i,j) + sigmas_uz * dt2rho;
        end
    end
    
    % Add volumetric force source term
    t = double(it-1)*dt;
    deg2rad = pi / 180.d0;
    a = pi*pi*f0*f0;
    %     source_term = factor * exp(-a*(t-t0)^2);                              % Gaussian
    %     source_term =  -factor*2.0*a*(t-t0)*exp(-a*(t-t0)^2);                 % First derivative of a Gaussian:
    source_term = -factor * (1.0 - 2.0*a*(t-t0)^2)*exp(-a*(t-t0)^2);        % Ricker source time function (second derivative of a Gaussian):
    force_x = sin(angle_force * deg2rad) * source_term;
    force_y = cos(angle_force * deg2rad) * source_term;
    
    ux(3, isrc, jsrc) = ux(3, isrc, jsrc) + force_x * dt^2.d0 / rho(isrc, jsrc);
    uz(3, isrc, jsrc) = uz(3, isrc, jsrc) + force_y * dt^2.d0 / rho(isrc, jsrc);
   
    % Start ABS
    lmargin = [abs_thick abs_thick];
    rmargin = lmargin;
    weights = ones(nx,nz);
    for iz = 2:nz
        for ix = 2:nx
            i = 0;
            j = 0;
            k = 0;
            if (ix < lmargin(1) + 1)
                i = lmargin(1) + 1 - ix;
            end
            if (iz < lmargin(2) + 1)
                k = lmargin(2) + 1 - iz;
            end
            if (nx - rmargin(1) < ix)
                i = ix - nx + rmargin(1);
            end
            if (nz - rmargin(2) < iz)
                k = iz - nz + rmargin(2);
            end
            if (i == 0 && j == 0 && k == 0)
                continue
            end
            rr = abs_rate * abs_rate * double(i*i + j*j + k*k );
            weights(ix,iz) = exp(-rr);
        end
    end
    ux(3,2:end-1,2:end-1) = squeeze(ux(3,2:end-1,2:end-1)) .* weights;
    ux(2,2:end-1,2:end-1) = squeeze(ux(2,2:end-1,2:end-1)) .* weights;
    ux(1,2:end-1,2:end-1) = squeeze(ux(1,2:end-1,2:end-1)) .* weights;
    uz(3,2:end-1,2:end-1) = squeeze(uz(3,2:end-1,2:end-1)) .* weights;
    uz(2,2:end-1,2:end-1) = squeeze(uz(2,2:end-1,2:end-1)) .* weights;
    uz(1,2:end-1,2:end-1) = squeeze(uz(1,2:end-1,2:end-1)) .* weights;
    % End of ABS
    
    % Exchange data between t-2 (1), t-1 (2) and t (3)
    ux(1,:,:) = ux(2,:,:);
    ux(2,:,:) = ux(3,:,:);
    uz(1,:,:) = uz(2,:,:);
    uz(2,:,:) = uz(3,:,:);
    
    % output information
    if mod(it,IT_DISPLAY) == 0
        tc = single((it-1)*dt);
        fprintf('Time step: %d \t %.4f s\n',it, single(tc));
        u=sqrt(ux(3,:,:).^2 + uz(3,:,:).^2); u = squeeze(u(1,:,:))';
        imagesc(nx_vec, nz_vec, u);  colorbar; colormap jet;
        title(['Step = ',num2str(it),' Time: ',sprintf('%.4f',tc),' sec']);
        xlabel('m'); ylabel('m');
        axis equal tight; 
        drawnow;
        
        
%     surf(u); shading interp; lighting phong; colormap hot; axis off; zlim([amp_src_min amp_src_max]);
%     set(gcf,'Color', [0 0 0], 'Name', sprintf('Tiny FDTD, step = %i', it));
%     title(sprintf('FDTD Time = %.2f msec',tc),'Color',[1 0 0],'FontSize', 22); 
%     view(-12, 65); colorbar; drawnow;

    end
end
disp('End');
