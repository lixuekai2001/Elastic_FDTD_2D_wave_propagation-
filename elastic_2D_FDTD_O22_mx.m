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
% --------------------------------------------------------------
% The code is intentionally writen in a single file
% to simplify start up.
%
% The program does not save any files, add such option manually if needed.
%
% The goal is to provide a simple example of wave propagation
% in elastic medium.
%
% --------------------------------------------------------------
% Oleg Ovcharenko and Vladimir Kazei, 2018
% oleg.ovcharenko@kaust.edu.sa
% vladimir.kazei@kaust.edu.sa
%
% King Abdullah University of Science and Technology
% Thuwal, Saudi Arabia
% --------------------------------------------------------------

close all;
clear all;

% Output periodicity in time steps
IT_DISPLAY = 100;

%% MODEL
% Model dimensions, [m]
% zmax = 2000.0;
% xmax = 2000.0;
% xmin = 0.0;
% zmin = 0.0;

% Grid dimensions
nx = 1001;
nz = 201;

dx = 10;
dz = dx;

% Elastic parameters
vp = 3300.0 * ones(nz, nx);         % velocity of compressional waves, [m/s]
vs = vp / 1.732;                    % velocity of shear waves, [m/s]
rho = 2800.0 * ones(size(vp));      % density, [kg/m3]

% % Uncomment below to load Marmousi II benchmark model
% load('./data/marmousi_ii_all_1_25m.mat')
% krefine = 0.1;
% vp = imresize(vp, krefine);
% vs = imresize(vs, krefine);
% rho = imresize(rho, krefine) * 1000;
% xmax = 1000 * xmax * krefine;
% zmax = 1000 * zmax * krefine;
% dx = 2 * 1000 * dx * krefine;
% dz = 2 * 1000 * dz * krefine;
% [nz, nx] = size(vp);

% Grid step, [m]
% dx = (xmax-xmin)/(nx-1);
% dz = (zmax-zmin)/(nz-1);

% Lame parameters
lam = rho.*(vp.^2 - 2*vs.^2);    % first Lame parameter
mu = rho.*vs.^2;                   % shear modulus, [N/m2]

%% TIME STEPPING
t_total = 1.;                      % [sec] recording duration
dt = 1/(max(vp(:)) * sqrt(1.0/dx^2 + 1.0/dz^2));
nt = round(t_total/dt);             % number of time steps
t = [0:nt]*dt;

%% SOURCE
f0 = 10.0;                          % dominant frequency of the wavelet
t0 = 1.20 / f0;                     % excitation time
factor = 1e10;                      % amplitude coefficient
angle_force = 90.0;                 % spatial orientation

jsrc = round(nz/2);                 % source location along OZ
isrc = round(nx/2);                 % source location along OX

deg2rad = pi / 180.d0;              % convet degrees to radians
a = pi*pi*f0*f0;
dt2rho_src = dt^2/rho(jsrc, isrc);
%     source_term = factor * exp(-a*(t-t0)^2);                             % Gaussian
%     source_term =  -factor*2.0*a*(t-t0)*exp(-a*(t-t0)^2);                % First derivative of a Gaussian:
source_term = -factor * (1.0 - 2.0*a*(t-t0).^2).*exp(-a*(t-t0).^2);        % Ricker source time function (second derivative of a Gaussian):

force_x = sin(angle_force * deg2rad) * source_term * dt2rho_src;
force_z = cos(angle_force * deg2rad) * source_term * dt2rho_src;

%% ABSORBING BOUNDARY (ABS)
abs_thick = min(floor(0.20*nx), floor(0.20*nz));         % thicknes of the layer
abs_rate = 0.3/abs_thick;      % decay rate

lmargin = [abs_thick abs_thick];
rmargin = lmargin;
weights = ones(nz+2,nx+2);
for iz = 1:nz+2
    for ix = 1:nx+2
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
        weights(iz,ix) = exp(-rr);
    end
end

%% SUMMARY
fprintf('#################################################\n');
fprintf('2D elastic FDTD seismic wave propagation in \ndisplacement formulation with Cerjan(1985) \nboundary conditions\n');
fprintf('#################################################\n');
fprintf('Model:\n\t%dx%d\tgrid\n\t%.1fx%.1f\tdx[m]\n',nx, nz, dx,dz);
fprintf('\t%.2fx%.2f\t[km]\n',nx*dx/1000, nz*dz/1000);
fprintf('\t%.2f...%.2f\t[km/s] vp\n', min(vp(:)), max(vp(:)));
fprintf('\t%.2f...%.2f\t[km/s] vs\n', min(vs(:)), max(vs(:)));
fprintf('\t%.2f...%.2f\t[kg/m3] rho\n', min(rho(:)), max(rho(:)));
fprintf('Other:\n\t%.2f\tCFL number\n', max(vp(:))*dt * sqrt(1.0/dx^2 + 1.0/dz^2));
fprintf('\t%.2f\t[m] shortest wavelength\n\t%d, %d\t PPWavelength OX, OZ\n', min(vp(:))/f0, floor(min(vp(:))/f0/dx), floor(min(vp(:))/f0/dz));

%% ALLOCATE MEMORY FOR WAVEFIELD
ux3 = zeros(nz+2,nx+2);            % Wavefields at t
uz3 = zeros(nz+2,nx+2);            
ux2 = zeros(nz+2,nx+2);            % Wavefields at t-1
uz2 = zeros(nz+2,nx+2);            
ux1 = zeros(nz+2,nx+2);            % Wavefields at t-2
uz1 = zeros(nz+2,nx+2);            

% Coefficients for derivatives
co_dxx = 1/dx^2;
co_dzz = 1/dz^2; 
co_dxz = 1/(4.0 * dx * dz);
co_dzx = 1/(4.0 * dx * dz);

dt2rho=(dt^2)./rho;
lam_2mu = lam + 2 * mu;

%% Loop over TIME
tic;
for it = 1:nt                       
    ux3 = zeros(size(ux2));              
    uz3 = zeros(size(uz2));
    % Second-order derivatives
    % Ux
    dux_dxx = co_dxx * (ux2(1:end-2,2:end-1) - 2*ux2(2:end-1,2:end-1) + ux2(3:end,2:end-1));
    dux_dzz = co_dzz * (ux2(2:end-1,1:end-2) - 2*ux2(2:end-1,2:end-1) + ux2(2:end-1,3:end));
    dux_dxz = co_dxz * (ux2(3:end,1:end-2) - ux2(3:end,3:end) ...
                      - ux2(1:end-2,1:end-2) + ux2(1:end-2,3:end));
    dux_dzx = co_dzx * (ux2(3:end,1:end-2) - ux2(3:end,3:end) ...
                      - ux2(1:end-2,1:end-2) + ux2(1:end-2,3:end));
    % Uz             
    duz_dxx = co_dxx * (uz2(1:end-2,2:end-1) - 2*uz2(2:end-1,2:end-1) + uz2(3:end,2:end-1));
    duz_dzz = co_dzz * (uz2(2:end-1,1:end-2) - 2*uz2(2:end-1,2:end-1) + uz2(2:end-1,3:end));
    duz_dxz = co_dxz * (uz2(3:end,1:end-2) - uz2(3:end,3:end) ...
                      - uz2(1:end-2,1:end-2) + uz2(1:end-2,3:end));
    duz_dzx = co_dzx * (uz2(3:end,1:end-2) - uz2(3:end,3:end) ...
                      - uz2(1:end-2,1:end-2) + uz2(1:end-2,3:end));
    % Stress G             
    sigmas_ux = lam_2mu .* dux_dxx + lam .* duz_dzx + mu .* (dux_dzz + duz_dxz);
    sigmas_uz = mu .* (dux_dzx + duz_dxx) + lam .* dux_dxz + lam_2mu .* duz_dzz;
    % U(t) = 2*U(t-1) - U(t-2) + G dt2/rho;
    ux3(2:end-1,2:end-1) = 2.0*ux2(2:end-1,2:end-1) - ux1(2:end-1,2:end-1) + sigmas_ux.*dt2rho;
    uz3(2:end-1,2:end-1) = 2.0*uz2(2:end-1,2:end-1) - uz1(2:end-1,2:end-1) + sigmas_uz.*dt2rho;
    % Add source term
    ux3(jsrc, isrc) = ux3(jsrc, isrc) + force_x(it);
    uz3(jsrc, isrc) = uz3(jsrc, isrc) + force_z(it);
    % Absorbing boundaries
    ux3 = ux3 .* weights;
    ux2 = ux2 .* weights;
    uz3 = uz3 .* weights;
    uz2 = uz2 .* weights;
    % Exchange data between t-2 (1), t-1 (2) and t (3)
    ux1 = ux2;
    ux2 = ux3;
    uz1 = uz2;
    uz2 = uz3;
    % Output
    if mod(it,IT_DISPLAY) == 0
        tc = single((it-1)*dt);
        fprintf('Time step: %d \t %.4f s\n',it, single(tc));
        u=sqrt(ux3.^2 + uz3.^2); imagesc(u);  
        %colorbar; colormap jet; 
%         title(['Step = ',num2str(it),' Time: ',sprintf('%.4f',tc),' sec']);
%         xlabel('m'); ylabel('m');
        axis equal tight; 
        drawnow;
    end
end
toc;
disp('End');
