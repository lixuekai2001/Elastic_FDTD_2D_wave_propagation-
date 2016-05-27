%% Calculate maximum Ricker amplitude. [maxbar, minbar] = Ricker_amplitude(f0, t0, DELTAT, time_vec)

function [maxbar, minbar] = Ricker_amplitude(f0, t0, factor, ANGLE_FORCE, DELTAT, time_vec, rhov)
    % value of PI
    PI = 3.141592653589793238462643d0;

    % conversion from degrees to radians
    DEGREES_TO_RADIANS = PI / 180.d0;

    a = pi*pi*f0*f0;
    t = time_vec;
%     tbar = -factor*a*(time_vec-t0).*exp(-a*(time_vec-t0).^2)*DELTAT^2/rho;
    source_term =  -factor * (1.d0 - 2.d0*a*(t-t0).^2).*exp(-a*(t-t0).^2);
    force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
    force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
    
    disp_x  = (force_x * DELTAT^2.0)/rhov;
    disp_y  = (force_y * DELTAT^2.0)/rhov;
    
    maxbar = max(disp_y);
    minbar = min(disp_y);
    
%     maxbar = abs(minbar)+abs(maxbar);
%     minbar = 0;
    
clearvars tbar;
    fprintf('\n');
    fprintf('Colorbar:\n max: %f\n min: %f\n', maxbar, minbar);
end