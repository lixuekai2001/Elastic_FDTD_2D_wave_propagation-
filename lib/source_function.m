%% 2D source function

function [force_x, force_y] = source_function(f0, t0, factor, ANGLE_FORCE, t)

    PI = 3.141592653589793238462643d0;     % value of PI
    DEGREES_TO_RADIANS = PI / 180.d0;       % conversion from degrees to radians

    a = pi*pi*f0*f0;

    % Gaussian:
    %       source_term = factor * exp(-a*(t-t0)^2);
      
    % First derivative of a Gaussian:
    %       source_term =  -factor*2.d0*a*(t-t0)*exp(-a*(t-t0)^2);
    
    % Ricker source time function (second derivative of a Gaussian):
    %       source_term = -factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
    
    source_term =  -factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
    force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
    force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
    
end