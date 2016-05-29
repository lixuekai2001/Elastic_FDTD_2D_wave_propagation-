%% Fill data matrix under given curve

function [model_cp, model_cs, model_rho, xdscr, zdscr] = make_model_under_curve(model_cp, model_cs, model_rho, val_cp, val_cs, val_rho, x_vec, z_vec, DELTAX, DELTAZ, NX, NZ)
    
    x_range = floor(x_vec/DELTAX+1);    % find nearest x nodes to each curve point
    x_range(x_range == 0) = 1;          % if there are 0s(matlab counts from 1) change them to 1s
    x_range(x_range > NX) = NX;         % if > NX change them to NX
    x_range_min = min(min(x_range));    % find the left border
    x_range_max = max(max(x_range));    % find the right border
   
    xdscr = [];
    zdscr = [];
    
    one_vec = ones(NZ,1);
    z_trial = [0:NZ-1]*DELTAZ;
    
    for i = x_range_min:x_range_max
        x_trial = (i-1)*DELTAX*one_vec;
        [xi,zi]=curveintersect(x_trial,z_trial,x_vec, z_vec);
        
        if ~isempty(zi)
            j = floor(zi/DELTAZ)+1;
            model_cp(i,1:j) = val_cp;
            model_cs(i,1:j) = val_cs;
            model_rho(i,1:j) = val_rho;
            
            xdscr = [xdscr xi];
            zdscr = [zdscr zi];
        end
    end
end
