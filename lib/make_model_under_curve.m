%% Fill data matrix under given curve

function [model_cp, model_cs, model_rho] = make_model_under_curve(model_cp, model_cs, model_rho, val_cp, val_cs, val_rho, x_vec, interface, DELTAX, DELTAZ, NZ)

    ii = 0;
    for i=x_vec/DELTAX+1
        ii = ii + 1;
        for j = 1:NZ
            y_test = DELTAZ * j;
            if y_test <= interface(ii)
               model_cp(i,j) = val_cp;
               model_cs(i,j) = val_cs;
               model_rho(i,j) = val_rho;
            end
        end
    end
end