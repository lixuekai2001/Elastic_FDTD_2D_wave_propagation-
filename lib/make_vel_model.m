function [model_cp, model_cs, model_rho] = make_vel_model(NX, NZ, x_vec, z_vec)
    zmax = max(z_vec);
    xmax = max(x_vec);
    z2_top =0*x_vec;
    % z2_top = 0.005*zmax*sin(2*pi*x_vec/xmax);
    z2_3_top = 0.003*zmax*sin(2*pi*x_vec/xmax);
    z3_top = 0.005*zmax*sin(2*pi*x_vec/xmax);
    z4_top = 0.008*zmax*sin(2*pi*x_vec/xmax + pi/2);
    z4_5_top = 0.04*zmax*sin(2*2*pi*x_vec/xmax+pi/2);
    z4_5_top(1:length(x_vec)/2) = 0;

    model_rho = 2200*ones(NZ,NX);

    model_cp = 1.5*ones(NZ,NX);
    for i=1:NX
        for j = 1:NZ
            if j>=NZ/4+1+z2_top(i)
                model_cp(i,j) = 2.5;
            end
        end
    end

    for i=1:NX
        for j = 1:NZ
            if j>=2*NZ/5+1+z2_3_top(i)
                model_cp(i,j) = 2.0;
            end
        end
    end

    for i=1:NX
        for j = 1:NZ
            if j>=2*NZ/4+1+z3_top(i)
                model_cp(i,j) = 3.5;
            end
        end
    end

    for i=1:NX
        for j = 1:NZ
            if j>=3*NZ/4+1+z4_top(i)
                model_cp(i,j) = 3.0;
            end
        end
    end

    for i=1:NX
        for j = 1:NZ
            if j>=NZ+1+z4_5_top(i)
                model_cp(i,j) = 4.5;
            end
        end
    end
    model_cs = model_cp/1.7320d0;

    
    imagesc(model_cs);
end