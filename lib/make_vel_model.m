function [model_cp, model_cs, model_rho, interface_list] = make_vel_model(NX, NZ, XMAX, XMIN, ZMAX, ZMIN)
    n_interfaces = 4;
    
    DELTAX=(XMAX-XMIN)/NX; %[m]
    DELTAZ=(ZMAX-ZMIN)/NZ; %[m]
    
    x_vec=[0:NX]*DELTAX;	%[m]
    z_vec=[0:NZ]*DELTAZ;
    
    zmax = max(z_vec);
    xmax = max(x_vec);
    
    int_1 = DELTAZ*(3*NZ/4+1) + 0*x_vec;
    int_2 = DELTAZ*(3*NZ/5+1) + 0.03*zmax*sin(2*pi*x_vec/xmax);
    int_3 = DELTAZ*(2*NZ/4+1) + 0.05*zmax*sin(2*pi*x_vec/xmax);
    int_4 = DELTAZ*(1*NZ/4+1) + 0.08*zmax*sin(2*pi*x_vec/xmax + pi/2);
    int_5 = 0.2*zmax*(1 - sin(2*2*pi*x_vec/xmax+pi/2));
    int_5(1:round(length(x_vec)/2)) = 0;

    model_rho = 2200*ones(NZ,NX);

    model_cp = 1.5*ones(NZ,NX);
    model_cs = zeros(NZ,NX);
    for i=1:NX
        for j = 1:NZ
            y_test = DELTAZ * j;
%             if j>=NZ/4+1+int_1(i)
            if y_test <= int_1(i)
                model_cp(i,j) = 2.5;
                model_cs(i,j) = model_cp(i,j)/1.7320d0;
            end
        end
    end
%     imagesc(model_cp); drawnow; input('?');
    
    for i=1:NX
        for j = 1:NZ
            y_test = DELTAZ * j;
%             if j>=2*NZ/5+1+int_2(i)
            if y_test <= int_2(i)
                model_cp(i,j) = 2.0;
                model_cs(i,j) = model_cp(i,j)/1.7320d0;
            end
        end
    end
%     imagesc(model_cp); drawnow; input('?');
    for i=1:NX
        for j = 1:NZ
            y_test = DELTAZ * j;
%             if j>=2*NZ/4+1+int_3(i)
            if y_test <= int_3(i)
                model_cp(i,j) = 3.5;
                model_cs(i,j) = model_cp(i,j)/1.7320d0;
            end
        end
    end
%     imagesc(model_cp); drawnow; input('?');
    for i=1:NX
        for j = 1:NZ
            y_test = DELTAZ * j;
%             if j>=3*NZ/4+1+int_4(i)
            if y_test <= int_4(i)
                model_cp(i,j) = 3.0;
                model_cs(i,j) = model_cp(i,j)/1.7320d0;
            end
        end
    end
%     imagesc(model_cp); drawnow; input('?');
    for i=1:NX
        for j = 1:NZ
            y_test = DELTAZ * j;
%             if j>=NZ+1+int_5(i)
            if y_test <= int_5(i)
                model_cp(i,j) = 4.5;
                model_cs(i,j) = model_cp(i,j)/1.7320d0;
            end
        end
    end

%         imagesc(model_cp); drawnow; input('?');
    % Fill up interfaces list
    interface_list = cell(n_interfaces,2);
    interface_list{1,1} = x_vec;
    interface_list{1,2} = int_1;
    
    interface_list{2,1} = x_vec;
    interface_list{2,2} = int_2;
    
    interface_list{3,1} = x_vec;
    interface_list{3,2} = int_3;
    
    interface_list{4,1} = x_vec;
    interface_list{4,2} = int_4;
    
    interface_list{5,1} = x_vec;
    interface_list{5,2} = int_5;

    
     imagesc(x_vec, z_vec,model_cs'); hold on;
     set(gca,'YDir','normal');
end