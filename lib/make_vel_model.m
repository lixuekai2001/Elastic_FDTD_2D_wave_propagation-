%% Create arbitrary velocity model of cp, cs, rho, with arbitrary number of interfaces of orbitrary length

function [model_cp, model_cs, model_rho, interface_list] = make_vel_model(NX, NZ, XMAX, XMIN, ZMAX, ZMIN)
    
    n_interfaces = 5;
    
    DELTAX=(XMAX-XMIN)/(NX-1);  %[m]
    DELTAZ=(ZMAX-ZMIN)/(NZ-1);  %[m]
    
%     x_vec=[0:NX-1]*DELTAX;	%[m]
%     z_vec=[0:NZ-1]*DELTAZ;    %[m

    x_vec = linspace(XMIN,XMAX,40*NX);
    
    x_vec1 = 40*x_vec;
    x_vec2 = 40*x_vec;
    x_vec3 = 40*x_vec;
    x_vec4 = 40*x_vec;
    x_vec5 = x_vec(round(length(x_vec)/2)+5:end-5);
    
    int_1 = DELTAZ*(3*NZ/4+1) + 0*x_vec1;
    int_2 = DELTAZ*(3*NZ/5+1) + 0.03*ZMAX*sin(2*pi*x_vec2/XMAX);
    int_3 = DELTAZ*(2*NZ/4+1) + 0.05*ZMAX*sin(2*pi*x_vec3/XMAX);
    int_4 = DELTAZ*(1*NZ/4+1) + 0.08*ZMAX*sin(2*pi*x_vec4/XMAX + pi/2);
    int_5 = 0.2*ZMAX*(1 - sin(2*2*pi*x_vec5/XMAX+pi/2));


    model_rho = 1000*ones(NZ,NX);
    model_cp = 1500.d0*ones(NZ,NX);
    model_cs = zeros(NZ,NX);

    % Interface 1
    [model_cp, model_cs, model_rho, x_vec1, int_1] = make_model_under_curve(model_cp, model_cs, model_rho, 2500.d0, 2500.d0/1.732d0, ... 
                                                            2200.d0, x_vec1, int_1, DELTAX, DELTAZ, NX, NZ);
    % Interface 2    
    [model_cp, model_cs, model_rho, x_vec2, int_2] = make_model_under_curve(model_cp, model_cs, model_rho, 2000.d0, 2000.d0/1.732d0, ... 
                                                            2000.d0, x_vec2, int_2, DELTAX, DELTAZ, NX, NZ);
    % Interface 3
    [model_cp, model_cs, model_rho, x_vec3, int_3] = make_model_under_curve(model_cp, model_cs, model_rho, 3500.d0, 3500.d0/1.732d0, ... 
                                                            2600.d0, x_vec3, int_3, DELTAX, DELTAZ, NX, NZ);
    % Interface 4
    [model_cp, model_cs, model_rho, x_vec4, int_4] = make_model_under_curve(model_cp, model_cs, model_rho, 3000.d0, 3000.d0/1.732d0, ... 
                                                            2300.d0, x_vec4, int_4, DELTAX, DELTAZ, NX, NZ);
    % Interface 5    
    [model_cp, model_cs, model_rho, x_vec5, int_5] = make_model_under_curve(model_cp, model_cs, model_rho, 4500.d0, 4500.d0/1.732d0, ... 
                                                        2800.d0, x_vec5, int_5, DELTAX, DELTAZ, NX, NZ);

    % Fill up interfaces list
    interface_list = cell(n_interfaces,2);

    interface_list{1,1} = x_vec1;
    interface_list{1,2} = int_1;
    
    interface_list{2,1} = x_vec2;
    interface_list{2,2} = int_2;
    
    interface_list{3,1} = x_vec3;
    interface_list{3,2} = int_3;
    
    interface_list{4,1} = x_vec4;
    interface_list{4,2} = int_4;
    
    interface_list{5,1} = x_vec5;
    interface_list{5,2} = int_5;

    
%      imagesc(x_vec, z_vec,model_cp'); hold on;
%      set(gca,'YDir','normal');
%      xlabel('m', 'fontSize',14, 'fontweight','b');
%      ylabel('m', 'fontSize',14, 'fontweight','b');
%      title('Model. Cp + Interfaces', 'fontSize',18);
%      colorbar();
end