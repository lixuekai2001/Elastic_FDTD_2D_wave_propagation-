%% Elastic VTI solver with matrix 4x4 as input

function [ux, uy] = solver_mx_VTI_elastic(ux, uy, DELTAT, coeffux, coeffuy,rho,C)

    ux(3,:,:)=0.d0;
    uy(3,:,:)=0.d0;
    
    [~, NX, NY] = size(ux);
    
    for i = 2:NX-1
        for j = 2:NY-1
            
            rhov=rho(i,j);
            
            A_ux=coeffux{i,j};
            value_dux_dxx=A_ux(1,1:3)*[ux(2,i-1,j); ux(2,i,j); ux(2,i+1,j)];
            value_dux_dyy=A_ux(2,1:3)*[ux(2,i,j-1); ux(2,i,j); ux(2,i,j+1)];
            value_dux_dxy=A_ux(3,:)*[ux(2,i+1,j+1); ux(2,i+1,j-1); ux(2,i-1,j+1); ux(2,i-1,j-1)];
            value_dux_dyx=A_ux(4,:)*[ux(2,i+1,j+1); ux(2,i+1,j-1); ux(2,i-1,j+1); ux(2,i-1,j-1)];
            
            A_uy=coeffuy{i,j};
            value_duy_dxx=A_uy(1,1:3)*[uy(2,i-1,j); uy(2,i,j); uy(2,i+1,j)];
            value_duy_dyy=A_uy(2,1:3)*[uy(2,i,j-1); uy(2,i,j); uy(2,i,j+1)];
            value_duy_dxy=A_uy(3,:)*[uy(2,i+1,j+1); uy(2,i+1,j-1); uy(2,i-1,j+1); uy(2,i-1,j-1)];
            value_duy_dyx=A_uy(4,:)*[uy(2,i+1,j+1); uy(2,i+1,j-1); uy(2,i-1,j+1); uy(2,i-1,j-1)];

            %--------------------------------------------------------------------------------------------------------------------
            
            dt2rho=(DELTAT^2.d0)/rhov;
% 
%           sigmas_ux= c11v * value_dux_dxx + c13v * value_duy_dyx + c44v * value_dux_dyy + c44v * value_duy_dxy;
%           sigmas_uy= c44v * value_dux_dyx + c44v * value_duy_dxx + c13v * value_dux_dxy + c33v * value_duy_dyy;

            sigmas_ux= value_dux_dxx + value_duy_dyx + value_dux_dyy + value_duy_dxy;
            sigmas_uy= value_dux_dyx + value_duy_dxx + value_dux_dxy + value_duy_dyy;

            ux(3,i,j) = 2.d0 * ux(2,i,j) - ux(1,i,j) + sigmas_ux * dt2rho;
            uy(3,i,j) = 2.d0 * uy(2,i,j) - uy(1,i,j) + sigmas_uy * dt2rho;

        end
    end
end