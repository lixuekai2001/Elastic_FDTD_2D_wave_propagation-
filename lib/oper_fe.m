% %FE_BOUNDARY Solver
% UI - 12x1 - column of u and it's derivatives
% nvx and nvy - components of normal vector at chosen region
% Ca, Cb - elastic moduli tensor
% rho(i,j+1),rho(i,j) - densities above and below the interface
% dx,dz - grid step
% eta0x, eta1x -sum=1, show where interface intersects vertical grid line :
% eta0y, eta1y -sum=1, show where interface intersects horizontal grid line ..
% ii - 1 - x, 2- y, to put in Cijkl as i

function [UJ, CJI]=oper_fe(i,j,UI,nvx,nvy,C,rho,dx,dz,eta0x,eta1x,eta0y,eta1y,direct)
    PLT=false;
    xe1ux=1.d0;
    ye1uy=1.d0;
    
    xe2ux_dx=nvy;
    xe2ux_dy=-nvx;
    ye2uy_dx=nvy;
    ye2uy_dy=-nvx;
    
    xe3ux_dxx=nvy^2.d0;
    xe3ux_dyy=nvx^2.d0;
    xe3ux_dxy=2.d0*nvx*nvy;
    ye3uy_dxx=nvy^2.d0;
    ye3uy_dyy=nvx^2.d0;
    ye3uy_dxy=2.d0*nvx*nvy;
    
    if strcmp(direct,'vert')
        %------------------------------------------------------------------
        %                  Vertical
        %------------------------------------------------------------------
        %Above
        xe4ux_dxa=nvx*C(i,j+1,1,1,1,1)+nvy*C(i,j+1,1,2,1,1);
        xe4ux_dya=nvx*C(i,j+1,1,1,1,2)+nvy*C(i,j+1,1,2,1,2);
        xe4uy_dxa=nvx*C(i,j+1,1,1,2,1)+nvy*C(i,j+1,1,2,2,1);
        xe4uy_dya=nvx*C(i,j+1,1,1,2,2)+nvy*C(i,j+1,1,2,2,2);

        ye4ux_dxa=nvx*C(i,j+1,2,1,1,1)+nvy*C(i,j+1,2,2,1,1);
        ye4ux_dya=nvx*C(i,j+1,2,1,1,2)+nvy*C(i,j+1,2,2,1,2);
        ye4uy_dxa=nvx*C(i,j+1,2,1,2,1)+nvy*C(i,j+1,2,2,2,1);
        ye4uy_dya=nvx*C(i,j+1,2,1,2,2)+nvy*C(i,j+1,2,2,2,2);

        %Below
        xe4ux_dxb=nvx*C(i,j,1,1,1,1)+nvy*C(i,j,1,2,1,1);
        xe4ux_dyb=nvx*C(i,j,1,1,1,2)+nvy*C(i,j,1,2,1,2);
        xe4uy_dxb=nvx*C(i,j,1,1,2,1)+nvy*C(i,j,1,2,2,1);
        xe4uy_dyb=nvx*C(i,j,1,1,2,2)+nvy*C(i,j,1,2,2,2);

        ye4ux_dxb=nvx*C(i,j,2,1,1,1)+nvy*C(i,j,2,2,1,1);
        ye4ux_dyb=nvx*C(i,j,2,1,1,2)+nvy*C(i,j,2,2,1,2);
        ye4uy_dxb=nvx*C(i,j,2,1,2,1)+nvy*C(i,j,2,2,2,1);
        ye4uy_dyb=nvx*C(i,j,2,1,2,2)+nvy*C(i,j,2,2,2,2);

        %---------------------------------
        %Above
        xe5ux_dxxa=nvx*nvy*C(i,j+1,1,1,1,1)+nvy^2.d0*C(i,j+1,1,2,1,1);
        xe5ux_dxya=nvx*nvy*C(i,j+1,1,1,1,2)+nvy^2.d0*C(i,j+1,1,2,1,2)-nvx^2.d0*C(i,j+1,1,1,1,1)-nvx*nvy*C(i,j+1,1,2,1,1);
        xe5uy_dxxa=nvx*nvy*C(i,j+1,1,1,2,1)+nvy^2.d0*C(i,j+1,1,2,2,1);
        xe5uy_dxya=nvx*nvy*C(i,j+1,1,1,2,2)+nvy^2.d0*C(i,j+1,1,2,2,2)-nvx^2.d0*C(i,j+1,1,1,2,1)-nvx*nvy*C(i,j+1,1,2,2,1);
        xe5ux_dyya=-1.d0*(nvx*nvy*C(i,j+1,1,2,1,2)+nvx^2.d0*C(i,j+1,1,1,1,2));
        xe5uy_dyya=-1.d0*(nvx*nvy*C(i,j+1,1,2,2,2)+nvx^2.d0*C(i,j+1,1,1,2,2));

        ye5ux_dxxa=nvx*nvy*C(i,j+1,2,1,1,1)+nvy^2.d0*C(i,j+1,2,2,1,1);
        ye5ux_dxya=nvx*nvy*C(i,j+1,2,1,1,2)+nvy^2.d0*C(i,j+1,2,2,1,2)-nvx^2.d0*C(i,j+1,2,1,1,1)-nvx*nvy*C(i,j+1,2,2,1,1);
        ye5uy_dxxa=nvx*nvy*C(i,j+1,2,1,2,1)+nvy^2.d0*C(i,j+1,2,2,2,1);
        ye5uy_dxya=nvx*nvy*C(i,j+1,2,1,2,2)+nvy^2.d0*C(i,j+1,2,2,2,2)-nvx^2.d0*C(i,j+1,2,1,2,1)-nvx*nvy*C(i,j+1,2,2,2,1);
        ye5ux_dyya=-1.d0*(nvx*nvy*C(i,j+1,2,2,1,2)+nvx^2.d0*C(i,j+1,2,1,1,2));
        ye5uy_dyya=-1.d0*(nvx*nvy*C(i,j+1,2,2,2,2)+nvx^2.d0*C(i,j+1,2,1,2,2));

        %Below
        xe5ux_dxxb=nvx*nvy*C(i,j,1,1,1,1)+nvy^2.d0*C(i,j,1,2,1,1);
        xe5uy_dxxb=nvx*nvy*C(i,j,1,1,2,1)+nvy^2.d0*C(i,j,1,2,2,1);
        xe5ux_dxyb=nvx*nvy*C(i,j,1,1,1,2)+nvy^2.d0*C(i,j,1,2,1,2)-nvx^2.d0*C(i,j,1,1,1,1)-nvx*nvy*C(i,j,1,2,1,1);
        xe5uy_dxyb=nvx*nvy*C(i,j,1,1,2,2)+nvy^2.d0*C(i,j,1,2,2,2)-nvx^2.d0*C(i,j,1,1,2,1)-nvx*nvy*C(i,j,1,2,2,1);
        xe5ux_dyyb=-1.d0*(nvx*nvy*C(i,j,1,2,1,2)+nvx^2.d0*C(i,j,1,1,1,2));
        xe5uy_dyyb=-1.d0*(nvx*nvy*C(i,j,1,2,2,2)+nvx^2.d0*C(i,j,1,1,2,2));

        ye5ux_dxxb=nvx*nvy*C(i,j,2,1,1,1)+nvy^2.d0*C(i,j,2,2,1,1);
        ye5ux_dxyb=nvx*nvy*C(i,j,2,1,1,2)+nvy^2.d0*C(i,j,2,2,1,2)-nvx^2.d0*C(i,j,2,1,1,1)-nvx*nvy*C(i,j,2,2,1,1);
        ye5uy_dxxb=nvx*nvy*C(i,j,2,1,2,1)+nvy^2.d0*C(i,j,2,2,2,1);
        ye5uy_dxyb=nvx*nvy*C(i,j,2,1,2,2)+nvy^2.d0*C(i,j,2,2,2,2)-nvx^2.d0*C(i,j,2,1,2,1)-nvx*nvy*C(i,j,2,2,2,1);
        ye5ux_dyyb=-1.d0*(nvx*nvy*C(i,j,2,2,1,2)+nvx^2.d0*C(i,j,2,1,1,2));
        ye5uy_dyyb=-1.d0*(nvx*nvy*C(i,j,2,2,2,2)+nvx^2.d0*C(i,j,2,1,2,2));

        %---------------------------------
        %Above
        xe6ux_dxxa=C(i,j+1,1,1,1,1)/rho(i,j+1);
        xe6ux_dxya=(C(i,j+1,1,1,1,2)+C(i,j+1,1,2,1,1))/rho(i,j+1);
        xe6uy_dxxa=C(i,j+1,1,1,2,1)/rho(i,j+1);
        xe6uy_dxya=(C(i,j+1,1,1,2,2)+C(i,j+1,1,2,2,1))/rho(i,j+1);
        xe6ux_dyya=C(i,j+1,1,2,1,2)/rho(i,j+1);
        xe6uy_dyya=C(i,j+1,1,2,2,2)/rho(i,j+1);

        ye6ux_dxxa=C(i,j+1,2,1,1,1)/rho(i,j+1);
        ye6ux_dxya=(C(i,j+1,2,1,1,2)+C(i,j+1,2,2,1,1))/rho(i,j+1);
        ye6uy_dxxa=C(i,j+1,2,1,2,1)/rho(i,j+1);
        ye6uy_dxya=(C(i,j+1,2,1,2,2)+C(i,j+1,2,2,2,1))/rho(i,j+1);
        ye6ux_dyya=C(i,j+1,2,2,1,2)/rho(i,j+1);
        ye6uy_dyya=C(i,j+1,2,2,2,2)/rho(i,j+1);

        %Below
        xe6ux_dxxb=C(i,j,1,1,1,1)/rho(i,j);
        xe6ux_dxyb=(C(i,j,1,1,1,2)+C(i,j,1,2,1,1))/rho(i,j);
        xe6uy_dxxb=C(i,j,1,1,2,1)/rho(i,j);
        xe6uy_dxyb=(C(i,j,1,1,2,2)+C(i,j,1,2,2,1))/rho(i,j);
        xe6ux_dyyb=C(i,j,1,2,1,2)/rho(i,j);
        xe6uy_dyyb=C(i,j,1,2,2,2)/rho(i,j);

        ye6ux_dxxb=C(i,j,2,1,1,1)/rho(i,j);
        ye6ux_dxyb=(C(i,j,2,1,1,2)+C(i,j,2,2,1,1))/rho(i,j);
        ye6uy_dxxb=C(i,j,2,1,2,1)/rho(i,j);
        ye6uy_dxyb=(C(i,j,2,1,2,2)+C(i,j,2,2,2,1))/rho(i,j);
        ye6ux_dyyb=C(i,j,2,2,1,2)/rho(i,j);
        ye6uy_dyyb=C(i,j,2,2,2,2)/rho(i,j);
    else
        %------------------------------------------------------------------
        %                  Horizontal
        %------------------------------------------------------------------
        %Above
        xe4ux_dxa=nvx*C(i+1,j,1,1,1,1)+nvy*C(i+1,j,1,2,1,1);
        xe4ux_dya=nvx*C(i+1,j,1,1,1,2)+nvy*C(i+1,j,1,2,1,2);
        xe4uy_dxa=nvx*C(i+1,j,1,1,2,1)+nvy*C(i+1,j,1,2,2,1);
        xe4uy_dya=nvx*C(i+1,j,1,1,2,2)+nvy*C(i+1,j,1,2,2,2);

        ye4ux_dxa=nvx*C(i+1,j,2,1,1,1)+nvy*C(i+1,j,2,2,1,1);
        ye4ux_dya=nvx*C(i+1,j,2,1,1,2)+nvy*C(i+1,j,2,2,1,2);
        ye4uy_dxa=nvx*C(i+1,j,2,1,2,1)+nvy*C(i+1,j,2,2,2,1);
        ye4uy_dya=nvx*C(i+1,j,2,1,2,2)+nvy*C(i+1,j,2,2,2,2);

        %Below
        xe4ux_dxb=nvx*C(i,j,1,1,1,1)+nvy*C(i,j,1,2,1,1);
        xe4ux_dyb=nvx*C(i,j,1,1,1,2)+nvy*C(i,j,1,2,1,2);
        xe4uy_dxb=nvx*C(i,j,1,1,2,1)+nvy*C(i,j,1,2,2,1);
        xe4uy_dyb=nvx*C(i,j,1,1,2,2)+nvy*C(i,j,1,2,2,2);

        ye4ux_dxb=nvx*C(i,j,2,1,1,1)+nvy*C(i,j,2,2,1,1);
        ye4ux_dyb=nvx*C(i,j,2,1,1,2)+nvy*C(i,j,2,2,1,2);
        ye4uy_dxb=nvx*C(i,j,2,1,2,1)+nvy*C(i,j,2,2,2,1);
        ye4uy_dyb=nvx*C(i,j,2,1,2,2)+nvy*C(i,j,2,2,2,2);

        %---------------------------------
        %Above
        xe5ux_dxxa=nvx*nvy*C(i+1,j,1,1,1,1)+nvy^2.d0*C(i+1,j,1,2,1,1);
        xe5ux_dxya=nvx*nvy*C(i+1,j,1,1,1,2)+nvy^2.d0*C(i+1,j,1,2,1,2)-nvx^2.d0*C(i+1,j,1,1,1,1)-nvx*nvy*C(i+1,j,1,2,1,1);
        xe5uy_dxxa=nvx*nvy*C(i+1,j,1,1,2,1)+nvy^2.d0*C(i+1,j,1,2,2,1);
        xe5uy_dxya=nvx*nvy*C(i+1,j,1,1,2,2)+nvy^2.d0*C(i+1,j,1,2,2,2)-nvx^2.d0*C(i+1,j,1,1,2,1)-nvx*nvy*C(i+1,j,1,2,2,1);
        xe5ux_dyya=-1.d0*(nvx*nvy*C(i+1,j,1,2,1,2)+nvx^2.d0*C(i+1,j,1,1,1,2));
        xe5uy_dyya=-1.d0*(nvx*nvy*C(i+1,j,1,2,2,2)+nvx^2.d0*C(i+1,j,1,1,2,2));

        ye5ux_dxxa=nvx*nvy*C(i+1,j,2,1,1,1)+nvy^2.d0*C(i+1,j,2,2,1,1);
        ye5ux_dxya=nvx*nvy*C(i+1,j,2,1,1,2)+nvy^2.d0*C(i+1,j,2,2,1,2)-nvx^2.d0*C(i+1,j,2,1,1,1)-nvx*nvy*C(i+1,j,2,2,1,1);
        ye5uy_dxxa=nvx*nvy*C(i+1,j,2,1,2,1)+nvy^2.d0*C(i+1,j,2,2,2,1);
        ye5uy_dxya=nvx*nvy*C(i+1,j,2,1,2,2)+nvy^2.d0*C(i+1,j,2,2,2,2)-nvx^2.d0*C(i+1,j,2,1,2,1)-nvx*nvy*C(i+1,j,2,2,2,1);
        ye5ux_dyya=-1.d0*(nvx*nvy*C(i+1,j,2,2,1,2)+nvx^2.d0*C(i+1,j,2,1,1,2));
        ye5uy_dyya=-1.d0*(nvx*nvy*C(i+1,j,2,2,2,2)+nvx^2.d0*C(i+1,j,2,1,2,2));

        %Below
        xe5ux_dxxb=nvx*nvy*C(i,j,1,1,1,1)+nvy^2.d0*C(i,j,1,2,1,1);
        xe5uy_dxxb=nvx*nvy*C(i,j,1,1,2,1)+nvy^2.d0*C(i,j,1,2,2,1);
        xe5ux_dxyb=nvx*nvy*C(i,j,1,1,1,2)+nvy^2.d0*C(i,j,1,2,1,2)-nvx^2.d0*C(i,j,1,1,1,1)-nvx*nvy*C(i,j,1,2,1,1);
        xe5uy_dxyb=nvx*nvy*C(i,j,1,1,2,2)+nvy^2.d0*C(i,j,1,2,2,2)-nvx^2.d0*C(i,j,1,1,2,1)-nvx*nvy*C(i,j,1,2,2,1);
        xe5ux_dyyb=-1.d0*(nvx*nvy*C(i,j,1,2,1,2)+nvx^2.d0*C(i,j,1,1,1,2));
        xe5uy_dyyb=-1.d0*(nvx*nvy*C(i,j,1,2,2,2)+nvx^2.d0*C(i,j,1,1,2,2));

        ye5ux_dxxb=nvx*nvy*C(i,j,2,1,1,1)+nvy^2.d0*C(i,j,2,2,1,1);
        ye5ux_dxyb=nvx*nvy*C(i,j,2,1,1,2)+nvy^2.d0*C(i,j,2,2,1,2)-nvx^2.d0*C(i,j,2,1,1,1)-nvx*nvy*C(i,j,2,2,1,1);
        ye5uy_dxxb=nvx*nvy*C(i,j,2,1,2,1)+nvy^2.d0*C(i,j,2,2,2,1);
        ye5uy_dxyb=nvx*nvy*C(i,j,2,1,2,2)+nvy^2.d0*C(i,j,2,2,2,2)-nvx^2.d0*C(i,j,2,1,2,1)-nvx*nvy*C(i,j,2,2,2,1);
        ye5ux_dyyb=-1.d0*(nvx*nvy*C(i,j,2,2,1,2)+nvx^2.d0*C(i,j,2,1,1,2));
        ye5uy_dyyb=-1.d0*(nvx*nvy*C(i,j,2,2,2,2)+nvx^2.d0*C(i,j,2,1,2,2));

        %---------------------------------
        %Above
        xe6ux_dxxa=C(i+1,j,1,1,1,1)/rho(i+1,j);
        xe6ux_dxya=(C(i+1,j,1,1,1,2)+C(i+1,j,1,2,1,1))/rho(i+1,j);
        xe6uy_dxxa=C(i+1,j,1,1,2,1)/rho(i+1,j);
        xe6uy_dxya=(C(i+1,j,1,1,2,2)+C(i+1,j,1,2,2,1))/rho(i+1,j);
        xe6ux_dyya=C(i+1,j,1,2,1,2)/rho(i+1,j);
        xe6uy_dyya=C(i+1,j,1,2,2,2)/rho(i+1,j);

        ye6ux_dxxa=C(i+1,j,2,1,1,1)/rho(i+1,j);
        ye6ux_dxya=(C(i+1,j,2,1,1,2)+C(i+1,j,2,2,1,1))/rho(i+1,j);
        ye6uy_dxxa=C(i+1,j,2,1,2,1)/rho(i+1,j);
        ye6uy_dxya=(C(i+1,j,2,1,2,2)+C(i+1,j,2,2,2,1))/rho(i+1,j);
        ye6ux_dyya=C(i+1,j,2,2,1,2)/rho(i+1,j);
        ye6uy_dyya=C(i+1,j,2,2,2,2)/rho(i+1,j);

        %Below
        xe6ux_dxxb=C(i,j,1,1,1,1)/rho(i,j);
        xe6ux_dxyb=(C(i,j,1,1,1,2)+C(i,j,1,2,1,1))/rho(i,j);
        xe6uy_dxxb=C(i,j,1,1,2,1)/rho(i,j);
        xe6uy_dxyb=(C(i,j,1,1,2,2)+C(i,j,1,2,2,1))/rho(i,j);
        xe6ux_dyyb=C(i,j,1,2,1,2)/rho(i,j);
        xe6uy_dyyb=C(i,j,1,2,2,2)/rho(i,j);

        ye6ux_dxxb=C(i,j,2,1,1,1)/rho(i,j);
        ye6ux_dxyb=(C(i,j,2,1,1,2)+C(i,j,2,2,1,1))/rho(i,j);
        ye6uy_dxxb=C(i,j,2,1,2,1)/rho(i,j);
        ye6uy_dxyb=(C(i,j,2,1,2,2)+C(i,j,2,2,2,1))/rho(i,j);
        ye6ux_dyyb=C(i,j,2,2,1,2)/rho(i,j);
        ye6uy_dyyb=C(i,j,2,2,2,2)/rho(i,j);
    end
    
    Bb=[xe1ux 0 0 0 0 0 0 0 0 0 0 0; ...
        0 xe2ux_dx xe2ux_dy 0 0 0 0 0 0 0 0 0; ...
        0 0 0 xe3ux_dxx xe3ux_dyy xe3ux_dxy 0 0 0 0 0 0; ...
        0 xe4ux_dxb xe4ux_dyb 0 0 0 0 xe4uy_dxb xe4uy_dyb 0 0 0; ...
        0 0 0 xe5ux_dxxb xe5ux_dyyb xe5ux_dxyb 0 0 0 xe5uy_dxxb xe5uy_dyyb xe5uy_dxyb; ...
        0 0 0 xe6ux_dxxb xe6ux_dyyb xe6ux_dxyb 0 0 0 xe6uy_dxxb xe6uy_dyyb xe6uy_dxyb; ...
        0 0 0 0 0 0 ye1uy 0 0 0 0 0; ...
        0 0 0 0 0 0 0 ye2uy_dx ye2uy_dy 0 0 0; ...
        0 0 0 0 0 0 0 0 0 ye3uy_dxx ye3uy_dyy ye3uy_dxy; ...
        0 ye4ux_dxb ye4ux_dyb 0 0 0 0 ye4uy_dxb ye4uy_dyb 0 0 0; ...
        0 0 0 ye5ux_dxxb ye5ux_dyyb ye5ux_dxyb 0 0 0 ye5uy_dxxb ye5uy_dyyb ye5uy_dxyb; ...
        0 0 0 ye6ux_dxxb ye6ux_dyyb ye6ux_dxyb 0 0 0 ye6uy_dxxb ye6uy_dyyb ye6uy_dxyb];

    r=eta0x;
    s=eta0y;

    Ab=[1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz 0 0 0 0 0 0; ...
        0 1.d0 0 r*dx 0 s*dz 0 0 0 0 0 0; ...
        0 0 1.d0 0 s*dz r*dx 0 0 0 0 0 0; ...
        0 0 0 1.d0 0 0 0 0 0 0 0 0; ...
        0 0 0 0 1.d0 0 0 0 0 0 0 0; ...
        0 0 0 0 0 1.d0 0 0 0 0 0 0;
        0 0 0 0 0 0 1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz;...
        0 0 0 0 0 0 0 1.d0 0 r*dx 0 s*dz;...
        0 0 0 0 0 0 0 0 1.d0 0 s*dz r*dx;...
        0 0 0 0 0 0 0 0 0 1.d0 0 0;...
        0 0 0 0 0 0 0 0 0 0 1.d0 0;...
        0 0 0 0 0 0 0 0 0 0 0 1.d0];
    
    if PLT 
        subplot(2,3,3);
        pcolor(Ab);
        colorbar();
        title(['Ab:' 'r=' num2str(r) ' s=' num2str(s)]);
        set(gca,'YDir','reverse');
    end

    Ba=[xe1ux 0 0 0 0 0 0 0 0 0 0 0; ...
        0 xe2ux_dx xe2ux_dy 0 0 0 0 0 0 0 0 0; ...
        0 0 0 xe3ux_dxx xe3ux_dyy xe3ux_dxy 0 0 0 0 0 0; ...
        0 xe4ux_dxa xe4ux_dya 0 0 0 0 xe4uy_dxa xe4uy_dya 0 0 0; ...
        0 0 0 xe5ux_dxxa xe5ux_dyya xe5ux_dxya 0 0 0 xe5uy_dxxa xe5uy_dyya xe5uy_dxya; ...
        0 0 0 xe6ux_dxxa xe6ux_dyya xe6ux_dxya 0 0 0 xe6uy_dxxa xe6uy_dyya xe6uy_dxya; ...
        0 0 0 0 0 0 ye1uy 0 0 0 0 0; ...
        0 0 0 0 0 0 0 ye2uy_dx ye2uy_dy 0 0 0; ...
        0 0 0 0 0 0 0 0 0 ye3uy_dxx ye3uy_dyy ye3uy_dxy; ...
        0 ye4ux_dxa ye4ux_dya 0 0 0 0 ye4uy_dxa ye4uy_dya 0 0 0; ...
        0 0 0 ye5ux_dxxa ye5ux_dyya ye5ux_dxya 0 0 0 ye5uy_dxxa ye5uy_dyya ye5uy_dxya; ...
        0 0 0 ye6ux_dxxa ye6ux_dyya ye6ux_dxya 0 0 0 ye6uy_dxxa ye6uy_dyya ye6uy_dxya];

    r=eta1x;
    s=eta1y;
    Aa=[1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz 0 0 0 0 0 0; ...
        0 1.d0 0 r*dx 0 s*dz 0 0 0 0 0 0; ...
        0 0 1.d0 0 s*dz r*dx 0 0 0 0 0 0; ...
        0 0 0 1.d0 0 0 0 0 0 0 0 0; ...
        0 0 0 0 1.d0 0 0 0 0 0 0 0; ...
        0 0 0 0 0 1.d0 0 0 0 0 0 0;
        0 0 0 0 0 0 1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz;...
        0 0 0 0 0 0 0 1.d0 0 r*dx 0 s*dz;...
        0 0 0 0 0 0 0 0 1.d0 0 s*dz r*dx;...
        0 0 0 0 0 0 0 0 0 1.d0 0 0;...
        0 0 0 0 0 0 0 0 0 0 1.d0 0;...
        0 0 0 0 0 0 0 0 0 0 0 1.d0];
    
    if PLT
        subplot(2,3,1);
        pcolor(Aa)
        colorbar();
        title(['Aa:' 'r=' num2str(r) ' s=' num2str(s)]);
        set(gca,'YDir','reverse');
    end
    
%     Ba=eye(12,12);
%     Bb=eye(12,12);
    
    CJI=svdinv(Ba*Aa)*Bb*Ab;
    %UJ=svdinv(Ba*Aa)*(Bb*Ab)*UI;
    UJ=CJI*UI;
    
    if PLT
        subplot(2,3,2);
        pcolor(Ba);
        colorbar();
        title('Ba');
        set(gca,'YDir','reverse');
        
        subplot(2,3,4);
        pcolor(Bb);
        colorbar();
        title('Bb');
        set(gca,'YDir','reverse');
        %set(gca,'XDir','reverse');

        subplot(2,3,5);
        pcolor(CJI);
        colorbar();
        title('CJI');
        set(gca,'YDir','reverse');
        %set(gca,'XDir','reverse');
        
        subplot(2,3,6);
%         pcolor(CJI-(svdinv(Aa)*Ab));
          pcolor(svdinv(Ba*eye(12))*Bb*eye(12));
         colorbar();
        title('Ba-1 * Bb');
        set(gca,'YDir','reverse');

        set(gcf,'Name',['xEta0=' num2str(eta0x) ' yeta0=' num2str(eta0y)]);
        set(gcf,'NumberTitle','off');
        drawnow;
        %input('Next?');
    end
end