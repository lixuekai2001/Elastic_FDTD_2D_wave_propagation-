% . ---|--> .

function [A0,B0,A1,B1]=A0B0A1B1(i,j,iia,jja,iib,jjb,nvx,nvy,C,rho,dx,dz,eta0x,eta1x,eta0y,eta1y)

%BC eq.1
    xe1ux=1.d0;
    ye1uy=1.d0;

%BC eq.2
    xe2ux_dx=nvy;
    xe2ux_dy=-nvx;
    ye2uy_dx=nvy;
    ye2uy_dy=-nvx;

%BC eq.3
    xe3ux_dxx=nvy^2.d0;
    xe3ux_dyy=nvx^2.d0;
    xe3ux_dxy=2.d0*nvx*nvy;
    ye3uy_dxx=nvy^2.d0;
    ye3uy_dyy=nvx^2.d0;
    ye3uy_dxy=2.d0*nvx*nvy;
    
%BC eq.4
    %Above
    xe4ux_dxa=nvx*C(i+iia,j+jja,1,1,1,1)+nvy*C(i+iia,j+jja,1,2,1,1);
    xe4ux_dya=nvx*C(i+iia,j+jja,1,1,1,2)+nvy*C(i+iia,j+jja,1,2,1,2);
    xe4uy_dxa=nvx*C(i+iia,j+jja,1,1,2,1)+nvy*C(i+iia,j+jja,1,2,2,1);
    xe4uy_dya=nvx*C(i+iia,j+jja,1,1,2,2)+nvy*C(i+iia,j+jja,1,2,2,2);

    ye4ux_dxa=nvx*C(i+iia,j+jja,2,1,1,1)+nvy*C(i+iia,j+jja,2,2,1,1);
    ye4ux_dya=nvx*C(i+iia,j+jja,2,1,1,2)+nvy*C(i+iia,j+jja,2,2,1,2);
    ye4uy_dxa=nvx*C(i+iia,j+jja,2,1,2,1)+nvy*C(i+iia,j+jja,2,2,2,1);
    ye4uy_dya=nvx*C(i+iia,j+jja,2,1,2,2)+nvy*C(i+iia,j+jja,2,2,2,2);

    %Below
    xe4ux_dxb=nvx*C(i+iib,j+jjb,1,1,1,1)+nvy*C(i+iib,j+jjb,1,2,1,1);
    xe4ux_dyb=nvx*C(i+iib,j+jjb,1,1,1,2)+nvy*C(i+iib,j+jjb,1,2,1,2);
    xe4uy_dxb=nvx*C(i+iib,j+jjb,1,1,2,1)+nvy*C(i+iib,j+jjb,1,2,2,1);
    xe4uy_dyb=nvx*C(i+iib,j+jjb,1,1,2,2)+nvy*C(i+iib,j+jjb,1,2,2,2);

    ye4ux_dxb=nvx*C(i+iib,j+jjb,2,1,1,1)+nvy*C(i+iib,j+jjb,2,2,1,1);
    ye4ux_dyb=nvx*C(i+iib,j+jjb,2,1,1,2)+nvy*C(i+iib,j+jjb,2,2,1,2);
    ye4uy_dxb=nvx*C(i+iib,j+jjb,2,1,2,1)+nvy*C(i+iib,j+jjb,2,2,2,1);
    ye4uy_dyb=nvx*C(i+iib,j+jjb,2,1,2,2)+nvy*C(i+iib,j+jjb,2,2,2,2);

%BC eq.5
    %Above
    xe5ux_dxxa=nvx*nvy*C(i+iia,j+jja,1,1,1,1)+nvy^2.d0*C(i+iia,j+jja,1,2,1,1);
    xe5ux_dxya=nvx*nvy*C(i+iia,j+jja,1,1,1,2)+nvy^2.d0*C(i+iia,j+jja,1,2,1,2)-nvx^2.d0*C(i+iia,j+jja,1,1,1,1)-nvx*nvy*C(i+iia,j+jja,1,2,1,1);
    xe5uy_dxxa=nvx*nvy*C(i+iia,j+jja,1,1,2,1)+nvy^2.d0*C(i+iia,j+jja,1,2,2,1);
    xe5uy_dxya=nvx*nvy*C(i+iia,j+jja,1,1,2,2)+nvy^2.d0*C(i+iia,j+jja,1,2,2,2)-nvx^2.d0*C(i+iia,j+jja,1,1,2,1)-nvx*nvy*C(i+iia,j+jja,1,2,2,1);
    xe5ux_dyya=-1.d0*(nvx*nvy*C(i+iia,j+jja,1,2,1,2)+nvx^2.d0*C(i+iia,j+jja,1,1,1,2));
    xe5uy_dyya=-1.d0*(nvx*nvy*C(i+iia,j+jja,1,2,2,2)+nvx^2.d0*C(i+iia,j+jja,1,1,2,2));

    ye5ux_dxxa=nvx*nvy*C(i+iia,j+jja,2,1,1,1)+nvy^2.d0*C(i+iia,j+jja,2,2,1,1);
    ye5ux_dxya=nvx*nvy*C(i+iia,j+jja,2,1,1,2)+nvy^2.d0*C(i+iia,j+jja,2,2,1,2)-nvx^2.d0*C(i+iia,j+jja,2,1,1,1)-nvx*nvy*C(i+iia,j+jja,2,2,1,1);
    ye5uy_dxxa=nvx*nvy*C(i+iia,j+jja,2,1,2,1)+nvy^2.d0*C(i+iia,j+jja,2,2,2,1);
    ye5uy_dxya=nvx*nvy*C(i+iia,j+jja,2,1,2,2)+nvy^2.d0*C(i+iia,j+jja,2,2,2,2)-nvx^2.d0*C(i+iia,j+jja,2,1,2,1)-nvx*nvy*C(i+iia,j+jja,2,2,2,1);
    ye5ux_dyya=-1.d0*(nvx*nvy*C(i+iia,j+jja,2,2,1,2)+nvx^2.d0*C(i+iia,j+jja,2,1,1,2));
    ye5uy_dyya=-1.d0*(nvx*nvy*C(i+iia,j+jja,2,2,2,2)+nvx^2.d0*C(i+iia,j+jja,2,1,2,2));

    %Below
    xe5ux_dxxb=nvx*nvy*C(i+iib,j+jjb,1,1,1,1)+nvy^2.d0*C(i+iib,j+jjb,1,2,1,1);
    xe5uy_dxxb=nvx*nvy*C(i+iib,j+jjb,1,1,2,1)+nvy^2.d0*C(i+iib,j+jjb,1,2,2,1);
    xe5ux_dxyb=nvx*nvy*C(i+iib,j+jjb,1,1,1,2)+nvy^2.d0*C(i+iib,j+jjb,1,2,1,2)-nvx^2.d0*C(i+iib,j+jjb,1,1,1,1)-nvx*nvy*C(i+iib,j+jjb,1,2,1,1);
    xe5uy_dxyb=nvx*nvy*C(i+iib,j+jjb,1,1,2,2)+nvy^2.d0*C(i+iib,j+jjb,1,2,2,2)-nvx^2.d0*C(i+iib,j+jjb,1,1,2,1)-nvx*nvy*C(i+iib,j+jjb,1,2,2,1);
    xe5ux_dyyb=-1.d0*(nvx*nvy*C(i+iib,j+jjb,1,2,1,2)+nvx^2.d0*C(i+iib,j+jjb,1,1,1,2));
    xe5uy_dyyb=-1.d0*(nvx*nvy*C(i+iib,j+jjb,1,2,2,2)+nvx^2.d0*C(i+iib,j+jjb,1,1,2,2));

    ye5ux_dxxb=nvx*nvy*C(i+iib,j+jjb,2,1,1,1)+nvy^2.d0*C(i+iib,j+jjb,2,2,1,1);
    ye5ux_dxyb=nvx*nvy*C(i+iib,j+jjb,2,1,1,2)+nvy^2.d0*C(i+iib,j+jjb,2,2,1,2)-nvx^2.d0*C(i+iib,j+jjb,2,1,1,1)-nvx*nvy*C(i+iib,j+jjb,2,2,1,1);
    ye5uy_dxxb=nvx*nvy*C(i+iib,j+jjb,2,1,2,1)+nvy^2.d0*C(i+iib,j+jjb,2,2,2,1);
    ye5uy_dxyb=nvx*nvy*C(i+iib,j+jjb,2,1,2,2)+nvy^2.d0*C(i+iib,j+jjb,2,2,2,2)-nvx^2.d0*C(i+iib,j+jjb,2,1,2,1)-nvx*nvy*C(i+iib,j+jjb,2,2,2,1);
    ye5ux_dyyb=-1.d0*(nvx*nvy*C(i+iib,j+jjb,2,2,1,2)+nvx^2.d0*C(i+iib,j+jjb,2,1,1,2));
    ye5uy_dyyb=-1.d0*(nvx*nvy*C(i+iib,j+jjb,2,2,2,2)+nvx^2.d0*C(i+iib,j+jjb,2,1,2,2));

%BC eq.6
    %Above
    xe6ux_dxxa=C(i+iia,j+jja,1,1,1,1)/rho(i+iia,j+jja);
    xe6ux_dxya=(C(i+iia,j+jja,1,1,1,2)+C(i+iia,j+jja,1,2,1,1))/rho(i+iia,j+jja);
    xe6uy_dxxa=C(i+iia,j+jja,1,1,2,1)/rho(i+iia,j+jja);
    xe6uy_dxya=(C(i+iia,j+jja,1,1,2,2)+C(i+iia,j+jja,1,2,2,1))/rho(i+iia,j+jja);
    xe6ux_dyya=C(i+iia,j+jja,1,2,1,2)/rho(i+iia,j+jja);
    xe6uy_dyya=C(i+iia,j+jja,1,2,2,2)/rho(i+iia,j+jja);

    ye6ux_dxxa=C(i+iia,j+jja,2,1,1,1)/rho(i+iia,j+jja);
    ye6ux_dxya=(C(i+iia,j+jja,2,1,1,2)+C(i+iia,j+jja,2,2,1,1))/rho(i+iia,j+jja);
    ye6uy_dxxa=C(i+iia,j+jja,2,1,2,1)/rho(i+iia,j+jja);
    ye6uy_dxya=(C(i+iia,j+jja,2,1,2,2)+C(i+iia,j+jja,2,2,2,1))/rho(i+iia,j+jja);
    ye6ux_dyya=C(i+iia,j+jja,2,2,1,2)/rho(i+iia,j+jja);
    ye6uy_dyya=C(i+iia,j+jja,2,2,2,2)/rho(i+iia,j+jja);

    %Below
    xe6ux_dxxb=C(i+iib,j+jjb,1,1,1,1)/rho(i+iib,j+jjb);
    xe6ux_dxyb=(C(i+iib,j+jjb,1,1,1,2)+C(i+iib,j+jjb,1,2,1,1))/rho(i+iib,j+jjb);
    xe6uy_dxxb=C(i+iib,j+jjb,1,1,2,1)/rho(i+iib,j+jjb);
    xe6uy_dxyb=(C(i+iib,j+jjb,1,1,2,2)+C(i+iib,j+jjb,1,2,2,1))/rho(i+iib,j+jjb);
    xe6ux_dyyb=C(i+iib,j+jjb,1,2,1,2)/rho(i+iib,j+jjb);
    xe6uy_dyyb=C(i+iib,j+jjb,1,2,2,2)/rho(i+iib,j+jjb);

    ye6ux_dxxb=C(i+iib,j+jjb,2,1,1,1)/rho(i+iib,j+jjb);
    ye6ux_dxyb=(C(i+iib,j+jjb,2,1,1,2)+C(i+iib,j+jjb,2,2,1,1))/rho(i+iib,j+jjb);
    ye6uy_dxxb=C(i+iib,j+jjb,2,1,2,1)/rho(i+iib,j+jjb);
    ye6uy_dxyb=(C(i+iib,j+jjb,2,1,2,2)+C(i+iib,j+jjb,2,2,2,1))/rho(i+iib,j+jjb);
    ye6ux_dyyb=C(i+iib,j+jjb,2,2,1,2)/rho(i+iib,j+jjb);
    ye6uy_dyyb=C(i+iib,j+jjb,2,2,2,2)/rho(i+iib,j+jjb);
    
    
%Full matrices for left and right points
    B0=[xe1ux 0 0 0 0 0 0 0 0 0 0 0; ...
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
    A0=[1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz 0 0 0 0 0 0; ...
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
    
    B1=[xe1ux 0 0 0 0 0 0 0 0 0 0 0; ...
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

    r=-eta1x;
    s=-eta1y;
    A1=[1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz 0 0 0 0 0 0; ...
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
end