% . ---|--> .
%v2
function [A0,B0,A1,B1] = A0B0A1B1triso2(i,j,iia,jja,iib,jjb,nvx,nvy,C,rho,dx,dz,eta0x,eta1x,eta0y,eta1y,DEBUG_MODE)


nvy2 = nvy^2.d0;
nvx2 = nvx^2.d0;
nvxnvy = nvx*nvy;

if nargin > 16
    SHOWALL=DEBUG_MODE;
    if SHOWALL
        fprintf('\n iia = %d\n jja = %d\n iib = %d\n jjb = %d\n',iia,jja,iib,jjb);
        fprintf(' eta0x = %f\n eta1x = %f\n eta0y = %f\n eta1y = %f\n',eta0x,eta1x,eta0y,eta1y);
        fprintf(' nvx=%f\n nvy=%f\n',nvx,nvy);
        fprintf('\n');
    end
else
    SHOWALL = 0;

%BC eq.1
    xe1ux = 1.d0;
    ye1uy = 1.d0;
   
%BC eq.2
    xe2ux_dx = nvy;
    xe2ux_dy = -nvx;
    ye2uy_dx = nvy;
    ye2uy_dy = -nvx;

%BC eq.3
    xe3ux_dxx = nvy2;
    xe3ux_dyy = nvx2;
    xe3ux_dxy = -2.d0*nvxnvy;
    ye3uy_dxx = nvy2;
    ye3uy_dyy = nvx2;
    ye3uy_dxy = -2.d0*nvxnvy;

%BC eq.4
    %Above
    xe4ux_dxa = nvx*C(i+iia,j+jja,1);
    xe4ux_dya = nvy*C(i+iia,j+jja,4);
    xe4uy_dxa = nvy*C(i+iia,j+jja,4);
    xe4uy_dya = nvx*C(i+iia,j+jja,2);
    
    ye4ux_dxa = nvy*C(i+iia,j+jja,2);
    ye4ux_dya = nvx*C(i+iia,j+jja,4);
    ye4uy_dxa = nvx*C(i+iia,j+jja,4);
    ye4uy_dya = nvy*C(i+iia,j+jja,3);

    %Below
    xe4ux_dxb = nvx*C(i+iib,j+jjb,1);
    xe4ux_dyb = nvy*C(i+iib,j+jjb,4);
    xe4uy_dxb = nvy*C(i+iib,j+jjb,4);
    xe4uy_dyb = nvx*C(i+iib,j+jjb,2);

    ye4ux_dxb = nvy*C(i+iib,j+jjb,2);
    ye4ux_dyb = nvx*C(i+iib,j+jjb,4);
    ye4uy_dxb = nvx*C(i+iib,j+jjb,4);
    ye4uy_dyb = nvy*C(i+iib,j+jjb,3);

%BC eq.5
    %Above
    xe5ux_dxxa = nvxnvy*C(i+iia,j+jja,1);
    xe5ux_dxya = nvy2*C(i+iia,j+jja,4)-nvx2*C(i+iia,j+jja,1);
    xe5uy_dxxa = nvy2*C(i+iia,j+jja,4);
    xe5uy_dxya = nvxnvy*C(i+iia,j+jja,2)-nvxnvy*C(i+iia,j+jja,4);
    xe5ux_dyya = - nvxnvy*C(i+iia,j+jja,4);
    xe5uy_dyya = - nvx2*C(i+iia,j+jja,2);

    ye5ux_dxxa = nvy2*C(i+iia,j+jja,2);
    ye5ux_dxya = nvxnvy*C(i+iia,j+jja,4)-nvxnvy*C(i+iia,j+jja,2);
    ye5uy_dxxa = nvxnvy*C(i+iia,j+jja,4);
    ye5uy_dxya = nvy2*C(i+iia,j+jja,3)-nvx2*C(i+iia,j+jja,4);
    ye5ux_dyya =  - nvx2*C(i+iia,j+jja,4);
    ye5uy_dyya =  - nvxnvy*C(i+iia,j+jja,3);

    %Below
    xe5ux_dxxb = nvxnvy*C(i+iib,j+jjb,1);
    xe5ux_dxyb = nvy2*C(i+iib,j+jjb,4)-nvx2*C(i+iib,j+jjb,1);
    xe5uy_dxxb = nvy2*C(i+iib,j+jjb,4);
    xe5uy_dxyb = nvxnvy*C(i+iib,j+jjb,2)-nvxnvy*C(i+iib,j+jjb,4);
    xe5ux_dyyb = - nvxnvy*C(i+iib,j+jjb,4);
    xe5uy_dyyb = - nvx2*C(i+iib,j+jjb,2);

    ye5ux_dxxb = nvy2*C(i+iib,j+jjb,2);
    ye5ux_dxyb = nvxnvy*C(i+iib,j+jjb,4)-nvxnvy*C(i+iib,j+jjb,2);
    ye5uy_dxxb = nvxnvy*C(i+iib,j+jjb,4);
    ye5uy_dxyb = nvy2*C(i+iib,j+jjb,3)-nvx2*C(i+iib,j+jjb,4);
    ye5ux_dyyb =  - nvx2*C(i+iib,j+jjb,4);
    ye5uy_dyyb =  - nvxnvy*C(i+iib,j+jjb,3);

%BC eq.6
    %Above
    xe6ux_dxxa = C(i+iia,j+jja,1)/rho(i+iia,j+jja);
    xe6uy_dxya = (C(i+iia,j+jja,2)+C(i+iia,j+jja,4))/rho(i+iia,j+jja);
    xe6ux_dyya = C(i+iia,j+jja,4)/rho(i+iia,j+jja);

    ye6ux_dxya = (C(i+iia,j+jja,4)+C(i+iia,j+jja,2))/rho(i+iia,j+jja);
    ye6uy_dxxa = C(i+iia,j+jja,4)/rho(i+iia,j+jja);
    ye6uy_dyya = C(i+iia,j+jja,3)/rho(i+iia,j+jja);

    %Below
    xe6ux_dxxb = C(i+iib,j+jjb,1)/rho(i+iib,j+jjb);
    xe6uy_dxyb = (C(i+iib,j+jjb,2)+C(i+iib,j+jjb,4))/rho(i+iib,j+jjb);
    xe6ux_dyyb = C(i+iib,j+jjb,4)/rho(i+iib,j+jjb);

    ye6ux_dxyb = (C(i+iib,j+jjb,4)+C(i+iib,j+jjb,2))/rho(i+iib,j+jjb);
    ye6uy_dxxb = C(i+iib,j+jjb,4)/rho(i+iib,j+jjb);
    ye6uy_dyyb = C(i+iib,j+jjb,3)/rho(i+iib,j+jjb);
    
    
%Full matrices for left and right points
    B0 = [xe1ux 0 0 0 0 0 0 0 0 0 0 0; ...
        0 xe2ux_dx xe2ux_dy 0 0 0 0 0 0 0 0 0; ...
        0 0 0 xe3ux_dxx xe3ux_dyy xe3ux_dxy 0 0 0 0 0 0; ...
        0 xe4ux_dxb xe4ux_dyb 0 0 0 0 xe4uy_dxb xe4uy_dyb 0 0 0; ...
        0 0 0 xe5ux_dxxb xe5ux_dyyb xe5ux_dxyb 0 0 0 xe5uy_dxxb xe5uy_dyyb xe5uy_dxyb; ...
        0 0 0 xe6ux_dxxb xe6ux_dyyb 0 0 0 0 0 0 xe6uy_dxyb; ...
        0 0 0 0 0 0 ye1uy 0 0 0 0 0; ...
        0 0 0 0 0 0 0 ye2uy_dx ye2uy_dy 0 0 0; ...
        0 0 0 0 0 0 0 0 0 ye3uy_dxx ye3uy_dyy ye3uy_dxy; ...
        0 ye4ux_dxb ye4ux_dyb 0 0 0 0 ye4uy_dxb ye4uy_dyb 0 0 0; ...
        0 0 0 ye5ux_dxxb ye5ux_dyyb ye5ux_dxyb 0 0 0 ye5uy_dxxb ye5uy_dyyb ye5uy_dxyb; ...
        0 0 0 0 0 ye6ux_dxyb 0 0 0 ye6uy_dxxb ye6uy_dyyb 0];    
    if SHOWALL
        subplot(233);
        pcolor(flipud(B0));
        colorbar();
        title(['B0 iib=' num2str(iib) ' jjb=' num2str(jjb)]);
    end

    r = eta0x;
    s = eta0y;
    
        A0=[1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz 0 0 0 0 0 0; ...
            0 1.d0 0 r*dx 0 s*dz 0 0 0 0 0 0; ...
            0 0 1.d0 0 s*dz r*dx 0 0 0 0 0 0; ...
            0 0 0 1.d0 0 0 0 0 0 0 0 0; ...
            0 0 0 0 1.d0 0 0 0 0 0 0 0; ...
            0 0 0 0 0 1.d0 0 0 0 0 0 0; ...
            0 0 0 0 0 0 1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz; ...
            0 0 0 0 0 0 0 1.d0 0 r*dx 0 s*dz; ...
            0 0 0 0 0 0 0 0 1.d0 0 s*dz r*dx; ...
            0 0 0 0 0 0 0 0 0 1.d0 0 0; ...
            0 0 0 0 0 0 0 0 0 0 1.d0 0; ...
            0 0 0 0 0 0 0 0 0 0 0 1.d0];
        
    if SHOWALL
        subplot(232);
        pcolor(flipud(A0));
        colorbar();
        title(['A0 r=' num2str(r) ' s=' num2str(s)]);
    end
    
    B1 = [xe1ux 0 0 0 0 0 0 0 0 0 0 0; ...
        0 xe2ux_dx xe2ux_dy 0 0 0 0 0 0 0 0 0; ...
        0 0 0 xe3ux_dxx xe3ux_dyy xe3ux_dxy 0 0 0 0 0 0; ...
        0 xe4ux_dxa xe4ux_dya 0 0 0 0 xe4uy_dxa xe4uy_dya 0 0 0; ...
        0 0 0 xe5ux_dxxa xe5ux_dyya xe5ux_dxya 0 0 0 xe5uy_dxxa xe5uy_dyya xe5uy_dxya; ...
        0 0 0 xe6ux_dxxa xe6ux_dyya 0 0 0 0 0 0 xe6uy_dxya; ...
        0 0 0 0 0 0 ye1uy 0 0 0 0 0; ...
        0 0 0 0 0 0 0 ye2uy_dx ye2uy_dy 0 0 0; ...
        0 0 0 0 0 0 0 0 0 ye3uy_dxx ye3uy_dyy ye3uy_dxy; ...
        0 ye4ux_dxa ye4ux_dya 0 0 0 0 ye4uy_dxa ye4uy_dya 0 0 0; ...
        0 0 0 ye5ux_dxxa ye5ux_dyya ye5ux_dxya 0 0 0 ye5uy_dxxa ye5uy_dyya ye5uy_dxya; ...
        0 0 0 0 0 ye6ux_dxya 0 0 0 ye6uy_dxxa ye6uy_dyya 0];    
    if SHOWALL
        subplot(236);
        pcolor(flipud(B1));
        colorbar();
        title(['B1 iia=' num2str(iia) ' jja=' num2str(jja)]);
    end

    r = eta1x;
    s = eta1y;
    
        A1=[1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz 0 0 0 0 0 0; ...
            0 1.d0 0 r*dx 0 s*dz 0 0 0 0 0 0; ...
            0 0 1.d0 0 s*dz r*dx 0 0 0 0 0 0; ...
            0 0 0 1.d0 0 0 0 0 0 0 0 0; ...
            0 0 0 0 1.d0 0 0 0 0 0 0 0; ...
            0 0 0 0 0 1.d0 0 0 0 0 0 0; ...
            0 0 0 0 0 0 1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz; ...
            0 0 0 0 0 0 0 1.d0 0 r*dx 0 s*dz; ...
            0 0 0 0 0 0 0 0 1.d0 0 s*dz r*dx; ...
            0 0 0 0 0 0 0 0 0 1.d0 0 0; ...
            0 0 0 0 0 0 0 0 0 0 1.d0 0; ...
            0 0 0 0 0 0 0 0 0 0 0 1.d0];
    if SHOWALL
        subplot(235);
        pcolor(flipud(A1));
        colorbar();
        title(['A1 r=' num2str(r) ' s=' num2str(s)]);
    end


    
%     if SHOWALL
% fprintf('eq1:\n');
%     fprintf(' xe1ux = %.2f\n ye1uy = %.2f\n\n',xe1ux,ye1uy);
% fprintf('eq2:\n'); 
%     fprintf(' xe2ux_dx = nvy = %.4f\n xe2ux_dy = - nvx = %.4f\n ye2uy_dx = nvy = %.4f\n ye2uy_dy = -nvx = %.4f\n\n',xe2ux_dx,xe2ux_dy,ye2uy_dx,ye2uy_dy);
% fprintf('eq3:\n');
%     fprintf(' xe3ux_dxx = nvy2 = %f\n xe3ux_dyy = nvx2 = %f\n xe3ux_dxy = -2.d0*nvxnv = %f\n',xe3ux_dxx,xe3ux_dyy,xe3ux_dxy);
%     fprintf(' ye3uy_dxx = nvy2 = %f\n ye3uy_dyy = nvx2 = %f\n ye3uy_dxy = -2.d0*nvxnv = %f\n\n',ye3uy_dxx,ye3uy_dyy,ye3uy_dxy);
%     
% fprintf('eq4:\n');
%     fprintf('B1:\n');
%     fprintf('xe4ux_dxa = nvx*C(i+iia,j+jja,1) = %e\n', xe4ux_dxa);
%     fprintf('xe4ux_dya = nvy*C(i+iia,j+jja,4) = %e\n',xe4ux_dya);
%     fprintf('xe4uy_dxa = nvy*C(i+iia,j+jja,4) = %e\n',xe4uy_dxa);
%     fprintf('xe4uy_dya = nvx*C(i+iia,j+jja,2) = %e\n\n',xe4uy_dya);
%     fprintf('ye4ux_dxa = nvx*C(i+iia,j+jja,2) = %e\n',ye4ux_dxa);
%     fprintf('ye4ux_dya = nvy*C(i+iia,j+jja,4) = %e\n',ye4ux_dya);
%     fprintf('ye4uy_dxa = nvy*C(i+iia,j+jja,4) = %e\n',ye4uy_dxa);
%     fprintf('ye4uy_dya = nvx*C(i+iia,j+jja,3) = %e\n\n',ye4uy_dya);
%     
%     fprintf('B0:\n');
%     fprintf('xe4ux_dxb = nvx*C(i+iib,j+jjb,1) = %e\n',xe4ux_dxb);
%     fprintf('xe4ux_dyb = nvy*C(i+iib,j+jjb,4) = %e\n',xe4ux_dyb);
%     fprintf('xe4uy_dxb = nvy*C(i+iib,j+jjb,4) = %e\n',xe4uy_dxb);
%     fprintf('xe4uy_dyb = nvx*C(i+iib,j+jjb,2) = %e\n\n',xe4uy_dyb);
%     fprintf('ye4ux_dxb = nvx*C(i+iib,j+jjb,2) = %e\n',ye4ux_dxb);
%     fprintf('ye4ux_dyb = nvy*C(i+iib,j+jjb,4) = %e\n',ye4ux_dyb);
%     fprintf('ye4uy_dxb = nvy*C(i+iib,j+jjb,4) = %e\n',ye4uy_dxb);
%     fprintf('ye4uy_dyb = nvx*C(i+iib,j+jjb,3) = %e\n\n',ye4uy_dyb);
%     
% fprintf('eq5:\n');    
%     fprintf('B1:\n');
%     fprintf('xe5ux_dxxa = nvxnvy*C(i+iia,j+jja,1) = %e\n',xe5ux_dxxa);
%     fprintf('xe5ux_dxya = nvy2*C(i+iia,j+jja,4)-nvx2*C(i+iia,j+jja,1) = %e\n',xe5ux_dxya);
%     fprintf('xe5uy_dxxa = nvy2*C(i+iia,j+jja,4) = %e\n',xe5uy_dxxa);
%     fprintf('xe5uy_dxya = nvxnvy*C(i+iia,j+jja,2)-nvxnvy*C(i+iia,j+jja,4) = %e\n',xe5uy_dxya);
%     fprintf('xe5ux_dyya = - nvxnvy*C(i+iia,j+jja,4) = %e\n',xe5ux_dyya);
%     fprintf('xe5uy_dyya = - nvx2*C(i+iia,j+jja,2) = %e\n\n',xe5uy_dyya);
%     fprintf('ye5ux_dxxa = nvy2*C(i+iia,j+jja,2) = %e\n', ye5ux_dxxa);
%     fprintf('ye5ux_dxya = nvxnvy*C(i+iia,j+jja,4)-nvxnvy*C(i+iia,j+jja,2) = %e\n',ye5ux_dxya);
%     fprintf('ye5uy_dxxa = nvxnvy*C(i+iia,j+jja,4) = %e\n',ye5uy_dxxa);
%     fprintf('ye5uy_dxya = nvy2*C(i+iia,j+jja,3)-nvx2*C(i+iia,j+jja,4) = %e\n',ye5uy_dxya);
%     fprintf('ye5ux_dyya =  - nvx2*C(i+iia,j+jja,4) = %e\n',ye5ux_dyya);
%     fprintf('ye5uy_dyya =  - nvxnvy*C(i+iia,j+jja,3) = %e\n\n',ye5uy_dyya);
% 
%     fprintf('B0:\n');
%     fprintf('xe5ux_dxxb = nvxnvy*C(i+iib,j+jjb,1) = %e\n',xe5ux_dxxb);
%     fprintf('xe5ux_dxyb = nvy2*C(i+iib,j+jjb,4)-nvx2*C(i+iib,j+jjb,1) = %e\n',xe5ux_dxyb);
%     fprintf('xe5uy_dxxb = nvy2*C(i+iib,j+jjb,4) = %e\n',xe5uy_dxxb);
%     fprintf('xe5uy_dxyb = nvxnvy*C(i+iib,j+jjb,2)-nvxnvy*C(i+iib,j+jjb,4) = %e\n',xe5uy_dxyb);
%     fprintf('xe5ux_dyyb = - nvxnvy*C(i+iib,j+jjb,4) = %e\n',xe5ux_dyyb);
%     fprintf('xe5uy_dyyb = - nvx2*C(i+iib,j+jjb,2) = %e\n\n',xe5uy_dyyb);
%     fprintf('ye5ux_dxxb = nvy2*C(i+iib,j+jjb,2) = %e\n',ye5ux_dxxb);
%     fprintf('ye5ux_dxyb = nvxnvy*C(i+iib,j+jjb,4)-nvxnvy*C(i+iib,j+jjb,2) = %e\n',ye5ux_dxyb);
%     fprintf('ye5uy_dxxb = nvxnvy*C(i+iib,j+jjb,4) = %e\n',ye5uy_dxxb);
%     fprintf('ye5uy_dxyb = nvy2*C(i+iib,j+jjb,3)-nvx2*C(i+iib,j+jjb,4) = %e\n',ye5ux_dxyb);
%     fprintf('ye5ux_dyyb =  - nvx2*C(i+iib,j+jjb,4) = %e\n',ye5ux_dyyb);
%     fprintf('ye5uy_dyyb =  - nvxnvy*C(i+iib,j+jjb,3) = %e\n\n',ye5uy_dyyb);
%     end
end