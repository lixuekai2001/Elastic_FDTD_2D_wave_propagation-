%% construct Zahradnik operators
function [Aux, Auy] = construct_Zahradnik_operators(i,j,C, DELTAX, DELTAY)

            dxdy4 = 4.d0*DELTAX*DELTAY;
            one_over_2dx2 = 1.d0/(2.d0*DELTAX^2.d0);
            one_over_2dy2 = 1.d0/(2.d0*DELTAY^2.d0);
            one_over_2dxdy4 = 1.d0/(2.d0*dxdy4);

            %ux_dx2
            ml = C(i-1,j,1)+C(i,j,1);                  %left point 
            mr = C(i,j,1)+C(i+1,j,1);                  %right point
            mc = C(i-1,j,1) +2.d0*C(i,j,1)+C(i+1,j,1); %middle point
            coeffux_dx2 = one_over_2dx2*[ml -mc mr];
            %uy_dx2
            ml = C(i-1,j,4)+C(i,j,4);                  %left point 
            mr = C(i,j,4)+C(i+1,j,4);                  %right point
            mc = C(i-1,j,4) +2.d0*C(i,j,4)+C(i+1,j,4); %middle point
            coeffuy_dx2 = one_over_2dx2*[ml -mc mr];
            %ux_dy2
            ml = C(i,j-1,4)+C(i,j,4);                  %left point 
            mr = C(i,j,4)+C(i,j+1,4);                  %right point
            mc = C(i,j-1,4) +2.d0*C(i,j,4)+C(i,j+1,4); %middle point           
            coeffux_dy2 = one_over_2dy2*[ml -mc mr];
            %uy_dy2
            ml = C(i,j-1,3)+C(i,j,3);                  %left point 
            mr = C(i,j,3)+C(i,j+1,3);                  %right point
            mc = C(i,j-1,3) +2.d0*C(i,j,3)+C(i,j+1,3); %middle point  
            coeffuy_dy2 = one_over_2dy2*[ml -mc mr];
            
%             ux_dxdy
            mp1p1 = C(i+1,j+1,2)+C(i,j,2);
            mp1m1 = C(i+1,j-1,2)+C(i,j,2);
            mm1p1 = C(i-1,j+1,2)+C(i,j,2);
            mm1m1 = C(i-1,j-1,2)+C(i,j,2);
            coeffux_dxdy= one_over_2dxdy4*[mp1p1 -mp1m1 -mm1p1 mm1m1];
            
%             uy_dydx
%             mp1p1 = C(i+1,j+1,2)+C(i,j,2);
%             mp1m1 = C(i+1,j-1,2)+C(i,j,2);
%             mm1p1 = C(i-1,j+1,2)+C(i,j,2);
%             mm1m1 = C(i-1,j-1,2)+C(i,j,2);
            coeffuy_dydx= one_over_2dxdy4*[mp1p1 -mp1m1 -mm1p1 mm1m1];
            
%             uy_dxdy
            mp1p1 = C(i+1,j+1,4)+C(i,j,4);
            mp1m1 = C(i+1,j-1,4)+C(i,j,4);
            mm1p1 = C(i-1,j+1,4)+C(i,j,4);
            mm1m1 = C(i-1,j-1,4)+C(i,j,4);
            coeffuy_dxdy= one_over_2dxdy4*[mp1p1 -mp1m1 -mm1p1 mm1m1];
            
%             ux_dydx
%             mp1p1 = C(i+1,j+1,4)+C(i,j,4);
%             mp1m1 = C(i+1,j-1,4)+C(i,j,4);
%             mm1p1 = C(i-1,j+1,4)+C(i,j,4);
%             mm1m1 = C(i-1,j-1,4)+C(i,j,4);
            coeffux_dydx= one_over_2dxdy4*[mp1p1 -mp1m1 -mm1p1 mm1m1];
            
%             coeffux_dx2 = C(i,j,1)*tmp_dx2;
%             coeffux_dy2 = C(i,j,4)*tmp_dy2;
%             coeffux_dxdy= C(i,j,2)*tmp_dxdy;
      
%             coeffuy_dx2 = C(i,j,4)*tmp_dx2;
%             coeffuy_dy2 = C(i,j,3)*tmp_dy2;
%             coeffuy_dxdy= C(i,j,4)*tmp_dxdy;
            Aux = [[coeffux_dx2 0]; [coeffux_dy2 0]; coeffux_dxdy; coeffux_dydx];
            Auy = [[coeffuy_dx2 0]; [coeffuy_dy2 0]; coeffuy_dxdy; coeffuy_dydx];
end