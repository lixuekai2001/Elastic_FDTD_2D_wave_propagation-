%% Construct_interface_operators [Ovcharenko et al., 2015]

function [Aux, Auy] = construct_interface_operators(i,j, gr_x, gr_y, xdscr, ydscr, C, rho)
        % Construct eta0 and eta1 arrays for each marked point
        
        dx = gr_x(i,j) - gr_x(i-1,j);
        dy = gr_y(i,j) - gr_y(i,j-1);
        
        ZERO = 0.d0;
        max_size_xdscr = max(size(xdscr));
        min_size_xdscr = min(size(xdscr));
        pt0x=gr_x(i,j);
        pt0y=gr_y(i,j);
        ctr=0;
        for ik=1:-1:-1 %over columns of points from right to left
            for jk=1:-1:-1  %over rows of points from top to bottom
                    pt1x=gr_x(i+ik,j);
                    pt1y=gr_y(i,j+jk);
                    ctr=ctr+1;
                    x_trial=linspace(pt0x,pt1x,20);
                    y_trial=linspace(pt0y,pt1y,20);
                    [xi,yi]=curveintersect(x_trial,y_trial,xdscr, ydscr);
                    
                    if ~isempty([xi,yi]) % check if there is an intersection
                        if size(xi,1)*size(xi,2)>1  %get rid of multiple intersections
                            xi=xi(1);
                            yi=yi(1);
                        end
                        delta_x=abs(pt1x-pt0x);
                        delta_y=abs(pt1y-pt0y);
                        if delta_x<eps || delta_y<eps
                            if delta_x<eps %vertical line
                                arr_eta0y(i,j,ctr)=abs(yi-pt0y)/delta_y;
                                arr_eta1y(i,j,ctr)=1.d0-arr_eta0y(i,j,ctr);
                                arr_eta0x(i,j,ctr)=ZERO;
                                arr_eta1x(i,j,ctr)=ZERO;                                    
                            end
                            if delta_y<eps %hoizontal line
                                arr_eta0x(i,j,ctr)=abs(xi-pt0x)/delta_x;
                                arr_eta1x(i,j,ctr)=1.d0-arr_eta0x(i,j,ctr);
                                arr_eta0y(i,j,ctr)=ZERO;
                                arr_eta1y(i,j,ctr)=ZERO;
                            end
                            if delta_x<eps && delta_y<eps %if single point
                                arr_eta0y(i,j,ctr)=ZERO;
                                arr_eta1y(i,j,ctr)=ZERO;
                                arr_eta0x(i,j,ctr)=ZERO;
                                arr_eta1x(i,j,ctr)=ZERO;
                            end
                        else %if evrything ok with deltas
                                arr_eta0y(i,j,ctr)=abs((yi-pt0y)/delta_y);
                                arr_eta1y(i,j,ctr)=1.d0-arr_eta0y(i,j,ctr);
                                arr_eta0x(i,j,ctr)=abs((xi-pt0x)/delta_x);
                                arr_eta1x(i,j,ctr)=1.d0-arr_eta0x(i,j,ctr);
                        end
%                         if length(xdscr)==49 && i==99 && j==2 && ctr ==3
%                             disp('Here');
%                         end
                        %Define normal in point
                        tmp=abs(xdscr-xi);
                        [~,idx]=min(tmp);
                        if idx==1 
                            idx=2;
                        elseif idx >= max_size_xdscr
                            idx = max_size_xdscr - 1;
                        elseif idx <= min_size_xdscr
                            idx = min_size_xdscr - 1;
                        end
                        p1x=xdscr(idx-1); p2x=xi; p3x=xdscr(idx+1);
                        p1y=ydscr(idx-1); p2y=yi; p3y=ydscr(idx+1);                 
%                         if ((p1x==p2x) && (p1y==p2y)) || ((p3x==p2x) && (p3y==p2y))
%                             p1x=xdscr(idx-2);
%                             p1y=ydscr(idx-2);
%                         end
                        s12 = sqrt((p2x-p1x)^2+(p2y-p1y)^2);
                        s23 = sqrt((p3x-p2x)^2+(p3y-p2y)^2);
                        dxds = (s23^2*(p2x-p1x)+s12^2*(p3x-p2x))/(s12*s23*(s12+s23));
                        dyds = (s23^2*(p2y-p1y)+s12^2*(p3y-p2y))/(s12*s23*(s12+s23));
                        if isnan(dxds) || isnan(dyds)
                           s13 = sqrt((p3x-p1x)^2+(p3y-p1y)^2);
                           dxds = (p3x - p1x)/s13;
                           dyds = (p3y - p1y)/s13;
                        end
%                         tvx=dxds; % tangent components
%                         tvy=dyds;
                        nvx=-dyds;  % normal components
                        nvy=dxds;
                    else  % if there is no intersection
                            arr_eta1x(i,j,ctr)=1.d0*abs(ik);
                            arr_eta1y(i,j,ctr)=1.d0*abs(jk);
                            arr_eta0x(i,j,ctr)=abs(1.d0*ik)*abs((1.d0-arr_eta1x(i,j,ctr)));
                            arr_eta0y(i,j,ctr)=abs(1.d0*jk)*abs((1.d0-arr_eta1y(i,j,ctr)));
                            nvx=0.d0;
                            nvy=0.d0;
                    end   % checking if there are intersections

                    %Apply Mizutani operators
                    eta0x = arr_eta0x(i,j,ctr);
                    eta1x = arr_eta1x(i,j,ctr);
                    eta0y = arr_eta0y(i,j,ctr);
                    eta1y = arr_eta1y(i,j,ctr);                    

                    [A0,B0,A1,B1]=A0B0A1B1triso2(i,j,ik,jk,0,0,nvx,nvy,C,rho,dx,dy, ik*eta0x,-ik*eta1x, jk*eta0y, -jk*eta1y); % + 
                    CJI=svdinv(B0*A0)*(B1*A1);
                    coeffAux(ctr,:)=CJI(1,1:6);
                    coeffAuy(ctr,:)=CJI(7,7:12);
            end      %end of jk loop
        end  %%end of ik loop

        pre_dx2=svdinv([coeffAux(2,[1,2,4]); coeffAux(5,[1,2,4]); coeffAux(8,[1,2,4])]);
        pre_dy2=svdinv([coeffAux(4,[1,3,5]); coeffAux(5,[1,3,5]); coeffAux(6,[1,3,5])]);
        pre_dxdy=svdinv([coeffAux(1,:); coeffAux(3,:); coeffAux(7,:); coeffAux(9,:)]);
        coeffux_dx2 = pre_dx2(3,:);
        coeffux_dy2 = pre_dy2(3,:);
        coeffux_dxdy= pre_dxdy(6,:);
        coeffux_dydx = coeffux_dxdy;

        pre_dx2=svdinv([coeffAuy(2,[1,2,4]); coeffAuy(5,[1,2,4]); coeffAuy(8,[1,2,4])]);
        pre_dy2=svdinv([coeffAuy(4,[1,3,5]); coeffAuy(5,[1,3,5]); coeffAuy(6,[1,3,5])]);
        pre_dxdy=svdinv([coeffAuy(1,:); coeffAuy(3,:); coeffAuy(7,:); coeffAuy(9,:)]);
        coeffuy_dx2 = pre_dx2(3,:);
        coeffuy_dy2 = pre_dy2(3,:);
        coeffuy_dxdy= pre_dxdy(6,:);
        coeffuy_dydx = coeffuy_dxdy;

        Aux=[[coeffux_dx2 0]; [coeffux_dy2 0]; coeffux_dxdy; coeffux_dydx];
        Auy=[[coeffuy_dx2 0]; [coeffuy_dy2 0]; coeffuy_dxdy; coeffuy_dydx];
end
