
function [markers, xt_dis, yt_dis, nvecx, nvecy, xmn, ymn] = func_p_find_closest_grid_nodes(nx,ny,nc,gr_x,gr_y ,x_topo, y_topo)
% PLT_ON=false;
% fprintf('Started looking for markers, x,y-t_dis, nvec, xmn,ymn with\n');
topo_szx=size(x_topo,2)-1;
tgrx=round(topo_szx/nx);
% DELTAX=gr_x(2,2)-gr_x(1,2);
% DELTAY=gr_y(2,2)-gr_y(2,1);
% xt_dis=[x_topo(1)]; %x of discretized topography
% yt_dis=[y_topo(1)]; %y of discretized topography
% nvecx=zeros(nx+1,ny+1);
% nvecy=zeros(nx+1,ny+1);
% xmn=zeros(nx+1,ny+1);
% ymn=zeros(nx+1,ny+1);
markers=zeros(nx+1,ny+1);
for i=1:nx
    x_trial=(1+(i-1)*tgrx):(i*tgrx); %vector of curve x values
    y_trial=y_topo(x_trial);    
    for j=1:ny
%         xs=[];
%         ys=[];
%         str=[num2str(i) ' ' num2str(j)];
        if ~isempty(find(y_trial<gr_y(i,j+1), 1)) && ~isempty(find(y_trial>gr_y(i,j), 1))
            markers(i,j)=1;
            markers(i,j+1)=1;
            markers(i+1,j)=1;
            markers(i+1,j+1)=1;
% %Verical lines | |          
%             yvert=linspace(gr_y(1,j)-0.01*DELTAY,gr_y(1,j+1)+0.01*DELTAY,size(x_trial,2));
%             xvert=gr_x(i,j)*ones(size(yvert));
%             
%             if PLT_ON
%                 plot(xvert,yvert); hold on;
%                 text(gr_x(i,j)+DELTAX/4,gr_y(i,j)+DELTAY/2,str);
%             end
%             [x0,y0]=curveintersect(xvert,yvert,x_topo, y_topo);
%             if ~isempty([x0,y0]) && isempty(find(xs==x0, 1)) && isempty(find(ys==y0, 1))
%                     xt_dis=[xt_dis, x0];
%                     yt_dis=[yt_dis, y0];
% %                     text(x0,y0,[num2str(x0) ' ' num2str(y0)]);
%             if PLT_ON
%                     scatter(x0,y0,'b');
%             end
%                     xs=[xs, x0];
%                     ys=[ys, y0];
%             end
%             yvert=linspace(gr_y(1,j)-0.01*DELTAY,gr_y(1,j+1)+0.01*DELTAY,size(x_trial,2));
%             xvert=gr_x(i+1,j)*ones(size(yvert));
%             if PLT_ON
%                 plot(xvert,yvert); hold on;
%             end
%             [x0,y0]=curveintersect(xvert,yvert,x_topo, y_topo);
%             if ~isempty([x0,y0]) && isempty(find(xs==x0, 1)) && isempty(find(ys==y0, 1))
%                     xt_dis=[xt_dis, x0];
%                     yt_dis=[yt_dis, y0];
% %                     text(x0,y0,[num2str(x0) ' ' num2str(y0)]);
%                     if PLT_ON
%                         scatter(x0,y0,'b');
%                     end
%                     xs=[xs, x0];
%                     ys=[ys, y0];
%             end
% % Horizontal lines =
%             xhor=linspace(gr_x(i,1)-0.01*DELTAX,gr_x(i+1,1)+0.01*DELTAX,size(x_trial,2));
%             yhor=gr_y(i,j)*ones(size(xhor));
%             if PLT_ON
%                 plot(xhor,yhor); hold on;
%             end
%             [x0,y0]=curveintersect(xhor,yhor,x_topo, y_topo);
%             if ~isempty([x0,y0]) && isempty(find(xs==x0, 1)) && isempty(find(ys==y0, 1))
%                     xt_dis=[xt_dis, x0];
%                     yt_dis=[yt_dis, y0]; 
%                     if PLT_ON
%                         scatter(x0,y0,'b');   
%                     end
%                     xs=[xs, x0];
%                     ys=[ys, y0];
%             end
%             xhor=linspace(gr_x(i,1)-0.01*DELTAX,gr_x(i+1,1)+0.01*DELTAX,size(x_trial,2));
%             yhor=gr_y(i,j+1)*ones(size(xhor));
%             if PLT_ON
%                 plot(xhor,yhor); hold on;
%             end
%             [x0,y0]=curveintersect(xhor,yhor,x_topo, y_topo);
%             if ~isempty([x0,y0]) && isempty(find(xs==x0, 1)) && isempty(find(ys==y0, 1))
%                     xt_dis=[xt_dis, x0];
%                     yt_dis=[yt_dis, y0]; 
%                     if PLT_ON
%                         scatter(x0,y0,'b');
%                     end
%                     xs=[xs, x0];
%                     ys=[ys, y0];
%             end
%             if length(xs)==2 && length(ys)==2
%                     [xs_sorted,xs_index]=sort(xs);
%                     ys_sorted=ys(xs_index);
%                     xs=xs_sorted;
%                     ys=ys_sorted;
%                     x1n=xs(1);
%                     x2n=xs(2);
%                     y1n=ys(1);
%                     y2n=ys(2);
%                     cxmn=(x1n+x2n)/2;
%                     cymn=(y1n+y2n)/2;
%                     xmn(i,j)=cxmn;
%                     ymn(i,j)=cymn;
%                     x0n=x2n-x1n;
%                     y0n=y2n-y1n;
%                     cnorm=[y0n, -x0n]/sqrt(x0n^2+y0n^2);
%                     nvecx(i,j)=cnorm*[1 0]';
%                     nvecy(i,j)=cnorm*[0 1]' ;
%                     if PLT_ON
%                         line([cxmn cxmn-DELTAX*nvecx(i,j)],[cymn cymn-DELTAY*nvecy(i,j)],'Color','m'); hold on;
%                         drawnow;
%                     end
%             else
%                disp('ERROR. Grid is too sparse! You have >2 normals in one cell');
%             end
         end
    end
end
fprintf('N of involved near-boundary grid points = %d  %.2f%%\n', nnz(markers), nnz(markers)*100/(nx*ny));
% %Sort xt_dis  and ut_dis vector to get rid of loops on curve
%  [xt_dis_sorted, xt_sortindex] = sort(xt_dis);
%   yt_dis_sorted = yt_dis(xt_sortindex);  
%  xt_dis=xt_dis_sorted;
%  yt_dis=yt_dis_sorted;
%  
%  if PLT_ON
%     plot(xt_dis,yt_dis); hold on;
%     axis([min(gr_x(:,1)) max(gr_x(:,1)) -inf inf])
%     title('Normals to the interface given at each cell');
%     xlabel('m');
%     ylabel('m');
%  end

xt_dis=0; yt_dis=0; nvecx=0; nvecy=0; xmn=0; ymn = 0;
% fprintf('Calculating of normals finished\n');
end