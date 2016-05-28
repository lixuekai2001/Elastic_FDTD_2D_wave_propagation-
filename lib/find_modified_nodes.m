
function markers = find_modified_nodes(nx,nz,gr_x,gr_z ,interfaces)

    DELTAX=gr_x(3,2) - gr_x(2,2);  %[m]
    DELTAZ=gr_z(2,3) - gr_z(2,2);  %[m]

    num_of_interfaces = size(interfaces,1);
    markers=zeros(nx+1,nz+1);

    
%     for inter = 1:num_of_interfaces
%         topo_szx=size(x_topo,2)-1;
%         tgrx=round(topo_szx/nx);
% 
%         for i=1:nx
%             x_trial=(1+(i-1)*tgrx):(i*tgrx);
%             y_trial=y_topo(x_trial);    
%             for j=1:nz
%                 if ~isempty(find(y_trial<gr_y(i,j+1), 1)) && ~isempty(find(y_trial>gr_y(i,j), 1))
%                     markers(i,j)=1;
%                     markers(i,j+1)=1;
%                     markers(i+1,j)=1;
%                     markers(i+1,j+1)=1;
%                  end
%             end
%         end
%         fprintf('N of involved near-boundary grid points = %d  %.2f%%\n', nnz(markers), nnz(markers)*100/(nx*nz));
%     end

%     for i=1:nx+1
%         for j=1:nz+1
%             scatter(gr_x(i,j),gr_z(i,j)); hold on;
%         end
%     end
%     drawnow;
    
    for k = 1:num_of_interfaces
          int_x = interfaces{k,1};
          int_z = interfaces{k,2};
%         topo_szx=size(x_topo,2)-1;
%         tgrx=round(topo_szx/nx);
          range_x = int_x/DELTAX + 1;
          range_z = int_z/DELTAZ + 1;
          range_z_min = floor(range_z);
          range_z_max = ceil(range_z);
          
          len_range_x = length(range_x);
          len_range_z_min = length(range_z_min);
          for i = 1:len_range_x
%               for j = 1:len_range_z_min
                  markers(range_x(i), range_z_min(i)) = k;
                  markers(range_x(i), range_z_max(i)) = k;
                  xcoord_min = gr_x(range_x(i),range_z_min(i));
                  zcoord_min = gr_z(range_x(i),range_z_min(i));
                  scatter(xcoord_min,zcoord_min,'filled','r'); hold on;
                  
                  xcoord_max = gr_x(range_x(i),range_z_max(i));
                  zcoord_max = gr_z(range_x(i),range_z_max(i));
                  scatter(xcoord_max,zcoord_max, 'filled','r'); hold on;
                  drawnow;
%               end
          end
%         for i=range_x
%             x_trial=(1+(i-1)*tgrx):(i*tgrx);
%             y_trial=y_topo(x_trial);    
%             for j=1:nz
%                  if ~isempty(find(y_trial<gr_y(i,j+1), 1)) && ~isempty(find(y_trial>gr_y(i,j), 1))
%                     markers(i,j)=k;
%                     markers(i,j+1)=k;
%                     markers(i+1,j)=k;
%                     markers(i+1,j+1)=k;
%                  end
%             end
%         end
    end
    fprintf('N of involved near-boundary grid points = %d  %.2f%%\n', nnz(markers), nnz(markers)*100/((nx+1)*(nz+1)));
%     imagesc(flipud(markers'));
end