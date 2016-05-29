
function markers = find_modified_nodes(NX,NZ,gr_x,gr_z ,interfaces)
    NX_max = NX + 1;
    NZ_max = NZ + 1;
    
    DELTAX=gr_x(3,2) - gr_x(2,2);  %[m]
    DELTAZ=gr_z(2,3) - gr_z(2,2);  %[m]

    num_of_interfaces = size(interfaces,1);
    markers=zeros(NX_max,NZ_max);

    hchannel = DELTAZ/3.d0;
    wchannel = DELTAX/3.d0;
    
    fprintf('Define nodes where modified operator has to be imposed\n');
    fprintf('Number of interfaces = %d\n', num_of_interfaces);
    fprintf('Channel height = %.4f\n', hchannel);
    fprintf('Channel width = %.4f\n', wchannel);
    for k = 1:num_of_interfaces
          % Construct channel around interface
          channel_z_up = interfaces{k,2} + hchannel;
          channel_z_down = interfaces{k,2} - hchannel;
          channel_x_left = interfaces{k,1} - wchannel;
          channel_x_right = interfaces{k,1} + wchannel;
          
          % Left side
          range_x_left = channel_x_left/DELTAX + 1;
          range_x_left_min = floor(range_x_left);
          range_x_left_max = ceil(range_x_left);
          range_x_left_min(range_x_left_min == 0) = 1;
          range_x_left_max(range_x_left_max > NX_max) = NX_max;
          
          % Right side
          range_x_right = channel_x_right/DELTAX + 1;
          range_x_right_min = floor(range_x_right);
          range_x_right_max = ceil(range_x_right);
          range_x_right_min(range_x_right_min == 0) = 1;
          range_x_right_max(range_x_right_max > NX_max) = NX_max;
          
          % Up side
          range_z_up = channel_z_up/DELTAZ + 1;
          range_z_up_min = floor(range_z_up);
          range_z_up_max = ceil(range_z_up);
          range_z_up_min(range_z_up_min == 0) = 1;
          range_z_up_max(range_z_up_max > NZ_max) = NZ_max;
          
          % Down side
          range_z_down = channel_z_down/DELTAZ + 1;
          range_z_down_min = floor(range_z_down);
          range_z_down_max = ceil(range_z_down);
          range_z_down_min(range_z_down_min == 0) = 1;
          range_z_down_max(range_z_down_max > NZ_max) = NZ_max;
          
          len_range_x = length(range_x_right);
          for i = 1:len_range_x
                  markers(range_x_left_min(i), range_z_up_min(i)) = k;
                  markers(range_x_left_min(i), range_z_up_max(i)) = k;
                  markers(range_x_left_min(i), range_z_down_min(i)) = k;
                  markers(range_x_left_min(i), range_z_down_max(i)) = k;
                  markers(range_x_left_max(i), range_z_up_min(i)) = k;
                  markers(range_x_left_max(i), range_z_up_max(i)) = k;
                  markers(range_x_left_max(i), range_z_down_min(i)) = k;
                  markers(range_x_left_max(i), range_z_down_max(i)) = k;
                  
                  markers(range_x_right_min(i), range_z_up_min(i)) = k;
                  markers(range_x_right_min(i), range_z_up_max(i)) = k;
                  markers(range_x_right_min(i), range_z_down_min(i)) = k;
                  markers(range_x_right_min(i), range_z_down_max(i)) = k;
                  markers(range_x_right_max(i), range_z_up_min(i)) = k;
                  markers(range_x_right_max(i), range_z_up_max(i)) = k;
                  markers(range_x_right_max(i), range_z_down_min(i)) = k;
                  markers(range_x_right_max(i), range_z_down_max(i)) = k;
                  
%                   xcoord_min = gr_x(range_x(i),range_z_up_min(i));
%                   zcoord_min = gr_z(range_x(i),range_z_up_min(i));
%                   scatter(xcoord_min,zcoord_min,'filled','r'); hold on;
%                   
%                   xcoord_max = gr_x(range_x(i),range_z_up_max(i));
%                   zcoord_max = gr_z(range_x(i),range_z_up_max(i));
%                   scatter(xcoord_max,zcoord_max, 'filled','r'); hold on;
%                   
%                   xcoord_min = gr_x(range_x(i),range_z_down_min(i));
%                   zcoord_min = gr_z(range_x(i),range_z_down_min(i));
%                   scatter(xcoord_min,zcoord_min,'filled','b'); hold on;
%                   
%                   xcoord_max = gr_x(range_x(i),range_z_down_max(i));
%                   zcoord_max = gr_z(range_x(i),range_z_down_max(i));
%                   scatter(xcoord_max,zcoord_max, 'filled','b'); hold on;
%                   drawnow;
%               end
          end
    end
    fprintf('%d involved points or %.2f%% ', nnz(markers), nnz(markers)*100/((NX_max)*(NZ_max)));
    fprintf('...OK\n\n');
%     imagesc(flipud(markers'));
end