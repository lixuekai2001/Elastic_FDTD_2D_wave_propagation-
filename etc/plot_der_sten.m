% As input 6 x 9 matrix required

function plot_der_sten(invJaco)
    %clc;
    %close all;
    % invJaco=invJacoc{10,10}
    figure;
    set(gcf,'Name','Derivative stencils');
    set(gcf,'NumberTitle','off');
    for n=1:6
        subplot(2,3,n);
        ctr=0;
         for ii=1.d0:-1.d0:-1.d0
            for jj=1.d0:-1.d0:-1.d0
                ctr=ctr+1;
                %fprintf('i=%d j=%d\n', ii, jj);
                scatter(ii,jj, 100*abs(invJaco(n, ctr))/abs(max(invJaco(n,:))),'filled','b'); hold on;  
                switch n
                    case 1
                        title('u');
                    case 2
                        title('u,x');
                    case 3
                        title('u,y');
                    case 4
                        title('u,xx');
                    case 5
                        title('u,yy');
                    case 6
                        title('u,xy');
                end
                drawnow;
            end
         end
    end
end