%% Cerjan 1985 absorbing boundary

function  r = Cerjan_BC_weights(rrate, lmargin, rmargin,nnx,nnz)

      % Cerjan boundary conditions (2D)
    %   double precision :: r(nnx,nnz)
    %   real*8, intent(in) :: rrate
    %   integer, dimension(3), intent(in) :: lmargin, rmargin
    %   integer ix, iy, iz
    %   integer i, j, k
    %   double precision :: rr
  r = ones(nnx,nnz);
  
  for iz = 1:nnz
     for ix = 1:nnx
        
        i = 0;
        j = 0;
        k = 0;
           
        if (ix < lmargin(1) + 1)
            i = lmargin(1) + 1 - ix;
        end
        %   if (iy < lmargin(2) + 1) j = lmargin(2) + 1 - iy
        if (iz < lmargin(2) + 1)
            k = lmargin(2) + 1 - iz;
        end
   
        if (nnx - rmargin(1) < ix)
            i = ix - nnx + rmargin(1);
        end
        %if (nny - rmargin(2) < iy) j = iy - nny + rmargin(2)
        if (nnz - rmargin(2) < iz)
            k = iz - nnz + rmargin(2);
        end
           
        if (i == 0 && j == 0 && k == 0)
            continue
        end
        
        rr = rrate * rrate * double( i * i + j * j + k * k );
        r(ix,iz) = exp(-rr);
        
        %if(r(ix,iz).ne.1.d0) then
        %   print *, ix,iz,r(ix,iz)
        %endif
        
        %print *, ix,iy,r
        %ux2(:,:) = ux2(:,:) * r
        %ux1(:,:) = ux1(:,:) * r
        %ux(:,:) = ux(:,:) * r

        %uz2(:,:) = uz2(:,:) * r
        %uz1(:,:) = uz1(:,:) * r
        %uz(:,:) = uz(:,:) * r
        
     end
  end
end