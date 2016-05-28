%% Cerjan 1985 absorbing boundary
% Fuji N. 2016 (Fortran)
% Ovcharenko O. 2016 (MATLAB)

function  [ux, uy] = Cerjan_absorbing_BC(ux, uy, rate, lmargin, rmargin)

  [~, NX, NZ] = size(ux);

  weights = ones(NX,NZ);
  
  for iz = 1:NZ
     for ix = 1:NX
        
        i = 0;
        j = 0;
        k = 0;
           
        if (ix < lmargin(1) + 1)
            i = lmargin(1) + 1 - ix;
        end

        if (iz < lmargin(2) + 1)
            k = lmargin(2) + 1 - iz;
        end
   
        if (NX - rmargin(1) < ix)
            i = ix - NX + rmargin(1);
        end

        if (NZ - rmargin(2) < iz)
            k = iz - NZ + rmargin(2);
        end
           
        if (i == 0 && j == 0 && k == 0)
            continue
        end
        
        rr = rate * rate * double(i*i + j*j + k*k );
        weights(ix,iz) = exp(-rr);
        
     end
  end
  
  for ix = 1:NX
      for iz = 1:NZ
          ux(3,ix,iz) = ux(3,ix,iz) * weights(ix,iz);
          ux(2,ix,iz) = ux(2,ix,iz) * weights(ix,iz);
          ux(1,ix,iz) = ux(1,ix,iz) * weights(ix,iz);

          uy(3,ix,iz) = uy(3,ix,iz) * weights(ix,iz);
          uy(2,ix,iz) = uy(2,ix,iz) * weights(ix,iz);
          uy(1,ix,iz) = uy(1,ix,iz) * weights(ix,iz);
      end
  end
  
end