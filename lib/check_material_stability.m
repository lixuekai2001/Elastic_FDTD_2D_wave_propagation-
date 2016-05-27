%% Check anisotropic material stability

function check_material_stability(C)

    [NX, NY] = size(C(:,:,1));

    cond1 = 0;
    cond2 = 0;
    cond3 = 0;
    cond4 = 0;

    for i = 1:NX
        for j = 1:NY
        c11 = C(i,j,1);
        c13 = C(i,j,2);
        c33 = C(i,j,3);
        c44 = C(i,j,4);

        %From Roland Martin code
        %from Becache et al., INRIA report, equation 7 page 5 http://hal.inria.fr/docs/00/07/22/83/PDF/RR-4304.pdf
          if(c11*c33 - c13*c13 <= 0.d0)
              cond1 = 1;
    %           disp('There is a problem in definition of orthotropic material');
              %break
          end

        %check intrinsic mathematical stability of PML model for an anisotropic material
        %from E. B\'ecache, S. Fauqueux and P. Joly, Stability of Perfectly Matched Layers, group
        %velocities and anisotropic waves, Journal of Computational Physics, 188(2), p. 399-433 (2003)
          aniso_stability_criterion = ((c13+c44)^2 - c11*(c33-c44)) * ((c13+c44)^2 + c44*(c33-c44));
    %       fprintf('PML anisotropy stability criterion from Becache et al. 2003 = %e\n', aniso_stability_criterion);
          if(aniso_stability_criterion > 0.d0) % && (USE_PML_XMIN  ||  USE_PML_XMAX  ||  USE_PML_YMIN  ||  USE_PML_YMAX)
              cond2 = 1;
    %           fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 1');
             %break
          end

          aniso2 = (c13 + 2*c44)^2 - c11*c33;
    %       fprintf('PML aniso2 stability criterion from Becache et al. 2003 = %e\n',aniso2);
          if(aniso2 > 0.d0) %&& (USE_PML_XMIN  ||  USE_PML_XMAX  ||  USE_PML_YMIN  ||  USE_PML_YMAX)
             cond3 = 1;
    %          fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 2');
             %break
          end

          aniso3 = (c13 + c44)^2 - c11*c33 - c44^2;
    %       fprintf('PML aniso3 stability criterion from Becache et al. 2003 = %e\n',aniso3);
          if(aniso3 > 0.d0) %&& (USE_PML_XMIN  ||  USE_PML_XMAX  ||  USE_PML_YMIN  ||  USE_PML_YMAX)
              cond4 = 1;
    %           fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 3');
             %break
          end
        end
    end
    if cond1
        disp('There is a problem in definition of orthotropic material');
    end

    if cond2
        fprintf('PML anisotropy stability criterion from Becache et al. 2003 = %e\n', aniso_stability_criterion);
        fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 1');
    end

    if cond3
        fprintf('PML aniso2 stability criterion from Becache et al. 2003 = %e\n',aniso2);
        fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 2');
    end

    if cond4
        fprintf('PML aniso3 stability criterion from Becache et al. 2003 = %e\n',aniso3);
        fprintf('WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 3');
    end
    
    if cond1 || cond2 || cond3 || cond4
        fprintf('...FAIL\n');
        fprintf('\n');
    else
        fprintf('...OK\n');
        fprintf('\n');
    end
end