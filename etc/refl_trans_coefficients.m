%% Reflection/Transmission coefficients
function [Refl_coef, Trans_coef] = refl_trans_coefficients(cp_a, cp_b, rho_a, rho_b)
    %Reflection and transition coefficients
    Refl_coef=(rho_b*cp_b-rho_a*cp_a)/(rho_b*cp_b+rho_a*cp_a);
    Trans_coef=2.d0*rho_b*cp_b/(rho_b*cp_b+rho_a*cp_a);
    if Refl_coef<ZERO
      tmps=', inverse polarity';
    else
      tmps='';
    end
    fprintf('Below --> Above:\n');
    fprintf('  R= %.2f - reflection%s\n  T= %.2f - transmition\n', Refl_coef, tmps ,Trans_coef);
    fprintf('\n');
end