%% Check CFL stability condition.

function check_CFL(cp, DELTAT, DELTAX, DELTAY)
% R. Courant et K. O. Friedrichs et H. Lewy (1928)

  if nargin<4
      DELTAY = 0.d0;
  end
  
  Courant_number = cp * DELTAT * sqrt(1.d0/DELTAX^2.d0 + 1.d0/DELTAY^2.d0);
  
  fprintf('Courant number = %.4f',Courant_number);

  if Courant_number > 1.d0 
      fprintf('...FAIL. Time step is too large, simulation will be unstable.');
      check = false;
  else
    fprintf('...OK\n');
    check = true;
  end
    fprintf('\n');
end