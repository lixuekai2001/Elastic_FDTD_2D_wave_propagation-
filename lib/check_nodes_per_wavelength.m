%% Check number of nodes per wavelength
% cp has to be minimum possible
% f0 has to be maximum possible
% Ovcharenko O. 2016

function check_nodes_per_wavelength(cp, f0, DELTAX, DELTAY)

    wavelength = cp/f0;
    nodes_per_wavelength_x = wavelength/DELTAX;
    nodes_per_wavelength_y = wavelength/DELTAY;
    fprintf(' Shortest wavelength = %.2f m\n Nodes per wavelength:\n \t %.2f OX\n \t %.2f OY\n', wavelength, nodes_per_wavelength_x, nodes_per_wavelength_y);

    if nodes_per_wavelength_x < 10 || nodes_per_wavelength_y < 10
        disp('Too few nodes per wavelength. Decrease f0 or increase NX and NY');
    end
    fprintf('\n');
end