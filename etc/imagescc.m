%% Plotter. Just faster then imagesc

function imagescc(x_vec, z_vec, u, ptitle, pxlabel, pzlabel, axisnormal)   

    if nargin < 7
        axisnormal = false;
    end
    
    imagesc(x_vec,z_vec,u); hold on;
    title(ptitle, 'FontWeight','b', 'FontSize', 16);
    xlabel(pxlabel, 'FontWeight','b', 'FontSize', 14);
    ylabel(pzlabel, 'FontWeight','b', 'FontSize', 14);
    if axisnormal
        set(gca,'YDir','normal');
    end
end