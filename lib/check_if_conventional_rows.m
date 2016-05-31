%% Check presence of conventional operator rows

function [Aux, Auy] = check_if_conventional_rows(Aux, Auy, i, j, C, DELTAX, DELTAY)

    conv_stecncil_a2 = [1.0 -2.0 1.0 0.0]/(DELTAX^2);
    conv_stecncil_ab = [1.0 -1.0 -1.0 1.0]/(4*DELTAX*DELTAY);
    
    Auxc(1,:) = Aux(1,:)/C(i,j,1);
    Auxc(2,:) = Aux(2,:)/C(i,j,4);
    Auxc(3,:) = Aux(3,:)/C(i,j,2);
    Auxc(4,:) = Aux(4,:)/C(i,j,4);

    Auyc(1,:) = Auy(1,:)/C(i,j,4);
    Auyc(2,:) = Auy(2,:)/C(i,j,3);
    Auyc(3,:) = Auy(3,:)/C(i,j,4);
    Auyc(4,:) = Auy(4,:)/C(i,j,2);
    
    Auxc = round(Auxc*DELTAX^2)/DELTAX^2;
    
    conv_rows_a2 = ismember(Auxc, conv_stecncil_a2, 'rows');
    conv_rows_ab = ismember(Auxc, conv_stecncil_ab, 'rows');

    if nnz(conv_rows_a2)>0 || nnz(conv_rows_ab)>0
        [Auxz, ~] = construct_Zahradnik_operators(i,j,C, DELTAX, DELTAY);
        ind_a2 = find(conv_rows_a2);
        ind_ab = find(conv_rows_ab);
        if ~isempty(ind_a2)
            Aux(ind_a2,:) = Auxz(ind_a2,:);
        end
        if ~isempty(ind_ab)
            Aux(ind_ab,:) = Auxz(ind_ab,:);
        end
%         disp('checked1');
    end
    
    conv_stecncil_a2 = [1.0 -2.0 1.0 0.0]/(DELTAY^2);
    Auyc = round(Auyc*DELTAY^2)/DELTAY^2;
    conv_rows_a2 = ismember(Auyc, conv_stecncil_a2, 'rows');
    conv_rows_ab = ismember(Auyc, conv_stecncil_ab, 'rows');

    if nnz(conv_rows_a2)>0 || nnz(conv_rows_ab)>0
        [~, Auyz] = construct_Zahradnik_operators(i,j,C, DELTAX, DELTAY);
        ind_a2 = find(conv_rows_a2);
        ind_ab = find(conv_rows_ab);
        Auy(ind_a2,:) = Auyz(ind_a2,:);
        Auy(ind_ab,:) = Auyz(ind_ab,:);
%         disp('checked2');
    end
    
end