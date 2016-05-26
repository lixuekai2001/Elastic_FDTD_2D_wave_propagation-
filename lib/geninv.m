%Generelized matrix inversion
% Least-squares solution
% (AT A)-1 * A

function B=geninv(A)
    B=(A'*A)\A';
end