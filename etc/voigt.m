% Voigt mapping
%Maps 4d tesor to 2d
function res=voigt(C)
C=squeeze(C);
order=max(size(C));
order2=order^2;
res=zeros(2*order);
for i=1:order
    for j=1:order
        for k=1:order
            for l=1:order
                p=i*dkr(i,j)+(1-dkr(i,j))*(order2-i-j);
                q=k*dkr(k,l)+(1-dkr(k,l))*(order2-k-l);
                res(p,q)=C(i,j,k,l);
            end
        end
    end
end

end