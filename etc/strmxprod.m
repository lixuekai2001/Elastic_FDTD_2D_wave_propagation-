function C=strmxprod(A,B)
SHOW_RES=false;

[nRow1,nCol1]=size(A);
[nRow2,nCol2]=size(B);
if(nCol1~=nRow2)  || nRow1==0 || nRow2==0
  C=[]; 
  disp('ERROR! Check sizes of input matrices');
  return
end

C=cell(nRow1,nCol2);
    for i=1:nRow1;
        for j=1:nCol2;
            ctr=0;
            for p=1:nCol1;
                if isempty(B{p,j}) || isempty(A{i,p})
                    continue
                end
                if strcmp('0',num2str(A{i,p})) || strcmp('0',num2str(B{p,j}))
                    continue
                end
                ctr=ctr+1;
                if ctr==1
                    C{i,j}=strcat(C{i,j},strcat(num2str(A{i,p}),'*',num2str(B{p,j})));
                else
                    str=[' + ' num2str(A{i,p}) '*' num2str(B{p,j})];
                    C{i,j}=strcat(C{i,j},str);
                end
            end            
        end
    end
    
if SHOW_RES
    [s1,s2]=size(C);
    for i=1:s1
        fprintf('\n');
        for j=1:s2
            fprintf('%s  ',C{i,j});
        end
    end
    fprintf('\n');
end    
    
    
end