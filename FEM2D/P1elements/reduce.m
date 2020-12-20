function  [nInner, nDirichlet, nNeumann] = reduce(markers)
    
    nInner = 0;
    nDirichlet = nInner;
    nNeumann = nDirichlet;
    
    for j=1:length(markers)
        if markers(j)==0
            nInner = nInner + 1;
        else
            if mod(markers(j),2)==0
                nNeumann = nNeumann + 1;
            else
                nDirichlet = nDirichlet + 1;
            end
        end
    end
    
end