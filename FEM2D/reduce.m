function  [nInner, nDirichlet, nNeumann] = reduce()
    
    global geom;
    global borders;
    nInner = 0;
    nDirichlet = nInner;
    nNeumann = nDirichlet;
    
    for j=1:length(geom.pivot.nodelist)
        if geom.pivot.nodelist(j)==0
            if size(geom.pivot.Ne,1) > 0 && (ismember(j, borders(:,2)) || ismember(j, borders(:,3)))
                nNeumann = nNeumann + 1;
            else
                nInner = nInner + 1;
            end
        else
            nDirichlet = nDirichlet + 1;
        end
    end
    
end