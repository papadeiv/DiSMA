function [h_max, A_max, Nt] = discretize(idx, subspace)

    ref_area = 0.2/2^idx;
    mesh(ref_area, subspace);

    global geom;
    global Ndof;
    
    global triangles;
    global nodes;
    global borders;
    
    Nt = geom.nelements.nTriangles;
    triangles = [];
    triangles(:,1:Ndof) = geom.elements.triangles; % global indeces of triangles vertices
    for e=1:Nt
        triangles(e,Ndof+1) = geom.support.TInfo(e).Area; % surface areas of triangles 
        triangles(e,Ndof+2) = geom.support.TInfo(e).CG(1); % x_CG of triangles
        triangles(e,Ndof+3) = geom.support.TInfo(e).CG(2); % y_CG of triangles
    end

    if Ndof==6
        borders = zeros(size(geom.pivot.Ne,1),7);
        counter = 0;
        for b = 1:geom.nelements.nBorders
            % check if it is a border of Neumann's boundary
            if ismember(b,geom.pivot.Ne(:,1))
                counter = counter + 1;
                % fill-in the border's index
                borders(counter, 1) = b;
                % fill-in border's vertices
                borders(counter, 2:4) = [geom.elements.borders(b,1), geom.elements.borders(b,5), geom.elements.borders(b,2)];
                % compute and fill-in border's length (L2-norm)
                Vb = [geom.elements.coordinates(borders(counter,2),1),geom.elements.coordinates(borders(counter,2),2)];
                Ve = [geom.elements.coordinates(borders(counter,4),1),geom.elements.coordinates(borders(counter,4),2)];
                borders(counter, 5) = norm(Vb-Ve,2);
                % fill-in border's BC marker
                idx = find(b==geom.pivot.Ne(:,1));
                borders(counter, 6) = find(geom.input.BC.Boundary.Values==geom.pivot.Ne(idx,2));
                % fill in the parent triangle for the boundary border
                if geom.elements.borders(borders(counter, 1),3) > 0
                    borders(counter, 7) = geom.elements.borders(borders(counter, 1), 3);
                else
                    borders(counter, 7) = geom.elements.borders(borders(counter, 1), 4);
                end
            end
            % save the borders' length
            L(b)=sqrt((geom.elements.coordinates(geom.elements.borders(b,1),1)-geom.elements.coordinates(geom.elements.borders(b,2),1))^2+(geom.elements.coordinates(geom.elements.borders(b,1),2)-geom.elements.coordinates(geom.elements.borders(b,2),2))^2);     
        end
        % derive mesh's longest border
        h_max = max(L);
        
    else
        borders = zeros(size(geom.pivot.Ne,1),5);
        counter = 0;
        for b = 1:geom.nelements.nBorders
            % check if it is a border of Neumann's boundary
            if ismember(b,geom.pivot.Ne(:,1))
                counter = counter + 1;
                % fill-in the border's index
                borders(counter, 1) = b;
                % fill-in border's vertices
                borders(counter, 2:3) = [geom.elements.borders(b,1), geom.elements.borders(b,2)];
                % compute and fill-in border's length (L2-norm)
                Vb = [geom.elements.coordinates(borders(counter,2),1),geom.elements.coordinates(borders(counter,2),2)];
                Ve = [geom.elements.coordinates(borders(counter,3),1),geom.elements.coordinates(borders(counter,3),2)];
                borders(counter, 4) = norm(Vb-Ve,2);
                % fill-in border's BC marker
                idx = find(b==geom.pivot.Ne(:,1));
                borders(counter, 5) = find(geom.input.BC.Boundary.Values==geom.pivot.Ne(idx,2));
            end
            % save the borders' length
            L(b)=sqrt((geom.elements.coordinates(geom.elements.borders(b,1),1)-geom.elements.coordinates(geom.elements.borders(b,2),1))^2+(geom.elements.coordinates(geom.elements.borders(b,1),2)-geom.elements.coordinates(geom.elements.borders(b,2),2))^2);     
        end
        % derive mesh's longest border
        h_max = max(L);
        
    end
    % derive mesh's largest (triangles') area
    A_max = max(triangles(:,4));
    
    nodes = [];
    nodes(:,1) = geom.pivot.pivot; % global indeces of DOFs and BCs
    nodes(:,2) = boundary(geom.pivot.nodelist, geom.input.BC.Boundary.Values); % BC border marker
    nodes(:, 3:4) = geom.elements.coordinates; % (x,y)s of each node in DOFs and BCs

end

