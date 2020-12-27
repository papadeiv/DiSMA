function f = P2Neumann(f)

    % import all the boundary functions on the borders of the domain
    global boundary_functions;
    % import geometric entities of the domain
    global nodes;
    global borders;
    % computing the contribution to the source vector for each Neumann's border
    for b=1:size(borders,1)
        % extract coordinates of vertices and the medium node of the b-th border
        Vb = borders(b,2);
        x_b = nodes(Vb,3);
        y_b = nodes(Vb,4);
        Vm = borders(b,3);
        Ve = borders(b,4);
        x_e = nodes(Ve,3);
        y_e = nodes(Ve,4);
        % evaluate the Neumann function on the vertices
        g_e = boundary_functions{2, borders(b,6)}(x_e, y_e);
        g_b = boundary_functions{2, borders(b,6)}(x_b, y_b);
        % add the contribution of the beginning vertex
        if nodes(Vb,1)>0
            f(nodes(Vb,1)) = f(nodes(Vb,1)) + 0.5*borders(b,5)*(g_b/3 + g_e/6);
        end
        % add the contribution of the medium node
        f(nodes(Vm,1)) = f(nodes(Vm,1)) + 0.5*borders(b,5)*(1/3)*(g_b + g_e);
        % add the contribution of the ending vertex
        if nodes(Ve,1)>0
            f(nodes(Ve,1)) = f(nodes(Ve,1)) + 0.5*borders(b,5)*(g_b/6 + g_e/3);
        end
    end

end