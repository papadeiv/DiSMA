function f = P1Neumann(f, timestep)

    % import all the boundary functions on the borders of the domain
    global boundary_functions;
    % import geometric entities of the domain
    global nodes;
    global borders;
    % computing the contribution to the source vector for each Neumann's border
    for b=1:size(borders,1)
        % extract coordinates of vertices of the b-th border
        Ve = borders(b,2);
        x_e = nodes(Ve,3);
        y_e = nodes(Ve,4);
        Vb = borders(b,3);
        x_b = nodes(Vb,3);
        y_b = nodes(Vb,4);
        % evaluate the Neumann function on the vertices
        g_e = boundary_functions{2, borders(b,5)}(x_e, y_e, timestep);
        g_b = boundary_functions{2, borders(b,5)}(x_b, y_b, timestep);
        % add the contirbute
        if nodes(Vb,1)>0
            f(nodes(Vb,1)) = f(nodes(Vb,1)) + borders(b,4)*(g_b/3 + g_e/6);
        end
        if nodes(Ve,1)>0
            f(nodes(Ve,1)) = f(nodes(Ve,1)) + borders(b,4)*(g_b/6 + g_e/3);
        end
    end

end