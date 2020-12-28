function f = P2Neumann(f)
    %import the basis functions for the reference triagngle
    [phi, grad_phi] = basis();
    % import all the boundary functions on the borders of the domain
    global boundary_functions;
    % import geometric entities of the domain
    global triangles;
    global nodes;
    global borders;
    % computing the contribution to the source vector for each Neumann's border
    for b=1:size(borders,1)
        
        % BEGINNING EDGE VERTEX
        % extract index and coordinates of the vertex
        Vb = borders(b,2);
        x_b = nodes(Vb,3);
        y_b = nodes(Vb,4);
        % exctract marker of the vertex
        if nodes(Vb,1)>0
            marker_b = 1;
        else
            marker_b = 0;
        end
        % evaluate the basis functions associated to each node of the border
        j_b = find(triangles(borders(b,7),1:6)==Vb);
        phi_b = phi{j_b}(x_b, y_b);
        % evaluate the Neumann function
        g_b = boundary_functions{2, borders(b,6)}(x_b, y_b);
        
        % MEDIUM EDGE VERTEX
        % evaluate the basis functions associated to each node of the border
        Vm = borders(b,3);
        x_m = nodes(Vm,3);
        y_m = nodes(Vm,4);
        % exctract marker of the vertex
        if nodes(Vm,1)>0
            marker_m = 1;
        else
            marker_m = 0;
        end
        % evaluate the basis functions associated to each node of the border
        j_m = find(triangles(borders(b,7),1:6)==Vm);
        phi_m = phi{j_m}(0.5, 0.5);
        % evaluate the Neumann function
        g_m = boundary_functions{2, borders(b,6)}(0.5, 0.5);
        
        % ENDING EDGE VERTEX
        % extract index and coordinates of the vertex
        Ve = borders(b,4);
        x_e = nodes(Ve,3);
        y_e = nodes(Ve,4);
        % exctract marker of the vertex
        if nodes(Ve,1)>0
            marker_e = 1;
        else
            marker_e = 0;
        end
        % evaluate the basis functions associated to each node of the border
        j_e = find(triangles(borders(b,7),1:6)==Ve);
        phi_e = phi{j_e}(x_e, y_e);
        % evaluate the Neumann function
        g_e = boundary_functions{2, borders(b,6)}(x_e, y_e);
        
        f(abs(nodes(Vb,1))) = f(abs(nodes(Vb,1))) + marker_b*g_m*phi_m*borders(b,5);
        f(abs(nodes(Vm,1))) = f(abs(nodes(Vm,1))) + marker_m*g_m*phi_m*borders(b,5);
        f(abs(nodes(Ve,1))) = f(abs(nodes(Ve,1))) + marker_e*g_m*phi_m*borders(b,5);
    end

end