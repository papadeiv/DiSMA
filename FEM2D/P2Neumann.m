function f = P2Neumann(f, timestep)
    
    % import all the boundary functions on the borders of the domain
    global boundary_functions;
    % import geometric entities of the domain
    global nodes;
    global borders;
    % computing the contribution to the source vector for each Neumann's border
    for b=1:size(borders,1)
        
        % extract boundary marker of the border
        marker = borders(b,6);
        
        % extract index and coordiantes of the vertices
        % BEGINNING EDGE VERTEX
        Vb = borders(b,2);
        x_b = nodes(Vb,3);
        y_b = nodes(Vb,4);
        % MEDIUM EDGE VERTEX
        Vm = borders(b,3);
        x_m = nodes(Vm,3);
        y_m = nodes(Vm,4);
        % ENDING EDGE VERTEX
        Ve = borders(b,4);
        x_e = nodes(Ve,3);
        y_e = nodes(Ve,4);
        
        % define the affine parametric map for the reference border
        gamma_x = @(that) x_b + (x_e - x_b)*that;
        gamma_y = @(that) y_b + (y_e - y_b)*that;
        
        % define the basis functions restriction on the vertices
        % BEGINNING EDGE VERTEX
        phi_b = @(that) 2*that^2-3*that+1;
        % MEDIUM EDGE VERTEX
        phi_m = @(that) -4*that^2+4*that;
        % ENDING EDGE VERTEX
        phi_e = @(that) 2*that^2-that;
        
        % evaluate teh Neumann function on the mapped coordinate of each vertex
        % BEGINNING EDGE VERTEX
        g_b = boundary_functions{2, marker}(gamma_x(0), gamma_y(0), timestep);
        % MEDIUM EDGE VERTEX
        g_m = boundary_functions{2, marker}(gamma_x(0.5), gamma_y(0.5), timestep);
        % ENDING EDGE VERTEX
        g_e = boundary_functions{2, marker}(gamma_x(1), gamma_y(1), timestep);
        
        % compute the contribution of each node to the source vector with numerical quadrature
        % BEGINNING EDGE VERTEX
        if nodes(Vb,1)>0
            cav_simpson = (1/6)*g_b*phi_b(0) + (2/3)*g_m*phi_b(0.5) + (1/6)*g_e*phi_b(1);
            f(nodes(Vb,1)) = f(nodes(Vb,1)) + cav_simpson*borders(b,5);
        end
        % MEDIUM EDGE VERTEX
        if nodes(Vm,1)>0
            cav_simpson = (1/6)*g_b*phi_m(0.5) + (2/3)*g_m*phi_m(0.5) + (1/6)*g_e*phi_m(1);
            f(nodes(Vm,1)) = f(nodes(Vm,1)) + cav_simpson*borders(b,5);
        end
        % ENDING EDGE VERTEX
        if nodes(Ve,1)>0
            cav_simpson = (1/6)*g_b*phi_e(0) + (2/3)*g_m*phi_e(0.5) + (1/6)*g_e*phi_e(1);
            f(nodes(Ve,1)) = f(nodes(Ve,1)) + cav_simpson*borders(b,5);
        end
        
    end

end