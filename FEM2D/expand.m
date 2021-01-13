function u_h = expand(u, N)

    global t;
    global nodes;
    % check if there are error in the vectors' dimensions
    if length(nodes()) ~= N
        error('Total number of nodes (N) is different from the number of elements in the expanded list (nodes(:,1))');
    end
    % initialise the solution vector
    u_h = zeros(N,1);
    % import boundary functions
    global boundary_functions;
    % loop over the mesh nodes' indices
    for j=1:N
        % assign value to the GLOBAL index
        J = nodes(j,1);
        % check if there's a DOF at the j-th global index
        if J > 0
            u_h(j) = u(J);
        % otherwise assign the associated boundary value
        else
            % convert the Dirichlet pivot into natural integer
            J = -J;
            marker = nodes(j,2);
            if marker == 0
               error('One DOF has been wrongfully stored as a boundary node');
            end
            u_h(j) = boundary_functions{1, marker}(nodes(j,3), nodes(j,4),t);
        end
    end
    
end