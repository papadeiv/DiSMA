function assignIC(N)

    % import the globally-scoped numerical solution u_h at previous time-step
    global u_h_prev_t;
    % import all the boundary functions on the borders of the domain
    global initial_condition;
    % import geometrical properties of the mesh
    global nodes;
    % loop over the nodes and initialise the IC vector
    for j=1:N
        % assign value to the GLOBAL index
        J = nodes(j,1);
        % check if there's a DOF at the j-th global index
        if J > 0
            u_h_prev_t(J) = initial_condition(nodes(J,3),nodes(J,4));
        end
    end
end