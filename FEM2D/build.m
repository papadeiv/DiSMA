function [A, f] = build(Nh, N)

    Nd = N - Nh;
    % import number of DOFs defined on each element (depending on the subspace of the basis functions)
    global Ndof;
    % assemble the linear terms according to the choice of the subspace
    if Ndof==3
        [A, Ad, f, gd] = P1solver(Nh, Nd);
        % adding Neumann's borders nodes contribution to the source vector
        f = P1Neumann(f);
    else
        [A, Ad, f, gd] = P2solver(Nh, Nd);
        % adding Neumann's borders nodes contribution to the source vector
        f = P2Neumann(f);
    end
    % assemble the RHS of the linear system by including the border contribution
    f = f-Ad*gd;
    
end