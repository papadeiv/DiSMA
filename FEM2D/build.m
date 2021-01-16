function [] = build(Nh, N)
    
    % import flag valued dummy variable
    global t;
    % define the source vector
    global f;
    % compute number of Dirichlet's nodes
    Nd = N - Nh;
    % import number of DOFs defined on each element (depending on the subspace of the basis functions)
    global Ndof;
    % assemble the linear terms according to the choice of the subspace
    if Ndof==3
        if t==0
            P1builder(Nh, Nd);
            % adding Neumann's borders nodes contribution to the source vector
            f = P1Neumann(f,t);
        else
            unsteadyP1(Nh,Nd);
        end
    else
        if t==0
            P2builder(Nh, Nd);
            % adding Neumann's borders nodes contribution to the source vector
            f = P2Neumann(f,t);
        else
            unsteadyP2(Nh, Nd);
        end
    end
    
end