function [] = main(Ns, subspace)
    
    % initialize the error vectors
    E_inf = zeros(Ns+1,1);
    E_L2 = E_inf;
    E_H1 = E_L2;
    % initialize the mesh size convergence vectors
    h = zeros(Ns+1,1); % longest edge
    Ar = zeros(Ns+1,1); % largest area
    Nh = zeros(Ns+1,1); % number of DOFs
    % create full-screen window for solution's subplots
    figure('units','normalized','outerposition',[0 0 1 1])
    % iterations over the grid refinement
    for m=0:Ns
        % derive a discrete grid object on top of the continuos domain
        [h(m+1),Ar(m+1),Nt] = discretize(m, subspace);
        % extract the number of DOFs and BCs out of all the nodes
        [Ni,Nd,Nn] = reduce();
        % number of DOFs
        Nh(m+1) = Ni + Nn;
        % number of total nodes(==geom.nelements.nVerterxes)
        N = Nh(m+1) + Nd;
        % print informations at screen in runtime
        disp(fprintf("Iteration %d\nMesh info:\n     Number of (triangular) elements:              %d\n     Total number of nodes:                        %d\n     Number of DOFs(of which Neumann's):           %d(%d)\n",...
            m+1, Nt, N, Nh(m+1), Nn));
        % assembling the linear system associated to the discrete grid
        [A, f] = build(Nh(m+1), N);
        % solving the linear system with ad-hoc algorithm
        u_tmp = solve(A, f);
        % reconstruct the whole solution by including the (known) values at the boundaries
        u_h = expand(u_tmp, N);
        % estimate error (u_h-u) in different norms
        [u,E_inf(m+1),E_L2(m+1),E_H1(m+1)] = estimation(u_h, N);
        % plot the reconstructed (discrete) solution
        multiplot(Ns, m, u_h, u);
    end
    % estimate and plot the errors' convergence
    errors(h, Ar, Nh, E_inf, E_L2, E_H1);
    
end



