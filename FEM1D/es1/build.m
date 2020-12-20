function [A,b,u] = build(partition,nodeslist,mesh_size,nDOFs)
    
    % define coefficient functions
    mu = @(x) 1;
    mu_x = @(x) 0;
    sigma = @(x) 1;
    beta = @(x) 1;
    % define extact solution (if known) and its prime derivative
    U = @(x) sin(pi*x);
    U_x = @(x) pi*cos(pi*x);
    U_xx = @(x) -pi^2*sin(pi*x);
    % define source function
    f = @(x) -(mu_x(x)*U_x(x)+mu(x)*U_xx(x)) + beta(x)*U_x(x) + sigma(x)*U(x);
    F = @(x) 2*sin(x) + cos(x);
    % initialise linear system and (exact) solution vector
    D = zeros(nDOFs,nDOFs);
    C = zeros(nDOFs,nDOFs);
    R = zeros(nDOFs,nDOFs);
    b = zeros(nDOFs,1);
    u = zeros(nDOFs,1);
    % loop on each finite element (sub-interval) of the grid (partition)
    nElements = length(partition(:,1));
    for e=1:nElements
        for j=1:2
            j_glb = partition(e,j);
            x = nodeslist(1,j_glb);
            mrk_j = nodeslist(2,j_glb);
            if mrk_j >= 0
                for k=1:2
                    k_glb = partition(e,k);
                    mrk_k = nodeslist(2,k_glb);
                    if mrk_k >= 0
                        assemble;
                    end
                end
                b(mrk_j) = f(x)*mesh_size;
                u(mrk_j) = U(x);
            end
        end
    end
    A = D + C + R;
    
end