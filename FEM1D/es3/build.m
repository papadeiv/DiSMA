function [A,b] = build(partition,nodeslist,mesh_size,nDOFs)
    
    % define coefficient functions
    mu = @(x) 1+x;
    mu_x = @(x) 1;
    sigma = @(x) exp(x);
    beta = @(x) 3*x;
    % define source function
    f = @(x) 1;
    % initialise linear system
    D = zeros(nDOFs,nDOFs);
    C = zeros(nDOFs,nDOFs);
    R = zeros(nDOFs,nDOFs);
    b = zeros(nDOFs,1);
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
            end
        end
    end
    A = D + C + R;
    b(1) = b(1) - 2*(mu(0)/mesh_size+sigma(0)*mesh_size);
    b(nDOFs) = b(nDOFs) - 5*(mu(2)/mesh_size+sigma(2)*mesh_size);
end