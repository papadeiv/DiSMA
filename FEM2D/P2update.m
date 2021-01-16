function [g_t,f_t] = P2update(Nh, Nd,timestep)
    
    % import all the necessary functions of the variational problem
    global boundary_functions;
    % import geometric entities of the domain
    global triangles;
    global nodes;
    % initialise the vectors
    g_t = zeros(Nd, 1);
    f_t = zeros(Nh, 1);
    % compute the vectors
    for e=1:length(triangles(:,1))
        % extract e-th triangle's geometric data
        area = triangles(e,7);
        for j=1:6
            % check if there exist a trial basis function (i.e. the j-th node of the e-th triangle is a DOF for the linear system)
            if nodes(triangles(e,j),1)>0
                for k=1:6
                    % extract global indices from local indices (j, k) and the triangle index e                                                         
                    j_g = nodes(triangles(e,j),1);
                    k_g = nodes(triangles(e,k),1);
                    % check if there exist a test basis function (i.e. the k-th node of the e-th triangle is a DOF for the linear system)
                    if nodes(triangles(e,k),1)<0
                        k_g = -k_g;
                        % extract the marker of the BC associated to the k_g-th node
                        marker = nodes(triangles(e,k),2);
                        % compute the boundary function value associated to the k_g node
                        g_t(k_g) = boundary_functions{1, marker}(nodes(triangles(e,k),3), nodes(triangles(e,k),4),timestep);
                    end
                end
                % compute the (j_g) entry of the source vector
                f_t(j_g) = f_t(j_g) + 2*area*quad_source(j,timestep);
            end
        end
    end
end