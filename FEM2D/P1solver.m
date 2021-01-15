function P1solver(Nh, Nd)
    
    % import time-step
    global t;
    % import all the necessary functions of the variational problem
    global coefficient_functions;
    global source_function;
    global boundary_functions;
    % import geometric entities of the domain
    global triangles;
    global nodes;
    % initialise the linear system
    global M;
    M = zeros(Nh, Nh);
    global A;
    A = zeros(Nh, Nh);
    global f;
    f = zeros(Nh, 1);
    global Ad;
    Ad = zeros(Nh, Nd);
    global gd;
    gd = zeros(Nd,1);
    % assemble the linear system
    for e=1:length(triangles(:,1))
        % extract e-th triangle's geometric data
        area = triangles(e,4);
        x_G = triangles(e,5);
        y_G = triangles(e,6);
        % compute the incremental quantities for the e-th triangle
        [dx,dy] = delta(e);
        for j=1:3
            % check if there exist a trial basis function (i.e. the j-th node of the e-th triangle is a DOF for the linear system)
            if nodes(triangles(e,j),1)>0
                for k=1:3
                    % extract global indices from local indices (j, k) and the triangle index e                                                         
                    j_g = nodes(triangles(e,j),1);
                    k_g = nodes(triangles(e,k),1);
                    % check if there exist a test basis function (i.e. the k-th node of the e-th triangle is a DOF for the linear system)
                    if nodes(triangles(e,k),1)>0
                        % compute the entry of A for the j_g-th trial basis function and the k_g-th test basis functions
                        A(j_g,k_g) = A(j_g,k_g) + ...
                            coefficient_functions{1}(x_G,y_G,t)*(dy(k)*dy(j)+dx(k)*dx(j))/(4*area) + ...
                            0 + ...
                            0;
                    else
                        % convert the Dirichlet pivot into natural integer
                        k_g = -k_g;
                        % extract the marker of the BC associated to the k_g-th node
                        marker = nodes(triangles(e,k),2);
                        if marker == 0
                            error('One DOF has been wrongfuly stored as a boundary node');
                        end
                        % compute the entry of Ad for the j_g-th trial basis function and the k_g-th test basis functions
                        Ad(j_g,k_g) = Ad(j_g,k_g) + ...
                            coefficient_functions{1}(x_G,y_G,t)*(dy(k)*dy(j)+dx(k)*dx(j))/(4*area) + ... % diffusive term
                            0 + ... % convective term
                            0; % radiative term
                        % compute the boundary function value associated to the k_g node
                        gd(k_g) = boundary_functions{1, marker}(nodes(triangles(e,k),3), nodes(triangles(e,k),4),t);
                    end
                end
                % compute the (j_g) entry of the source vector
                f(j_g) = f(j_g) + source_function(x_G,y_G,t)*(area/3);
            end
        end
    end

end