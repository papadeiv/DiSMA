function [A, f] = build(Nh, N, triangles, borders, nodes)

    Nd = N - Nh;
    % import all the necessary functions of the variational problem
    global coefficient_functions;
    global boundary_functions;
    global source_function;
    % initialise linear system
    A = zeros(Nh, Nh);
    f = zeros(Nh, 1);
    Ad = zeros(Nh, Nd);
    gd = zeros(Nd,1);
    for e=1:length(triangles(:,1))
        % extract e-th triangle's geometric data
        area = triangles(e,4);
        x_G = triangles(e,5);
        y_G = triangles(e,6);
        % compute the incremental quantities for the e-th triangle
        [dx,dy] = delta(nodes,triangles(e,1:3));
        for j=1:3
            % check if there exist a trial basis function (i.e. the j-th node of the e-th triangle is a DOF for the linear system)
            if nodes(triangles(e,j),1)>0
                for k=1:3
                    % extract global indices from local indices (j, k) and the triangle index e                                                         
                    [j_g, k_g] = map(nodes(triangles(e,j),1), nodes(triangles(e,k),1));
                    % check if there exist a test basis function (i.e. the k-th node of the e-th triangle is a DOF for the linear system)
                    if nodes(triangles(e,k),1)>0
                        % compute the entry of A for the j_g-th trial basis function and the k_g-th test basis functions
                        A(j_g,k_g) = A(j_g,k_g) + ...
                            coefficient_functions{1}(x_G,y_G)*(dy(k)*dy(j)+dx(k)*dx(j))/(4*area) + ...
                            0 + ...
                            0;
                    else
                        % convert the Dirichlet pivot into natural integer
                        k_g = -k_g;
                        % extract the marker of the BC associated to the k_g-th node
                        marker = nodes(triangles(e,k),5);
                        if marker == 0
                            error('One DOF has been wrongfuly stored as a boundary node');
                        end
                        % compute the entry of Ad for the j_g-th trial basis function and the k_g-th test basis functions
                        Ad(j_g,k_g) = Ad(j_g,k_g) + ...
                            coefficient_functions{1}(x_G,y_G)*(dy(k)*dy(j)+dx(k)*dx(j))/(4*area) + ...
                            0 + ...
                            0;
                        % compute the boundary function value associated to the k_g node
                        gd(k_g) = boundary_functions{1, marker}(nodes(triangles(e,k),3), nodes(triangles(e,k),4));
                    end
                end
                % compute the (j_g) entry of the source vector
                f(j_g) = f(j_g) + source_function(x_G,y_G)*(area/3);
            end
        end
    end
    % adding Neumann's borders nodes contribution to the source vector
    for b=1:size(borders,1)
        % extract coordinates of vertices of the b-th border
        Ve = borders(b,2);
        x_e = nodes(Ve,3);
        y_e = nodes(Ve,4);
        Vb = borders(b,3);
        x_b = nodes(Vb,3);
        y_b = nodes(Vb,4);
        % evaluate the Neumann function on the vertices
        g_e = boundary_functions{2, borders(b,5)}(x_e, y_e);
        g_b = boundary_functions{2, borders(b,5)}(x_b, y_b);
        % add the contirbute
        if nodes(Vb,1)>0
            f(nodes(Vb,1)) = f(nodes(Vb,1)) + borders(b,4)*(g_b/3 + g_e/6);
        end
        if nodes(Ve,1)>0
            f(nodes(Ve,1)) = f(nodes(Ve,1)) + borders(b,4)*(g_b/6 + g_e/3);
        end
    end
    % assemble the RHS of the linear system by including the border contribution
    f = f-Ad*gd;
    
end