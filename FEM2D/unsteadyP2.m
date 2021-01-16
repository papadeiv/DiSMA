function unsteadyP2(Nh, Nd)
    
    % import geometric entities of the domain
    global triangles;
    global nodes;
    % initialise linear system
    global M;
    M = zeros(Nh, Nh);
    global Md;
    Md = zeros(Nh, Nd);
    global A;
    A = zeros(Nh, Nh);
    global Ad;
    Ad = zeros(Nh, Nd);
    % assemble the linear system
    for e=1:length(triangles(:,1))
        % extract e-th triangle's area
        area = triangles(e,7);
        % compute the incremental quantities for the e-th triangle
        [dx,dy] = delta(e);
        % assemble the jacobian of the mapping
        B = [dx(3) -dx(2); -dy(3) dy(2)];
        diffB = inv(B)*inv(transpose(B));
        for j=1:6
            % check if there exist a trial basis function (i.e. the j-th node of the e-th triangle is a DOF for the linear system)
            if nodes(triangles(e,j),1)>0
                for k=1:6
                    % extract global indices from local indices (j, k) and the triangle index e
                    j_g = nodes(triangles(e,j),1);
                    k_g = nodes(triangles(e,k),1);
                    % check if there exist a test basis function (i.e. the k-th node of the e-th triangle is a DOF for the linear system)
                    if nodes(triangles(e,k),1)>0
                        % compute the entry of A for the j_g-th trial basis function and the k_g-th test basis functions
                        M(j_g,k_g) = M(j_g,k_g) + 0;
                        % compute the entry of A for the j_g-th trial basis function and the k_g-th test basis functions
                        A(j_g,k_g) = A(j_g,k_g) + ...
                            2*area*quad_diffusive(j,k,diffB) + ...
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
                        
                        % compute the entry of A for the j_g-th trial basis function and the k_g-th test basis functions
                        Md(j_g,k_g) = Md(j_g,k_g) + 0;
                        % compute the entry of Ad for the j_g-th trial basis function and the k_g-th test basis functions
                        Ad(j_g,k_g) = Ad(j_g,k_g) + ...
                            2*area*quad_diffusive(j,k,diffB) + ...
                            0 + ...
                            0;
                    end
                end
            end
        end
    end

end