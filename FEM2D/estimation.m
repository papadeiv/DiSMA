function [u,e_inf,e_L2,e_H1] = estimation(u_h, N)

    global triangles;
    global nodes;
    % import exact solution and gradient
    global exact_functions;
    % define analytical solution vector
    u = zeros(N,1);
    % compute analytical solution
    for j=1:N
        u(j) = exact_functions{1}(nodes(j,3),nodes(j,4));
    end
    % compute L_inf norm error
    e_inf = max(abs(u-u_h));
    % import reference triangle nodes and quadrature weights
    [xhat, yhat, omega] = quadrature_weights();
    % define number of quadrature nodes
    Nq = length(xhat);
    % define accumulators for the sum
    acc_L2 = 0;
    acc_H1 = 0;
    % loop over mesh's triangles
    for e=1:length(triangles(:,1))
        % retrieve triangle's area
        area = triangles(e,4);
        % compute the incremental quantities for the e-th triangle
        [dx,dy] = delta(nodes, triangles(e,1:3));
        v1 = nodes(triangles(e,1),3:4);
        % define the error vector of the e-th triangle
        error_L2 = zeros(Nq,1);
        error_H1 = zeros(Nq,1);
        % loop over e-th triangle quadrature nodes
        for q=1:Nq
            % compute new coordinates associated to q-th node through affine map
            mapped_x = v1(1) + dx(3)*xhat(q) - dx(2)*yhat(q);
            mapped_y = v1(2) - dy(3)*xhat(q) + dy(2)*yhat(q);
            % compute basis function on q-th quadrature node
            phi = [1-xhat(q)-yhat(q); xhat(q); yhat(q)];
            % assign (constant) value to basis function's partial derivative
            phi_x = [-1; 1; 0];
            phi_y = [-1; 0; 1];
            % loop over e-th triangle verteces
            sum_L2 = 0;
            sum_H1 = zeros(2,1);
            for j=1:3
                sum_L2 = sum_L2 + u_h(triangles(e,j))*phi(j);
                sum_H1 = sum_H1 + 1/(2*area)*u_h(triangles(e,j))*[dy(j);dx(j)];
            end
            % fill the q-th component of error vector of the e-th triangle
            error_L2(q) = omega(q)*(exact_functions{1}(mapped_x, mapped_y) - sum_L2)^2;
            u_diff = [exact_functions{2}(mapped_x, mapped_y); exact_functions{3}(mapped_x, mapped_y)] - sum_H1;
            error_H1(q) = omega(q)*(transpose(u_diff)*u_diff);
        end
        acc_L2 = acc_L2 + 2*area*sum(error_L2); 
        acc_H1 = acc_H1 + 2*area*sum(error_H1);
    end
    e_L2 = sqrt(acc_L2);
    e_H1 = sqrt(acc_H1);
    
end