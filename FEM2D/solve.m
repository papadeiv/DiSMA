function u = solve(A, f)
    % define accuracy and maximum number of iterations of the solving algorithms
    toll = 1e-9;
    iter = 1000;
    % initialise the (temporary) solution vector
    Nh = length(f);
    u = zeros(Nh);
    if Nh < 1e7
        % solving the linear system DIRECTLY
        u = A\f;
    else
        % solving the linear system ITERATIVELY
        if issymmetric(A)==1
            % using conjugate gradient descent method (A is symmetric)
            u = pcg(A,f,toll,iter);
        else
            % using minimum squares method (A is not symmetric)
            u = gmres(A,f,[],toll,iter);
        end
    end
    
end