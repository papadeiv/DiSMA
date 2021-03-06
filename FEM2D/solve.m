function u = solve()
    
    % import the assembled RHS and LHS terms 
    global LHS;
    global RHS;
    % define accuracy and maximum number of iterations of the solving algorithms
    toll = 1e-9;
    iter = 1000;
    % initialise the (temporary) solution vector
    Nh = length(RHS);
    if Nh < 1e6
        % solving the linear system DIRECTLY
        u = LHS\RHS;
    else
        disp('\n**** Using iterative methods for the solution of the linear system ****\n');
        % solving the linear system ITERATIVELY
        if issymmetric(LHS)==1
            % using conjugate gradient descent method (A is symmetric)
            u = pcg(LHS,RHS,toll,iter);
        else
            % using minimum squares method (A is not symmetric)
            u = gmres(LHS,RHS,[],toll,iter);
        end
    end
    
end