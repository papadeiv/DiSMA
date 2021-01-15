function assemble(scheme, dt)
    
    % import current time-step
    global t;
    % define previous time-step
    prev_t = t - dt;
    % import all the terms of the built linear system
    global M;
    global A;
    global f;
    global Ad;
    global gd;
    % define the LHS (unknown's coefficient matrix) and the RHS (known's coefficients vector) 
    global LHS;
    global RHS;
    if dt > 0 % there is time-dependence
        if scheme == "1st-order"
            a = prev_t;
        else
            a = prev_t;
        end
        
        f = a*f;
        LHS = M + A;
        RHS = f;
        
    else % there is no time-dependence
        % update the source vector by including Dirichlet's nodes contribution
        f = f-Ad*gd;
        % compute the left-hand side of the final linear system
        LHS = M + A;
        % compute the right-hand side of the final linear system
        RHS = f;
    end
    
end