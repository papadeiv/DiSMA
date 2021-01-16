function assemble(Nh, Nd, scheme, dt)
    
    % import current time-step
    global t;
    % import number of DOFs defined on each element (depending on the subspace of the basis functions)
    global Ndof;
    % import all the terms of the built linear system
    global A;
    global Ad;
    global f;
    global gd;
    % define the LHS (unknown's coefficient matrix) and the RHS (known's coefficients vector) 
    global LHS;
    global RHS;
    
    if dt==0 % NO TIME_DEPENDENCE
        
       % compute the left-hand side of the final linear system
       LHS = A;
       % compute the right-hand side of the final linear system
       RHS = f-Ad*gd;
       
    else % TIME_DEPENDENCE
        
        % import numerical solution at previous time-step
        global u_h_prev_t;
        % import time-dependent linear system terms
        global M;
        global Md;
        % check against the time-discretisation scheme 
        if scheme == "1st-order"
            % compute the source and Dirichlet's nodes vector
            if Ndof==3
                [gd_prev_t, f_prev_t] = P1update(Nh, Nd, t-dt);
                f_prev_t = P1Neumann(f_prev_t, t-dt);
                [gd_t, f_t] = P1update(Nh, Nd, t);
                f_t = P1Neumann(f_t, t);
            else
                [gd_prev_t, f_prev_t] = P2update(Nh, Nd, t-dt);
                f_prev_t = P2Neumann(f_prev_t, t-dt);
                [gd_t, f_t] = P2update(Nh, Nd, t);
                f_t = P2Neumann(f_t, t);
            end
            % compute left and right hand side of the final linear system
            LHS = M + dt*A;
            RHS = M*u_h_prev_t - Md*(gd_t - gd_prev_t) + dt*(f_t - Ad*gd_t);
        else
            % compute the source and Dirichlet's nodes vector
            if Ndof==3
                [gd_prev_t, f_prev_t] = P1update(Nh, Nd, t-dt);
                f_prev_t = P1Neumann(f_prev_t, t-dt);
                [gd_t, f_t] = P1update(Nh, Nd, t);
                f_t = P1Neumann(f_t, t);
            else
                [gd_prev_t, f_prev_t] = P2update(Nh, Nd, t-dt);
                f_prev_t = P2Neumann(f_prev_t, t-dt);
                [gd_t, f_t] = P2update(Nh, Nd, t);
                f_t = P2Neumann(f_t, t);
            end
            % compute left and right hand side of the final linear system
            LHS = M + (dt/2)*A;
            % b = (1/2)*deltaT*(G_k+G_k_1) + (1/2)*deltaT*(F_k+F_k_1) - (Md*(uhd_k_1-uhd_k)')' - ((1/2)*deltaT*Ad*(uhd_k_1+uhd_k)')' + ((M-(1/2)*deltaT*A)*uh0_k')';
            RHS = M*u_h_prev_t - Md*(gd_t - gd_prev_t) + (dt/2)*((f_t + f_prev_t) - A*u_h_prev_t - Ad*(gd_t + gd_prev_t));
        end 
    end
    
end