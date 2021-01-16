function [] = main(time, multi, subspace, scheme, Ns)
    
    % create full-screen window for solution's subplots
    figure('units','normalized','outerposition',[0 0 1 1])
    % create dummy variable for time-step and time increment
    global t;
    t = 0;
    dt = 0;
    % check if the model is unstationary or steady
    if time==0
        % check if grid-convergence is enabled
        if multi=='Y'
            s=0;
        else
            s=Ns;
        end
        % define the error and mesh size vectors length
        l = (Ns - s) + 1;
        % initialize the error vectors
        E_inf = zeros(l,1);
        E_L2 = E_inf;
        E_H1 = E_L2;
        % initialize the mesh size convergence vectors
        h = zeros(l,1); % longest edge
        Ar = zeros(l,1); % largest area
        Nh = zeros(l,1); % number of DOFs    
        % iterations over the grid refinement
        for m=0:l-1
            % derive a discrete grid object on top of the continuos domain
            if multi=='Y'
                [h(m+1),Ar(m+1),Nt] = discretize(m, subspace);
            else
                [h(m+1),Ar(m+1),Nt] = discretize(Ns, subspace);
            end
            % extract the number of DOFs and BCs out of all the nodes
            [Ni,Nd,Nn] = reduce();
            % number of DOFs
            Nh(m+1) = Ni + Nn;
            % number of total nodes
            N = Nh(m+1) + Nd;
            % print informations at screen in runtime
            disp(sprintf("\n Grid-convergence iteration: %d",m+1));
            disp(sprintf(" Mesh info:\n     Number of (triangular) elements:              %d\n     Total number of nodes:                        %d\n     Number of DOFs(of which Neumann's):           %d(%d)\n",...
                Nt, N, Nh(m+1), Nn));
            % build the linear system terms
            build(Nh(m+1), N);
            % assembling the linear system associated to the discrete grid
            assemble(Nh, Nd, scheme,dt);
            % solving the linear system with ad-hoc algorithm
            u_tmp = solve();
            % reconstruct the whole solution by including the (known) values at the boundaries
            u_h = expand(u_tmp, N);
            % estimate error (u_h-u) in different norms
            [u,E_inf(m+1),E_L2(m+1),E_H1(m+1)] = estimation(u_h, N);
            % plot the reconstructed (discrete) solution
            multiplot(time, multi, Ns, m, u_h, u);
        end
        if multi=='Y'
            % estimate and plot the errors for the grid convergence
            errors(h, Ar, Nh, E_inf, E_L2, E_H1);
        end
    else
        % assign flag variable for dummy t
        t = 1;
        % print information on time-scheme
        disp(sprintf(" Unsteady simulation selected: including time-dependent schemes \n ******** using %s scheme\n",scheme));
        % define the error vector for time-convergence order estimation
        E_t = zeros(10,1);
        % derive a discrete grid object on top of the continuos domain
        [h,Ar,Nt] = discretize(Ns, subspace);
        % extract the number of DOFs and BCs out of all the nodes
        [Ni,Nd,Nn] = reduce();
        % number of DOFs
        Nh = Ni + Nn;
        % number of total nodes
        N = Nh + Nd;
        % define and initialise the IC vector to be used also for time-step progression
        global u_h_prev_t;
        u_h_prev_t = zeros(Nh,1);
        assignIC(N);
        % build the linear system terms
        build(Nh, N);
        % iterate over discretisation steps for time-convergence
        for k=0:0%10
            % compute number of time-steps at k-th iteration
            steps = 100 + k*30;
            % print informations at screen in runtime
            disp(sprintf("\n\n Time-convergence iteration: %d\n Number of time-steps: %d ", k+1, steps));
            % compute the iteration steps to which the solution is printed
            print_iter = ceil(steps/8);
            print_counter = 0;
            % compute the time increment
            dt = time/steps;
            disp(sprintf("\n ******** Starting time-loop ******** \n"));
            % iterate over time-loop
            for iter=0:steps
                % compute the t-th time-step and update t
                t = iter*dt;
                disp(sprintf("    Time-step: %f", t));
                % assembling the linear system associated to the discrete grid
                assemble(Nh, Nd, scheme, dt);
                % solving the linear system at t-th time-step
                u_tmp = solve();
                u_h_prev_t = u_tmp;
                % reconstruct the whole solution by including the (known) values at the boundaries
                u_h = expand(u_tmp, N);
                % estimate error (u_h-u) in different norms
                [u,e_inf,e_L2,e_H1] = estimation(u_h, N);
                % plot solution only at certain iterations to not saturate memory
                if rem(iter,print_iter) == 0
                    print_counter = print_counter + 1;
                    multiplot(time,'N', 8, print_counter, u_h, u);
                end
            end
            newline;
        end
    end
end



