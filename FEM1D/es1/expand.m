function [u,u_h] = expand(u_tmp,u_h_tmp)

    Nh = length(u_tmp);
    N = Nh + 2;
    u = zeros(N,1);
    u_x = u;
    u_h = u_x;
    % construct the exact solution vector
    u(1) = 0;
    u(2:N-1) = u_tmp;
    u(N) = u(1);
    % construct the apporximated solution vector
    u_h(1) = 0;
    u_h(2:N-1) = u_h_tmp;
    u_h(N) = u_h(1);
    
end