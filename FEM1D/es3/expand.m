function [u_h] = expand(u_h_tmp)

    Nh = length(u_h_tmp);
    N = Nh + 2;
    u_h = zeros(N,1);
    % construct the apporximated solution vector
    u_h(1) = 2;
    u_h(2:N-1) = u_h_tmp;
    u_h(N) = 5;
    
end