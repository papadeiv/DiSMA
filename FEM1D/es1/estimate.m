function [e_inf,e_L2] = estimate(h,u,u_x,u_h,u_h_x)
    
    % compute the error-norms between solutions
    e_inf = zeros(2,1);
    e_L2 = zeros(2,1);
    diff = u-u_h;
    e_inf(1,1) = max(abs(diff));
    e_L2(1,1) = sqrt(h*sum(diff.^2));
    % compute the error-norms between first derivatives of the solutions
    clear diff;
    diff = u_x - u_h_x;
    e_inf(2,1) = max(abs(diff));
    e_L2(2,1) = sqrt(h*sum(diff.^2));
    
end