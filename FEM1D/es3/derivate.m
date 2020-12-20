function [u_h_x,x_h] = derivate(xh,h,u_h)

    % define the nodes in which the derivative will be computed (midpoint
    % of each interval)
    N_x = length(xh)-1;
    x_h = zeros(1,N_x);
    u_h_x = zeros(1,N_x);
    for j=1:N_x
        x_h(j)=xh(j)+h/2;
        % compute approximated solution's first derivative
        u_h_x(j) = (u_h(j+1)-u_h(j))/h;
    end
    figure(1)
    hold on
    plot(x_h,u_h_x,'r--o','LineWidth',2)
    xlabel('x')
    ylabel('u_x(x)')
    title('Approximated derivative')
    
end