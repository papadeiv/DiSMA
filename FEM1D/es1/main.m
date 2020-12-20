clear all
close all
clc
% define computational domain size
L = pi;
% define number of simultaions
Ns = 4;
% initialise the error vectors
E_inf = zeros(2,Ns);
E_L2 = zeros(2,Ns);
E_H1 = zeros(2,Ns);
for n=1:Ns
    % create the discrete grid over the computational domain
    [Ih,xh,h,Nh] = discretise(n,L);
    % assembling the linear system
    [A,b,u_tmp] = build(Ih,xh,h,Nh);
    % solve for u_h
    u_h_tmp = A\b;
    % reconstruct the full-order solution by including the (known) boundary values
    [u,u_h] = expand(u_tmp,u_h_tmp);
    % compute the solutions' first derivatives for error-estimation purposes
    [u_h_x,u_x,x_h] = derivate(xh(1,:),h,u_h);
    % error estimation in different norms
    [E_inf(:,n),E_L2(:,n)] = estimate(h,u,u_x,u_h,u_h_x);
    % plotting the approximated solution
    figure(2)
    hold on
    plot(xh(1,:),u_h,'r--o','LineWidth',2)
    legend('u_{h}')
end
% plotting the analytical first derivative over the derivative of the approximated solution
figure(1)
hold on
plot(x_h,u_x,'b','LineWidth',1.5,'DisplayName','u_x')
xlabel('x')
ylabel('u_x(x)')
legend
title('Approximated (first) derivative')
% plotting the analytical exact solution over the approximated solution
figure(2)
hold on
plot(xh(1,:),u,'b','LineWidth',1.5,'DisplayName','u')
xlabel('x')
ylabel('u(x)')
legend
title('Approximated solution')
% plotting the error-convergence of the solution
steps = linspace(1,Ns,Ns);
figure(3)
plot(steps,E_inf(1,:),'r-*',steps,E_L2(1,:),'b-o')
xlabel('Simulation step')
ylabel('Error')
legend('L_{\infty} - norm','L_2 - norm')
title('Solution error estimation')
% plotting the error-convergence of the first derivatives
figure(4)
plot(steps,E_inf(2,:),'r-*',steps,E_L2(2,:),'b-o')
xlabel('Simulation step')
ylabel('Error')
legend('L_{\infty} - norm','L_2 - norm')
title('First derivatives error estimation')