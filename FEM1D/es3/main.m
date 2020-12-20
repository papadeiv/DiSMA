clear all
close all
clc
% define computational domain size
L = 2;
% define number of simultaions
Ns = 10;
% define vector of iterative values of apporx. solution in node x=x0
u0 = zeros(Ns,1);
for n=1:Ns
    % create the discrete grid over the computational domain
    [Ih,xh,h,Nh] = discretise(n,L);
    % assembling the linear system
    [A,b] = build(Ih,xh,h,Nh);
    % solve for u_h
    u_h_tmp = A\b;
    % reconstruct the full-order solution by including the (known) boundary values
    [u_h] = expand(u_h_tmp);
    % compute the solutions' first derivatives for error-estimation purposes
    [u_h_x,x_h] = derivate(xh(1,:),h,u_h);
    % extract the approx. solution's value in x=1
    u0(n) = evaluate(u_h,xh(1,:));
    % plotting the approximated solution
    figure(2)
    hold on
    plot(xh(1,:),u_h,'r--o','LineWidth',2)
    xlabel('x')
    ylabel('u(x)')
    title('Approximated solution')
end
figure(1)
hold on
plot(x_h,u_h_x,'b','LineWidth',1.5)
figure(2)
hold on
plot(xh(1,:),u_h,'b','LineWidth',1.5)
steps = linspace(1,Ns,Ns);
figure(3)
plot(steps,u0)
xlabel('Simulation steps')
ylabel('Value of u_h(x) at x=1')
title('Values of u_h(x=1) at different iterations')
