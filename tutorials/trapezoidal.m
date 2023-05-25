% prepare the simulation
clear all
clc
disp(sprintf("*******************************************\n****************** DiSMA ******************"));
disp("Differential Solver for Matlab Applications");
disp("********* an open-source project **********")
disp(sprintf("*******************************************"));
disp(sprintf("\n\n__________________________________________"))
disp(sprintf("\n            SIMULATION STARTING"))
disp(sprintf("__________________________________________"))

% include the triangulator and solver library functions
disp(sprintf("\n\n 1-dimensional simulation selected: adding FEM1D to path"));
addpath('../FEM1D')

% computational domain
a = -1;
b = 1;

% equation's parameters
mu = @(x) 1;
beta = @(x) 1;
sigma = @(x) 1;

% extact solution and derivatives
U = @(x) 1 - x^2;
U1 = @(x) -2*x;
U2 = @(x) -2;
X = linspace(a,b,1000);
u = zeros(1000,1);
for n=1:1000
  u(n) = U(X(n));
endfor

% forcing function
f = @(x) mu(x)*U2(x) + beta(x)*U1(x) + sigma(x)*U(x);

% boundary conditions
Dirichlet = [U(a), U(b)]; % label = +1
Neumann = [U1(a), U1(b)]; % label = -1

% mesh's parameters
N_nodes = 20;
N_dofs = N_nodes -2;
N_elements = N_nodes -1;
size = (b-a)/N_elements;

% compute mesh's nodes
mesh = linspace(a,b,N_nodes);
labels = zeros(1,N_nodes);
% assign DOFs labels
labels(1) = 1;
labels(N_nodes) = 1;

% compute mesh's elements
eles = zeros(2,N_elements);
for j=1:N_elements
  eles(:,j)= [j, j+1];
endfor

% assemble the linear system
A = zeros(N_nodes,N_nodes);
s = zeros(N_nodes,1);
uh = zeros(N_nodes,1);

% loop over each element
for e=1:N_elements
  % loop over each node (1st idx)
  for j=1:2
    j_glb = eles(j,e);
    x = mesh(j_glb);
    % check the label of the node
    if labels(j_glb)==0
      % DOF branch
      for k=1:2
        k_glb = eles(k,e);
        if j_glb==k_glb % main diagonal
          A(j_glb,k_glb) = -2*mu(x)/size + ...     diffusion
                            0 + ...                convection
                            (2/3)*sigma(x)*size; % rection
        elseif k_glb==j_glb-1 % lower diagonal
          A(j_glb,k_glb) = mu(x)/size + ...        diffusion
                          -0.5*beta(x) + ...       convection
                          (size*sigma(x))/6;     % reaction
        elseif k_glb==j_glb+1 % upper diagonal
          A(j_glb,k_glb) = mu(x)/size + ...        diffusion
                           0.5*beta(x) + ...        convection
                           (size*sigma(x))/6;     % reaction
        endif
      endfor
      s(j_glb) = f(x)*size;
    else
      % boundary branch
      if j_glb==1
        % left (lower) boundary
        if labels(j_glb)>0
          % Dirichlet branch
          A(j_glb, j_glb) = 1;
          A(j_glb + 1, j_glb) = mu(x)/size - 0.5*beta(x) + (size*sigma(x))/6;
          s(j_glb) = Dirichlet(1);
        else
          % Neumann branch
          A(j_glb, j_glb) = -mu(x)/size - 0.5*beta(x) + (1/3)*sigma(x)*size;
          A(j_glb + 1, j_glb) = mu(x)/size - 0.5*beta(x) + (size*sigma(x))/6;
          A(j_glb, j_glb + 1) = mu(x)/size + 0.5*beta(x) + (size*sigma(x))/6;
          s(j_glb) = f(x)*size/2 - mu(x)*Neumann(1);
        endif
      else
        % right (upper) boundary
        if labels(j_glb)>0
          % Dirichlet branch
          A(j_glb, j_glb) = 1;
          A(j_glb - 1, j_glb) = mu(x)/size + 0.5*beta(x) + (size*sigma(x))/6;
          s(j_glb) = Dirichlet(2);
        else
          % Neumann branch
          A(j_glb, j_glb) = -mu(x)/size + 0.5*beta(x) + (1/3)*sigma(x)*size;
          A(j_glb - 1, j_glb) = mu(x)/size - 0.5*beta(x) + (size*sigma(x))/6;
          A(j_glb, j_glb - 1) = mu(x)/size + 0.5*beta(x) + (size*sigma(x))/6;
          s(j_glb) = f(x)*size/2 - mu(x)*Neumann(2);
        endif
      endif
    endif
  endfor
endfor

% solve and plot the numerical solution
uh = A\s;
plot(mesh(1,:),uh,'r--o','LineWidth',2, 'DisplayName', 'u_h(x)')
grid on
hold on
plot(X,u,'b','LineWidth',1.5, 'DisplayName', '1-x^2')
xlabel('x', 'interpreter','latex')
ylabel('y', 'interpreter','latex')
legend('FontSize', 16)
title('N=20 nodes', 'interpreter','latex')
