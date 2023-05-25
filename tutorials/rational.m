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
a = -2;
b = 1;

% equation's parameters
mu = @(x) 1;
beta = @(x) 0;
sigma = @(x) 1;

% extact solution and derivatives
U = @(x) (x^3 - 1)/(7 - x^2);
U1 = @(x) (3*x^2*(7 - x^2) - (x^3 - 1)*(-2*x))/(7 - x^2)^2;
U2 = @(x) 1;

% forcing function
f = @(x) (63 - x*(294 + x*(8 + x*(63 + x*(x^3 - 14*x - 1)) )))/(x^2 - 7)^3;

% boundary conditions
left = +1;  % Dirichlet = +1, Neumann = -1
right = +1; % Dirichlet = +1, Neumann = -1

% quadrature samples in [0,1]
[nodes, weights] = samples(3, 0); % MTQR quadrature samples for the first element

% boundary conditions
Dirichlet = [U(a), U(b)];
Neumann = [U1(a), U1(b)];

% mesh parametrisation
N_nodes = 15;
lambda = 1; % 0 -> no refinement; 1 -> refinement @ x=a; 2 -> refinement @ x=b

% create the mesh
if lambda==0 % no refinement: uniform partition
  mesh = linspace(a,b,N_nodes);
  size = (b-a)/N_elements;
elseif lambda==1 % refinement around left/lower bound
  mesh = linspace(a,b,N_nodes);
elseif lambda==2 % refinement around right/upper bound
  mesh = linspace(a,b,N_nodes);
endif

% assign nodes labels
labels = zeros(1,N_nodes);
labels(1) = left;
labels(N_nodes) = right;

% build elements datastructure
N_elements = N_nodes -1;
eles = zeros(3,N_elements);
for j=1:N_elements
  eles(:,j)= [j, j+1, mesh(j+1)-mesh(j)];
endfor

% define (piecewise) linear basis function
phi = @(x) 1 - x;
psi = @(x) x;

% initialise the discretisation
A = zeros(N_nodes,N_nodes);
s = zeros(N_nodes,1);
uh = zeros(N_nodes,1);

% assemble the linear system
for e=1:N_elements
  size = eles(3,e);
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
      % compute LHS integral of the element via the specified quadrature
      quadrature = 0.0;
      if j==1 % phi_k(x) = 1 - x
        for n=1:length(nodes)
          quadrature = quadrature + weights(n)*f(size*nodes(n) + mesh(eles(1,e)))*phi(nodes(n));
        endfor
      else % phi_k(x) = x
        for n=1:length(nodes)
          quadrature = quadrature + weights(n)*f(size*nodes(n) + mesh(eles(1,e)))*psi(nodes(n));
        endfor
      endif
      % update the source term
      s(j_glb) = s(j_glb) + size*quadrature;
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
grid on
hold on
plot(mesh(1,:),uh,'r--o','LineWidth',2, 'DisplayName', 'u_h(x)')
% compute and plot the exact solution
X = linspace(a,b,1000);
u = zeros(1000,1);
for n=1:1000
  u(n) = U(X(n));
endfor
plot(X, u,'b','LineWidth',2, 'DisplayName', '(x^3 - 1)/(7 - x^2)')
xlabel('x', 'interpreter','latex')
ylabel('y', 'interpreter','latex')
legend('FontSize', 16)
title('N=15 nodes', 'interpreter','latex')
hold off

% compute error in L2-norm
N_dofs = N_nodes - 2;
error_vector = zeros(N_dofs, 1);
for n=1:N_dofs
  error_vector(n) = uh(n+1) - U(mesh(n+1));
endfor
error = norm(error_vector, 2);
