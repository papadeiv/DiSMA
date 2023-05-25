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
a = 0;
b = 3*pi;

% equation's parameters
mu = @(x) x^2;
beta = @(x) x;
sigma = @(x) x^2;

% order of Bessel's function
nu = 0.5;

% extact solution and derivatives
U = @(x) besselj(nu, x);
U1 = @(x) 1;
U2 = @(x) 1;

% forcing function
f = @(x) (nu^2)*U(x);

% boundary conditions
left = +1;  % Dirichlet = +1, Neumann = -1
right = +1; % Dirichlet = +1, Neumann = -1

% quadrature samples in [0,1]
[mtqr_nodes, mtqr_weights] = samples(3, 0); % MTQR quadrature samples for the first element
[gl_nodes, gl_weights] = samples(3, 0); % G-L quadrature samples for the last element

% boundary conditions
Dirichlet = [U(a), U(b)];
Neumann = [U1(a), U1(b)];

% mesh parametrisation
N_nodes = 100;
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
      % initialise LHS contribution
      quadrature = 0.0;
      if e==1
        % if the element is the first one than use MTQR samples
        sampled_x = mtqr_nodes;
        sampled_w = mtqr_weights;
      else
         % otherwise use G-L quadrature samples for all the remaining elements
         sampled_x = gl_nodes;
         sampled_w = gl_weights;
      endif
      % compute LHS integral of the element via the specified quadrature
      if j==1 % phi_k(x) = 1 - x
        for n=1:length(sampled_x)
          quadrature = quadrature + sampled_w(n)*f(size*sampled_x(n) + mesh(eles(1,e)))*phi(sampled_x(n));
        endfor
      else % phi_k(x) = x
        for n=1:length(sampled_x)
          quadrature = quadrature + sampled_w(n)*f(size*sampled_x(n) + mesh(eles(1,e)))*psi(sampled_x(n));
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
plot(X, U(X),'b','LineWidth',2, 'DisplayName', 'J_v (x), v=1/2')
xlabel('x', 'interpreter','latex')
ylabel('y', 'interpreter','latex')
legend('FontSize', 16)
title('N=100 nodes (mixed MTQR-GL quadrature)', 'interpreter','latex')
hold off

% compute error in L2-norm
N_dofs = N_nodes - 2;
error_vector = zeros(N_dofs, 1);
for n=1:N_dofs
  error_vector(n) = abs(uh(n+1) - U(mesh(n+1)));
endfor
error_L2 = norm(error_vector, 2);
error_Li = norm(error_vector, Inf);
