% clean stuff
clear all
close all
clc

% computational domain
a = % put the left/lower bound here
b = % put the right/upper bound here

% equation's parameters
mu = @(x) % put the diffusion function here
beta = @(x) % put the advection function here
sigma = @(x) % put the reaction function here

% extact solution and derivatives
U = @(x) % put the exact solution here
U1 = @(x) % put the solution's first derivative here
U2 = @(x) % put the solution's second derivative here

% forcing function
f = @(x) (nu^2)*besselj(nu, x);

% boundary conditions
left =  % Dirichlet = +1; Neumann = -1
right = % Dirichlet = +1; Neumann = -1

% quadrature samples in [0,1]
N_samples = % any natural number between 1 to 40
[nodes, weights] = samples(N_samples, 0); % G-L = 0; MTQR = 1

% boundary conditions
Dirichlet = [U(a), U(b)];
Neumann = [U1(a), U1(b)];

% mesh parametrisation
N_nodes = % any natural number greater than 1
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
      % compute LHS integral via Gaussian quadrature
      quadrature = 0.0;
      if j==1 % phi_k(x) = 1 - x
        for n=1:N_samples
          quadrature = quadrature + weights(n)*f(size*nodes(n) + mesh(eles(1,e)))*phi(nodes(n));
        endfor
      else % phi_k(x) = x
        for n=1:N_samples
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
ylim manual
ylim([-1.5 5])
grid on
hold on
plot(mesh(1,:),uh,'r--o','LineWidth',2, 'DisplayName', 'u_h(x)')
% compute and plot the exact solution
X = linspace(a,b,1000);
% IF PLOT OF THE EXACT SOLUTION DOESN'T SHOW UNCOMMENT THE LINES BELOW
%u = zeros(1000,1);
%for n=1:1000
%  u(n) = U(X(n));
%endfor
%plot(X, u,'b','LineWidth',2, 'DisplayName', 'u(x)')
plot(X, U(X),'b','LineWidth',2, 'DisplayName', 'u(x)')
xlabel('x', 'interpreter','latex')
ylabel('y', 'interpreter','latex')
legend('FontSize', 16)
title('Piecewise Linear Finite Element Solution', 'interpreter','latex')
hold off

% compute error in L2-norm
N_dofs = N_nodes - 2;
error_vector = zeros(N_dofs, 1);
for n=1:N_dofs
  error_vector(n) = uh(n+1) - U(mesh(n+1));
endfor
error = norm(error_vector, 2);
