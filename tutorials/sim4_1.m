
% prepare the simulation
clear all
clc

% include the triangulator and solver library functions
addpath('../FEM2D')
disp('FEM2D added to the path')
addpath('../FEM2D/bbtr30')
disp('bbtr30 added to the path')

% define coefficient functions of the PDE
global coefficient_functions;
coefficient_functions = {
    @(x,y) 1.0;...  % diffusion function (nu)
    @(x,y) 0.0;...  % convection function (beta)
    @(x,y) 0.0};    % reaction function (sigma)

% define analytical solution and its derivatives
global exact_functions;
exact_functions = {
    @(x,y) exp(-x^2 - y^2);...          % exact solution (U)
    @(x,y) -2*x*exp(-x^2 - y^2);...     % x-partial derivative (U_x)};
    @(x,y) -2*y*exp(-x^2 - y^2)};       % y-partial derivative (U_y)};

% define source function
global source_function;
source_function = @(x,y) -coefficient_functions{1}(x,y)*4*(x^2 + y^2 -1)*exact_functions{1}(x,y);

% construct boundary functions cell array
global boundary_functions;
boundary_functions = {
    % first row -> Dirichlet BCs
    @(x,y) exp(-x^2),...      % border 1
    @(x,y) exp(-(y^2+1)),...  % border 2
    @(x,y) exp(-(x^2+1)),...  % border 3
    @(x,y) exp(-y^2);         % border 4
    % second row -> Neumann BCs
    @(x,y) 2*y*exact_functions{1}(x,y),...              % border 1
    @(x,y) -2*x*exact_functions{1}(x,y),...             % border 2
    @(x,y) -2*y*exact_functions{1}(x,y),...             % border 3
    @(x,y) 2*x*exact_functions{1}(x,y)};                % border 4

% define computational domain 
% BY CONVENTION: border 1 goes from first node (v1) to second node (v2); border 2 goes from (v2) to (v3); border 3 from (v3) to (v4) and so on
v1 = [0, 0];
v2 = [1, 0];
v3 = [1, 1];
v4 = [0, 1];

global vertices;
vertices = [v1; v2; v3; v4];
clear v1 v2 v3 v4;

% imposition of the BCs markers on each boundary of the domain
% BY CONVENTION: odd markers are for Dirichlet BCs while even markers are for Neumann BCs
b1 = 1;
b2 = 2;
b3 = 4;
b4 = 3;

global boundaries;
boundaries = [b1 b2 b3 b4];
clear b1 b2 b3 b4;

% imposition of the BCs markers on each input vertex of the domain
% BY CONVENTION: odd markers are for Dirichlet BCs while even markers are for Neumann BCs
global inputs;
inputs = [5 7 0 9];

% define basis functions subspaces in string array
subspace = ["P1", "P2"];

% define choices for grid-convergent simulations
grid_convergence = ["N", "Y"];

% define choices for levels of increasing grid density [0 -> very coarse mesh; 9 -> very refined mesh]
Ns = [0, 1, 3, 5, 7, 9];

% launch the simulation
main(grid_convergence(1), Ns(5), subspace(1));