
% prepare the simulation
clear all
clc
disp(sprintf("*******************************************\n****************** DiSMA ******************"));
disp("Differential Solver for Matlab Applications");
disp("********* an open-source project **********")
disp(sprintf("*******************************************"));
disp(sprintf("\n\n __________________________________________"))
disp(sprintf("\n            SIMULATION STARTING"))
disp(sprintf(" __________________________________________"))

% include the triangulator and solver library functions
disp(sprintf("\n\n 2-dimensional simulation selected: adding FEM2D to path"));
addpath('../FEM2D')
disp(sprintf(' ******** FEM2D added to the path'));
disp(sprintf('\n Triangular elements selected: adding bbtr30 added to the path'));
addpath('../FEM2D/bbtr30')
disp(sprintf(' ******** bbtr30 added to the path\n\n'));

% define coefficient functions of the PDE
global coefficient_functions;
coefficient_functions = {
    @(x,y,t) 1.0;...  % diffusion function (nu)
    @(x,y,t) 0.0;...  % convection function (beta)
    @(x,y,t) 0.0};    % reaction function (sigma)

% define analytical solution and its derivatives
global exact_functions;
exact_functions = {
    @(x,y,t) exp(-x^2 - y^2);...          % exact solution (U)
    @(x,y,t) -2*x*exp(-x^2 - y^2);...     % x-partial derivative (U_x)};
    @(x,y,t) -2*y*exp(-x^2 - y^2)};       % y-partial derivative (U_y)};

% define source function
global source_function;
source_function = @(x,y,t) -coefficient_functions{1}(x,y,t)*4*(x^2 + y^2 -1)*exact_functions{1}(x,y,t);

% construct boundary functions cell array
global boundary_functions;
boundary_functions = {
    % first row -> Dirichlet BCs
    @(x,y,t) exp(-x^2),...      % border 1
    @(x,y,t) exp(-(y^2+1)),...  % border 2
    @(x,y,t) exp(-(x^2+1)),...  % border 3
    @(x,y,t) exp(-y^2);         % border 4
    % second row -> Neumann BCs
    @(x,y,t) 2*y*exact_functions{1}(x,y),...              % border 1
    @(x,y,t) -2*x*exact_functions{1}(x,y),...             % border 2
    @(x,y,t) -2*y*exact_functions{1}(x,y),...             % border 3
    @(x,y,t) 2*x*exact_functions{1}(x,y)};                % border 4

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
b2 = 3;
b3 = 5;
b4 = 7;

global boundaries;
boundaries = [b1 b2 b3 b4];
clear b1 b2 b3 b4;

% imposition of the BCs markers on each input vertex of the domain
% BY CONVENTION: odd markers are for Dirichlet BCs while even markers are for Neumann BCs
global inputs;
inputs = [9 11 13 15];

% define choices for time-dependence
time = [0, pi/2];

% define choices for time-discretisation schemes
time_scheme = ["1st-order", "2nd-order"];

% define choices for grid-convergent simulations in string array
grid_convergence = ["N", "Y"];

% define choices for basis functions subspaces in string array
subspace = ["P1", "P2"];

% define choices for levels of increasing grid density [0 -> very coarse mesh; 9 -> very refined mesh]
Ns = [0, 1, 3, 5, 7, 9];

% launch the simulation
main(time(1), grid_convergence(1), subspace(2), time_scheme(1), Ns(5));