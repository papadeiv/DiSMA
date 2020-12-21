%{
******************************************************************************************************
                                      HOMOGENEOUS DIRICHLET BCs
******************************************************************************************************
%}

% prepare the simulation
clear all
close all
clc

% including the triangulator library and solver library functions according to the choice of the subspace
addpath('bbtr30')
disp('bbtr30 added to the path')
addpath('P1elements')
disp(fprintf('P1elements added to the path'));

% define coefficient functions of the PDE
global coefficient_functions;
coefficient_functions = {
    @(x,y) 1.0;...  % diffusion function (nu)
    @(x,y) 0.0;...  % convection function (beta)
    @(x,y) 0.0};    % reaction function (sigma)

% define analytical solution and its derivatives
global exact_functions;
exact_functions = {
    @(x,y) 16*x*(1-x)*y*(1-y);...     % exact solution (U)
    @(x,y) 16*(2*x-1)*(y-1)*y;...     % x-partial derivative (U_x)};
    @(x,y) 16*(2*y-1)*(x-1)*x};       % y-partial derivative (U_y)};

% define source function
global source_function;
source_function = @(x,y) 32*(x*(1-x) + y*(1-y));

% construct boundary functions cell array
global boundary_functions;
boundary_functions = {
    % first row -> Dirichlet BCs
    @(x,y) 0.0,...      % border 1
    @(x,y) 0.0,...      % border 2
    @(x,y) 0.0,...      % border 3
    @(x,y) 0.0;         % border 4
    % second row -> Neumann BCs
    @(x,y) 0.0,...      % border 1
    @(x,y) 0.0,...      % border 2
    @(x,y) 0.0,...      % border 3
    @(x,y) 0.0};        % border 4

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
inputs = boundaries;

% define number of iterations for convergence (mesh refinement)
Ns = 7;

% launch the simulation
main(Ns);