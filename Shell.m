%% General Loading
% This is a more streamlined version of "main_shell_v2.m". It only performs
% the symbolic computations.

clc; clearvars;
diary('Log.txt')
diary on

%% Path
path = 'E:\Dokumente\Studium\Master\4th_Semester\Projektarbeit_Strukturdynamik\Matlab\'; % path to project folder

%% 1. Symbolic Computations
%% 1.1 Parameter Values CORRECT
par = struct;
par.a = 0.1;    %width in [m]
par.b = 0.1;    %length in [m]
par.h = 0.001;  %thickness in [m]
par.Rx = 1;     %radius of curvature in [m]
par.Ry = 1;
par.E = 206e9;   %Young's modulus in [Pa = N/m^2] %FIXED
par.rho = 7800;  %mass density in [kg/m^3]
par.nu = 0.3;    %Poisson's ration
par.A11 = 0; %geometric imperfection

%% 1.2 DoF
% n_dof = 22;
n_dof = 9;

if n_dof == 22
    M = 5; N = 5;
    Mw = 3; Nw = 3;
    M0 = 1; N0 = 1;
else
    M = 3; N = 3;
    Mw = 1; Nw = 1;
    M0 = 1; N0 = 1;
end

n_modes = ceil(M/2)*ceil(N/2);

%% 1.3 Damping coefficients

% natural frequencies for Rx = Ry = 1
ome11 = 952.26; 
ome13 = 2575.9; 
ome33 = 4472.3; 

par.eig = [ome11; ome13; ome33];

par.zeta = 0.004;

T11 = 1/ome11;
par.T = T11;

%% 1.5 Perform Symbolic Calculations

tic
shell = sym_ama(par,bc,n_dof,path);
t_sym = toc