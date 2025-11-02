% This script provides the geometric and material parameters of the shell,
% information about the boundary conditions and creates a model based on
% the shell Theory of Donnell by calling sym_Shell(). Afterwards, the
% external loading is defined and the non-dimensionalized EOM are solved in
% two intervals with different solvers and tolerances.
%% Global constants

global t_old    % needed for the waitbar
global tend

%% Path

path = pwd; % path to project folder

mkdir Elastic_Forces Donnell
addpath(genpath('Elastic_Forces'))
mkdir External_Forces Donnell
addpath(genpath('External_Forces'))
mkdir Jacobians Donnell
addpath(genpath('Jacobians'))
mkdir Kinetic_Energy Donnell
addpath(genpath('Kinetic_Energy'))
mkdir Sim_Res Donnell
addpath(genpath('Sim_Res'))

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
    M = 5; N = 5;       % maximum modal indices in u- and v-direction
    Mw = 3; Nw = 3;     % maximum modal indices in w-direction
    M0 = 1; N0 = 1;     % maximum modal indices of geometric imperfection
else
    M = 3; N = 3;
    Mw = 1; Nw = 1;
    M0 = 1; N0 = 1;
end

n_mod = ceil(M/2)*ceil(N/2); % number of modes (=4 for the 9-dof model and =9 for the 22-dof model)

%% 1.3 Damping coefficients

% natural frequencies for Rx = Ry = 1 in [Hz]
% The eigenfrequencies are needed to calculate the damping coefficients.
ome11 = 952.26; 
ome13 = 2575.9; 
ome33 = 4472.3; 

par.eig = zeros(n_mod,1);
par.eig(1:4) = [ome11 ome13 ome13 ome33]';
par.eig(5:end) = par.eig(4);    % for the high-order modes (e.g. (5,5)) the eigenfrequency ome33 is used

zeta11 = 0.004;
par.zeta = zeta11*ones(n_mod,1);% zeta11 is used for all modes

T_11 = 1/ome11; % the period of the first eigenfrequency is used for non-dimensionalization
par.T = T_11;

%% 1.4 Boundary Conditions
% This decides whether the simplified or the full additional terms for the
% modal extension are taken.

bc = 'simplified';
% bc = 'full';

%% 1.5 Perform Symbolic Computations

tic
shell = sym_Shell(par,bc,n_dof,path); % sym_Don performs the symbolic calculations for Donnell's shell theory
t_sym = toc

%% 2. External loading

% You can choose whether you want a ramping or a constant normal force.
amplitude = 'ramp';
% amplitude = 'const';

% excitation force
f_start = 0;    % initial force amplitude
df = 1;         % step size
f_end = 30;     % final value

t_start = 0;    % starting time of ramp
t_end = 30;     % finishing time of ramp

tau_0 = floor(t_start/T_11); % non-dimensionalized starting time of ramp
tau_end = floor(t_end/T_11); % non-dimensionalized finishing time of ramp


if strcmp(amplitude,'ramp') % this is used for saving the simulation results
    f_N = f_end;
else
    f_N = 0;    % change for constant force
end

freq_N = 0.8*ome11; % frequency of normal force in [Hz]
ome_N = 2*pi*freq_N;% frequency of normal force in [rad/s]

%% 3. Numerical Solution
% The following sections only treat the numerical solution. The modelling
% part is done.
%% 3.1 Non Dimensionalized
%% 3.1.1 Waitbar (got this from Ludwig Krumm, Institute of Mechanics and Ocean Engineering, TUHH)

%Die Funktionen "waitbar" und "solver_output" generiern einen Fortschritts-
%balken, der den Fortschritt des Integrators in Relation zur Endzeit der
%Simulation zeigt und nach Ende der Simulation geschlossen wird.

t0 = 0;        % initial time
tend = 10;     % duration of simulation in [s]
Tau_end = floor(tend/T_11); % non-dimensionalized final time
tend = Tau_end;

spc = 50;       % number of samples per cycle
dtau = 1/spc;

tau_span= [0:dtau:Tau_end]; % non-dimensionalized time vector

%% 3.1.2 Numerical Integration of non-dimensional EOM

% Initial Conditions
x0 = zeros(n_dof,1);
w0 = -[1.8 -0.1 -0.1 0.04]';   % values are taken from the time series in the paper; multiples of shell thickness
x0(2*n_mod+1:n_dof) = w0(1:ceil(Mw/2)*ceil(Nw/2));
v0 = zeros(n_dof,1);
y0 = [x0;v0];

%% First Interval

Tau_1 = floor(5/T_11);
tend = Tau_1;
tau1 = [0:dtau:Tau_1];

waitbar_instance=waitbar(0,'Simulation at tau = 0','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
cleanup=onCleanup(@()delete(waitbar_instance));
t_old = -1;
tic;

solver_output_fcn=@(t,Y,flag) solver_output(t,Y,flag,waitbar_instance,cleanup);

if strcmp(amplitude,'ramp') % different input arguments to the Jacobian and EOM are required for ramp and constant force
    opt23s = odeset('Jacobian',@(t,x) shell.Jf_nd(x,ramp_fun(f_start,df,f_N,tau_0,tau_end,t),ome_N,t),'OutputFcn',solver_output_fcn,'RelTol',1e-3,'AbsTol',1e-6,'Stats','On');
    [T1,Y1] = ode23s(@(t,x) EOM(t,x,ramp_fun(f_start,df,f_N,tau_0,tau_end,t),ome_N,2*n_dof,shell.dU_nd,shell.Q_nd),tau1,y0,opt23s);
else
    opt23s = odeset('Jacobian',@(t,x) shell.Jf_nd(x,f_N,ome_N,t),'OutputFcn',solver_output_fcn,'RelTol',1e-3,'AbsTol',1e-6,'Stats','On');
    [T1,Y1] = ode23s(@(t,x) EOM(t,x,f_N,ome_N,2*n_dof,shell.dU_nd,shell.Q_nd),tau1,y0,opt23s);
end

%% Second Interval

tend = Tau_end;
tau2 = [Tau_1:dtau:Tau_end];
y1 = Y1(end,:)';

waitbar_instance=waitbar(0,'Simulation at tau = 0','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
cleanup=onCleanup(@()delete(waitbar_instance));
t_old = -1;
tic;

solver_output_fcn=@(t,Y,flag) solver_output(t,Y,flag,waitbar_instance,cleanup);

if strcmp(amplitude,'ramp') % different input arguments to the Jacobian and EOM are required for ramp and constant force
    opt15s = odeset('Jacobian',@(t,x) shell.Jf_nd(x,ramp_fun(f_start,df,f_N,tau_0,tau_end,t),ome_N,t),'OutputFcn',solver_output_fcn,'RelTol',1e-6,'AbsTol',1e-12,'BDF','on','MaxOrder',5,'Stats','on');
    [T2,Y2] = ode15s(@(t,x) EOM(t,x,ramp_fun(f_start,df,f_N,tau_0,tau_end,t),ome_N,2*n_dof,shell.dU_nd,shell.Q_nd),tau2,y1,opt15s);
else
    opt15s = odeset('Jacobian',@(t,x) shell.Jf_nd(x,f_N,ome_N,t),'OutputFcn',solver_output_fcn,'RelTol',1e-6,'AbsTol',1e-12,'BDF','on','MaxOrder',5,'Stats','on');
    [T2,Y2] = ode15s(@(t,x) EOM(t,x,f_N,ome_N,2*n_dof,shell.dU_nd,shell.Q_nd),tau2,y1,opt15s);
end

%% Combine Results from both Intervals

T = [T1(1:end-1,1);T2];
Y = [Y1(1:end-1,:);Y2];

% Save Simulation Results
Res_path = strcat(path,'\Sim_Res\Donnell\');
Y_filename = sprintf(['Y_nd_%ddof_fn%d_' amplitude '_omef%d_10w0%d_1000z11%d_' bc '_%dp.mat'],n_dof,floor(f_N),floor(freq_N),abs(10*w0(1)),1000*zeta11,Tau_end);
proceed = 'N';

while(isequal(proceed,'N'))
    if (exist(Y_filename) == 2)
        OW = input('Resolution file with this name already exists. Shall I overwrite it? (Y/N)','s');
        if isequal(OW,'Y')
            save(strcat(Res_path,Y_filename),'T','Y','bc','solver','t_sim');
            disp('Existing file was overwritten.')
        else
            disp('Existing file was NOT overwritten.')
        end       
    else
        save(strcat(Res_path,Y_filename),'T','Y','bc','solver','t_sim','t_start','t_end','f_start','f_end');
        disp('Solution saved to file.')
    end
    proceed = input('Do you want to proceed? Y/N','s');
end
proceed = 'N';

%% 3.1.3 Overview

figure(1)           
plot(T,Y(:,2*n_mod+1))  % time series of w11
xlabel('$\tau$','Interpreter','latex')
ylabel(['$' cast(shell.q(2*n_mod+1),'char') '/h$'],'Interpreter','latex')


figure(2)           
plot(T,Y(:,1))  % time series of u11
xlabel('$\tau$','Interpreter','latex')
ylabel(['$' cast(shell.q(1),'char') '/h$'],'Interpreter','latex')


%% 3.1.4 Re-Dimensionalization

% t_skip = find(T > 0.5,1,'first')+1;
t_skip = max([length(T)-500000,0]);     % only keep last 500000 samples
Y = Y(1+t_skip:end,:);
Y(:,n_dof+1:end) = 1/(2*pi)*Y(:,n_dof+1:end);
T = T_11*T(1+t_skip:end);