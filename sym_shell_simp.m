function [shell] = sym_shell_simp(par,bc,n_dof,path)
% This version of the symbolic computations is the most recent one for the
% model without inertial coupling. It uses
% less simplifications and also computes the force and moment resultants at
% the boundaries to check satisfaction of the boundary conditions.

shell = struct;
% This function will store results from intensive symbolic integration and
% symbolic functions under the respective paths.

Us_path = strcat(path,'Elastic_Forces\Donnell\'); % path to Elastic_Forces folder
Q_path = strcat(path,'External_Forces\Donnell\');
T_path = strcat(path,'Kinetic_Energy\Donnell\');
J_path = strcat(path,'Jacobians\Donnell\');
EOM_path = strcat(path,'EOM\Donnell\');

%% Symbolic Parameters
Rx = sym('Rx'); Ry = sym('Ry');                    % radii of curvature
t = sym('t');                                      % time
x = sym('x');  y = sym('y');                       % coordinates along the edges
psi = sym('psi'); theta = sym('theta');            % angular coordinates
a = sym('a'); b = sym('b'); h = sym('h');          % curvilinear dimensions of the shell
E = sym('E'); nu = sym('nu'); rho = sym('rho');    % material properties (Young's modulus, Poisson's ratio and density)

% x = Rx * psi and y = Ry * theta
%% DoFs

if n_dof == 22
    M = 5; N = 5;
    Mw = 3; Nw = 3;
    M0 = 1; N0 = 1;
else
    M = 3; N = 3;
    Mw = 1; Nw = 1;
    M0 = 1; N0 = 1;
end

% Modal extension
U = sym('u%d%d', [M N]);
V = sym('v%d%d', [M N]);
W = sym('w%d%d', [Mw Nw]);      % modal amplitudes (time-dependent)
W0 = sym('A%d%d', [M0 N0]);     % geometric imperfections

% U = [ u11 u12 ... ]           % the entries with even indices are not
%     [ u21 u22 ... ]             needed
%     [ u31 u32 u33 ]

dU = sym('du%d%d', [M N]);
dV = sym('dv%d%d', [M N]);
dW = sym('dw%d%d', [Mw Nw]);    % modal velocities

u = 0; v = 0; w = 0; w0 = 0;


% modal extension for simply supported edges
% The formulation might seem unneccessarily complicated. This was done
% to bring the dofs into a nice order for plotting. The vector of modal
% amplitudes will look like: q = [u11 v11 u13 v13 ... u55 v55 w11 ...
% w33]'
    
q = sym('q%d', [n_dof 1]);
dq = sym('dq%d',[n_dof 1]);     % vectors of dofs and their derivatives
    
k = 2*ceil(Mw/2)*ceil(Nw/2) + 1;   % number of modes in w-direction plus 1
j = 1; i = 1;
 for m = 1:2:M
     for n = 1:2:N
         if n <= Nw && m <= Mw
            u = u + U(m,n)*cos(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT
            v = v + V(m,n)*sin(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
            w = w + W(m,n)*sin(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT

            q(i) = U(m,n);
            q(i+1) = V(m,n);
            q(2*ceil(M/2)*ceil(N/2)+j) = W(m,n);    % modes in w-direction are stored at the end of the vector

            dq(i) = dU(m,n);
            dq(i+1) = dV(m,n);
            dq(2*ceil(M/2)*ceil(N/2)+j) = dW(m,n);

            i = i + 2;
            j = j + 1;
        else
            u = u + U(m,n)*cos(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT
            v = v + V(m,n)*sin(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
                
            q(k) = U(m,n);
            q(k+1) = V(m,n);
                
            dq(k) = dU(m,n);
            dq(k+1) = dV(m,n);
                
            k = k + 2;
        end     
    end
end            
disp('Modal extension for movable edges calculated.')

u = u
v = v
w = w

% geometric imperfections
for m = 1:2:M0
    for n = 1:2:N0
        w0 = w0 + W0(m,n)*sin(m*pi*x/a)*sin(n*pi*y/b);  % CORRECT
    end
end

q = q
shell.q = q;
shell.dq = dq;

%% Additional terms for BCs

% simply supported shell with movable edges

uh = 0; vh = 0; u_sk = 0; v_sk = 0; u_ij = 0; v_ij = 0;

if (strcmp(bc,'simplified')) 
    for m = 1:2:1
        for n = 1:2:1
            for j = 1:2:N0
                for i = 1:2:M0
                    u_ij = u_ij + i/(m+i)*W0(i,j)*sin(j*pi*y/b)*sin((m+i)*pi*x/a);
                    v_ij = v_ij + j/(n+j)*W0(i,j)*sin(i*pi*x/a)*sin((n+j)*pi*y/b);
                end
            end
            uh = uh - 1/8*(m*pi/a*W(m,n)^2 - m*pi/a*W(m,n)^2*cos(2*n*pi*y/b))*sin(2*m*pi*x/a) - m*pi/a*W(m,n)*sin(n*pi*y/b)*u_ij; % FIXED
            vh = vh - 1/8*(n*pi/b*W(m,n)^2 - n*pi/b*W(m,n)^2*cos(2*m*pi*x/a))*sin(2*n*pi*y/b) - n*pi/b*W(m,n)*sin(m*pi*x/a)*v_ij; % FIXED  
        end
    end
    disp('Additional SIMPLIFIED expansion terms for Boundary Conditions added.')
else
    for m = 1:2:Mw
        for k = 1:2:Nw
            for s = 1:2:Mw
                u_sk = u_sk + s/(m+s)*W(s,k)*sin(k*pi*y/b)*sin((m+s)*pi*x/a);       % CORRECT
            end
        end
        for j = 1:2:N0
            for i = 1:2:M0
                u_ij = u_ij + i/(m+i)*W0(i,j)*sin(j*pi*y/b)*sin((m+i)*pi*x/a);      % CORRECT
            end
        end
        for n = 1:2:Nw
            for k = 1:2:Nw
                for s = 1:2:Mw
                    v_sk = v_sk + k/(n+k)*W(s,k)*sin(s*pi*x/a)*sin((n+k)*pi*y/b);       % CORRECT
                end
            end
            for j = 1:2:N0
                for i = 1:2:M0
                    v_ij = v_ij + j/(n+j)*W0(i,j)*sin(i*pi*x/a)*sin((n+j)*pi*y/b);      % CORRECT
                end
            end
            uh = uh - m*pi/a*(1/2*W(m,n)*sin(n*pi*y/b)*u_sk + W(m,n)*sin(n*pi*y/b)*u_ij); % CORRECT
            vh = vh - n*pi/b*(1/2*W(m,n)*sin(m*pi*x/a)*v_sk + W(m,n)*sin(m*pi*x/a)*v_ij);	% CORRECT  
            v_sk = 0; v_ij = 0;
        end
        u_sk = 0; u_ij = 0;
    end
    disp('Additional FULL expansion terms for Boundary Conditions added.')
end

uh = uh
vh = vh

% additional terms are added to the modal extension
u = u + uh; % FIXED
v = v + vh; % FIXED
%%
% plot_path = strcat(path,'Plots\Functions\');
% 
% u_plot = subs(u,[a b W0(1,1)],[par.a par.b par.A11]);
% v_plot = subs(v,[a b W0(1,1)],[par.a par.b par.A11]);
% w_plot = subs(w,[a b W0(1,1)],[par.a par.b par.A11]);
% 
% u_fname = sprintf('u_plot_%i.m',n_dof);
% v_fname = sprintf('v_plot_%i.m',n_dof);
% w_fname = sprintf('w_plot_%i.m',n_dof);
% 
% if(exist(u_fname,'file') == 2)
%     % do nothing
% else
%     matlabFunction(u_plot,'file',u_fname,'Vars',{q x y},'File',strcat(plot_path,u_fname));
%     matlabFunction(v_plot,'file',v_fname,'Vars',{q x y},'File',strcat(plot_path,v_fname));
%     matlabFunction(w_plot,'file',w_fname,'Vars',{q x y},'File',strcat(plot_path,w_fname));
% end

%% Strain-Displacement (Donnell)

% x and y have to be replaced by psi and theta for differentiation.
u = subs(u, [x y], [Rx*psi Ry*theta]);      % CORRECT
v = subs(v, [x y], [Rx*psi Ry*theta]);      % CORRECT
w = subs(w, [x y], [Rx*psi Ry*theta]);      % CORRECT
w0 = subs(w0, [x y], [Rx*psi Ry*theta]);    % CORRECT

% Donnell (much simpler)
eps_x0 = diff(u,psi)/Rx + w/Rx + 1/2*(diff(w,psi)/Rx)^2 + diff(w,psi)/Rx*diff(w0,psi)/Rx;           % CORRECT x2

eps_y0 = diff(v,theta)/Ry + w/Ry + 1/2*(diff(w,theta)/Ry)^2 + diff(w,theta)/Ry*diff(w0,theta)/Ry;   % CORRECT x2

gam_xy0 = diff(u,theta)/Ry + diff(v,psi)/Rx + diff(w,psi)/Rx*diff(w,theta)/Ry ...
          + diff(w,psi)/Rx*diff(w0,theta)/Ry + diff(w0,psi)/Rx*diff(w,theta)/Ry;                    % CORRECT x2

k_x = -diff(w,psi,2)/(Rx^2);            % CORRECT x2

k_y = -diff(w,theta,2)/(Ry^2);          % CORRECT x2

k_xy = -2*diff(w,psi,theta)/(Rx*Ry);    % CORRECT x2

% reverse substitution
u = subs(u, [psi theta], [x/Rx y/Ry]);      % CORRECT
v = subs(v, [psi theta], [x/Rx y/Ry]);      % CORRECT
w = subs(w, [psi theta], [x/Rx y/Ry]);      % CORRECT
w0 = subs(w0, [psi theta], [x/Rx y/Ry]);    % CORRECT

% psi and theta are replaced by x and y again.
eps_x0  = subs(eps_x0,  [psi theta], [x/Rx y/Ry]);  % CORRECT
eps_y0  = subs(eps_y0,  [psi theta], [x/Rx y/Ry]);  % CORRECT
gam_xy0 = subs(gam_xy0, [psi theta], [x/Rx y/Ry]);  % CORRECT

k_x  = subs(k_x,  [psi theta], [x/Rx y/Ry]);    % CORRECT
k_y  = subs(k_y,  [psi theta], [x/Rx y/Ry]);    % CORRECT
k_xy = subs(k_xy, [psi theta], [x/Rx y/Ry]);    % CORRECT

%% Force and Moment Resultants

N_x = E*h/(1-nu^2)*(eps_x0 + nu*eps_y0);
N_y = E*h/(1-nu^2)*(eps_y0 + nu*eps_x0);
M_x = E*h^3/(12*(1-nu^2))*(k_x + nu*k_y);
M_y = E*h^3/(12*(1-nu^2))*(k_y + nu*k_x);

disp('All force and moment resultants must be zero:')
N_xa = simplify(subs(N_x,x,a),'Steps',10)
N_x0 = simplify(subs(N_x,x,0),'Steps',10)
N_yb = simplify(subs(N_y,y,b),'Steps',10)
N_y0 = simplify(subs(N_y,y,0),'Steps',10)

M_xa = subs(M_x,x,a)
M_x0 = subs(M_x,x,0)
M_yb = subs(M_y,y,b)
M_y0 = subs(M_y,y,0)

%% Strain Energy
U_s = sym('U_s');

% For Donnell's theory the term U_cou that couples bending and stretching
% engery has to be neglected.

Us_fname = sprintf(['Us_don_%ddof_' bc '.mat'],n_dof); % filename contains number of dofs 

if(exist(Us_fname,'file') == 2)
    load(Us_fname)       % load Us from file to save calculations
else
    % stretching energy
    t1 = expand(eps_x0^2);                  % expanding saves calculation time
    tic
    t1 = int(t1,x,0,a);
    t1 = int(t1,y,0,b);
    t_x0 = toc
%     t1 = simplify(t1,'Seconds',10);         % CORRECT x2

    t2 = expand(eps_y0^2);
    tic
    t2 = int(t2,y,0,b);         
    t2 = int(t2,x,0,a);
    t_y0 = toc
%     t2 = simplify(t2,'Seconds',10);         % CORRECT x2

    t3 = expand(eps_x0*eps_y0);
    tic
    t3 = int(t3,x,0,a);
    t3 = 2*nu*int(t3,y,0,b);
    t_xy0 = toc
%     t3 = simplify(t3,'Seconds',10);         % CORRECT x2

    t4 = expand(gam_xy0^2);
    tic
    t4 = int(t4,x,0,a);
    t4 = (1 - nu)/2*int(t4,y,0,b);
    t_gam = toc
%     t4 = simplify(t4,'Seconds',10);         % CORRECT x2

    U_str = 1/2*E*h/(1 - nu^2)*(t1 + t2 + t3 + t4); % -C-O-R-R-E-C-T- % FIXED

    % bending energy
    tic
    t5 = expand(k_x^2);
    t5 = int(t5,x,0,a);
    t5 = int(t5,y,0,b);
%     t5 = simplify(t5,'Seconds',10);     % CORRECT x2

    t6 = expand(k_y^2);
    t6 = int(t6,x,0,a);
    t6 = int(t6,y,0,b);
%     t6 = simplify(t6,'Seconds',10);     % CORRECT x2

    t7 = expand(k_x*k_y);
    t7 = int(t7,x,0,a);
    t7 = 2*nu*int(t7,y,0,b);
%     t7 = simplify(t7,'Seconds',10);     % FIXED x2

    t8 = expand(k_xy^2);
    t8 = int(t8,x,0,a);
    t8 = (1 - nu)/2*int(t8,y,0,b);
%     t8 = simplify(t8,'Seconds',10);     % CORRECT x2
    t_ben = toc
 
    U_ben = E*h^3/(24*(1 - nu^2))*(t5 + t6 + t7 + t8);  % FIXED 

        
    U_s = U_str + U_ben;
    U_s = simplify(U_s,'Steps',15);
% The simplify-commands have been dropped to avoid errors.
end

% Save Us to file
save(strcat(Us_path,Us_fname),'U_s');

%% Kinetic Energy (for plotting)
% The kinetic energy is not really needed for the equations of motion,
% because the derivation in the Lagrangian approach will simplify it to a
% large degree. It is, however, nice for plotting in post-processing.
% diff(diff(T,dq),t) = rho * h * a * b / 4 * ddq

% T_fname = sprintf(['T_don_%ddof_10ab%d%d_1000h%d_R%d%d' bc '.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry); % filename contains number of dofs 
% T_name = sprintf(['T_don_%ddof.m'],n_dof);
% if(exist(T_fname,'file') == 2)
%     T_fun = str2func(T_fname(1:end-2));
%     load(T_name)
%     disp('T_fun loaded from file.')
% else
%     tic
%     T_s = expand(du^2 + dv^2 + dw^2);
%     T_s = int(T_s,x,0,a);
%     T_s = 1/2*rho*h*int(T_s,y,0,b)
%     save(strcat(T_path,T_name),'T_S');
%     T = subs(T_s,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11});
%     t_T = toc
%     
%     T_fun = matlabFunction(T,'file','T_fun','Vars',{dq},'File',strcat(T_path,T_fname));
%     disp('T_fun created with matlabFunction().')
% end
% 
% shell.T = T_fun;


%% Damping
% modal damping for simply supported edges
% The complexity of the formulation is again due to the order of dofs.
% The vector of damping coefficients looks like:
% c = [c11 c13 c31 c33 c35 c53 c55]

c = par.zeta*par.rho*par.h*par.a*par.b.*par.eig/2; % numerical values of the damping coefficients
F = 0;

for m = 1:2:M
    for n = 1:2:N
        if m == 1 && n == 1
            F = F + c(1)*(dW(m,n)^2 + dU(m,n)^2 + dV(m,n)^2); 
        elseif (m == 1 && n == 3) || (m == 3 && n == 1)             
            F = F + c(2)*(dU(m,n)^2 + dV(m,n)^2);             
        else
            F = F + c(3)*(dU(m,n)^2 + dV(m,n)^2);
        end
    end
end

F = simplify(a*b/8*F)

%% External Loading

% normal force
f_N = sym('f_N');       % concentrated normal force
ome_N = sym('ome_f');   % frequency in [rad/s]

W_f = f_N*cos(ome_N*t)*subs(w,[x y],[a/2 b/2]);

W_ext = W_f;

%% Lagrange Equations

dUs = sym('dUs%d',[n_dof 1]);
Q   = sym('Q%d',[n_dof 1]);
J = sym('J_%d%d',[2*n_dof,2*n_dof]);

% EOM: ddq = 4/(rho*h*a*b)*[diff(W_ext,q) - diff(F,dq) - diff(Us,q)]
for i = 1:n_dof
    dUs(i) = 4/(rho*h*a*b)*diff(U_s,q(i));                         % CORRECT
    Q(i)   = 4/(rho*h*a*b)*(-diff(F,dq(i)) + diff(W_ext,q(i)));   % CORRECT
end

dUs = dUs
Q = Q

X = [q;dq];       % left-hand side of the EOM  
EOM = [dq;Q-dUs]; % right-hand side of the EOM
EOM = simplify(EOM,'Steps',100);

for i = 1:2*n_dof
    J(:,i) = diff(EOM,X(i));
end
J = simplify(J,'Steps',100);
% J = [diff(rhs,u11) diff(rhs,v11) ... diff(rhs,w33) diff(rhs,du11) ... diff(rhs,dw33)]
shell.J = J;

%% Non-Dimensionalization

xi = sym('xi%d', [n_dof 1]);
dxi = sym('dxi%d',[n_dof 1]);   % non-dimensionalized vectors of dofs and their derivatives
T_11 = sym('T_11');             % period of first eigenfrequency
tau = sym('tau');               % non-dimensionalized time
J_nd = sym('J_nd_%d%d',[2*n_dof,2*n_dof]);

dUs_nd = T_11^2/h*subs(dUs,q,h*xi)
Q_nd = subs(Q,dq,h/T_11*dxi);
Q_nd = T_11^2/h*subs(Q_nd,t,tau*T_11)

Xi = [xi;dxi];
EOM_nd = [dxi;Q_nd-dUs_nd];
EOM_nd = simplify(EOM_nd,'Steps',100);

% Jacobian
for i = 1:2*n_dof
    J_nd(:,i) = diff(EOM_nd,Xi(i));
end
J_nd = simplify(J_nd,'Steps',100);
shell.J_nd = J_nd;


%% Matlab functions

% substitute all constants by numerical values
dUs = subs(dUs,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11});
% Q = subs(Q, c, c_val);
Q = subs(Q,{rho a b h},{par.rho par.a par.b par.h});
% EOM = subs(rhs,c, c_val);
EOM = subs(EOM,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11});
% J = subs(J, c, c_val);
J = subs(J,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11});

% The filenames are quite long as they contain most parameters that might
% be subject to change.
Q_fname = sprintf(['Q_don_%ddof_10ab%d%d_1000h%d_R%d%d_1000z11%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry,1000*par.zeta(1));
dUs_fname = sprintf(['dUs_don_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);
EOM_fname = sprintf(['EOM_don_%ddof_10ab%d%d_1000h%d_R%d%d_1000z11%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry,1000*par.zeta(1));
J_fname = sprintf(['J_don_%ddof_10ab%d%d_1000h%d_R%d%d_1000z11%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry,1000*par.zeta(1));

if (exist(Q_fname) == 2)
    Q_fun = str2func(Q_fname(1:end-2));
    disp('Q_fun loaded from file.')
else
    tic
    Q_fun = matlabFunction(Q,'file','Q_fun','Vars',{X [f_N ome_N] t},'File',strcat(Q_path,Q_fname));
    t_Q=toc
    disp('Q_fun created with matlabFunction().')
end

if (exist(dUs_fname) == 2)
    dUs_fun = str2func(dUs_fname(1:end-2));
    disp('dUs_fun loaded from file.')
else
    tic
    dUs_fun = matlabFunction(dUs,'file','dU_fun','Vars',{X t},'File',strcat(Us_path,dUs_fname));
    t_dUs=toc
    disp('dUs_fun created with matlabFunction().')
end

if (exist(EOM_fname) == 2)
    EOM_fun = str2func(EOM_fname(1:end-2));
    disp('J_fun loaded from file.')
else
    tic
    EOM_fun = matlabFunction(EOM,'file','EOM_fun','Vars',{X [f_N ome_N] t},'File',strcat(EOM_path,EOM_fname));
    t_EOM=toc
    disp('J_fun created with matlabFunction().')
end

if (exist(J_fname) == 2)
    J_fun = str2func(J_fname(1:end-2));
    disp('J_fun loaded from file.')
else
    tic
    J_fun = matlabFunction(J,'file','J_fun','Vars',{X t},'File',strcat(J_path,J_fname));
    t_J=toc
    disp('J_fun created with matlabFunction().')
end

%% Matlab functions for non-dimensionalized EOM

% substitute all constants by numerical values
dUs_nd = subs(dUs_nd,{rho a b h E Rx Ry nu W0(1,1) T_11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T});
% Q_nd = subs(Q_nd, c, c_val);
Q_nd = subs(Q_nd,{rho a b h T_11},{par.rho par.a par.b par.h par.T});
% J_nd = subs(J_nd, c, c_val);
J_nd = subs(J_nd,{rho a b h E Rx Ry nu W0(1,1) T_11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T});
EOM_nd = subs(EOM_nd,{rho a b h E Rx Ry nu W0(1,1) T_11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T});
% EOM_nd = subs(EOM_nd, c, c_val);

% names for function files
Q_nd_fname = sprintf(['Q_don_nd_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);
dUs_nd_fname = sprintf(['dUs_don_nd_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);
J_nd_fname = sprintf(['J_don_nd_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);
EOM_nd_fname = sprintf(['EOM_don_nd_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);

if (exist(Q_nd_fname) == 2)
    Q_nd_fun = str2func(Q_nd_fname(1:end-2));
    disp('Q_nd_fun loaded from file.')
else
    tic
    Q_nd_fun = matlabFunction(Q_nd,'file','Q_nd_fun','Vars',{Xi [f_N ome_N] tau},'File',strcat(Q_path,Q_nd_fname));
    t_Q=toc
    disp('Q_nd_fun created with matlabFunction().')
end
 
if (exist(dUs_nd_fname) == 2)
    dUs_nd_fun = str2func(dUs_nd_fname(1:end-2));
    disp('dUs_nd_fun loaded from file.')
else
    tic
    dUs_nd_fun = matlabFunction(dUs_nd,'file','dU_nd_fun','Vars',{Xi tau},'File',strcat(Us_path,dUs_nd_fname));
    t_dUs=toc
    disp('dUs_nd_fun created with matlabFunction().')
end

if (exist(J_nd_fname) == 2)
    J_nd_fun = str2func(J_nd_fname(1:end-2));
    disp('J_nd_fun loaded from file.')
else
    tic
    J_nd_fun = matlabFunction(J_nd,'file','J_nd_fun','Vars',{Xi tau},'File',strcat(J_path,J_nd_fname));
    t_J=toc
    disp('J_nd_fun created with matlabFunction().')
end

if (exist(EOM_nd_fname) == 2)
    EOM_nd_fun = str2func(EOM_nd_fname(1:end-2));
    disp('EOM_nd_fun loaded from file.')
else
    tic
    EOM_nd_fun = matlabFunction(EOM_nd,'file','EOM_nd_fun','Vars',{Xi, [f_N ome_N], tau},'File',strcat(EOM_path,EOM_nd_fname));
    t_EOM=toc
    disp('EOM_fun created with matlabFunction().')
end


%% Return function handles

shell.Q = Q_fun;
shell.dU = dUs_fun;
shell.EOM = EOM_fun; 
shell.Jf = J_fun;
shell.dU_nd = dUs_nd_fun;
shell.Q_nd = Q_nd_fun;
shell.Jf_nd = J_nd_fun;
shell.EOM_nd = EOM_nd_fun;

%%
disp('Symbolic Computations Completed!')
end
