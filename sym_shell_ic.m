function [shell] = sym_shell_ic(par,n_dof,path)
% This function performs the symbolic computations for the shell with inertial
% coupling.
%% Paths

U_path = strcat(path,'Elastic_Forces\Donnell\'); % path to Elastic_Forces folder
Q_path = strcat(path,'External_Forces\Donnell\');
T_path = strcat(path,'Kinetic_Energy\Donnell\');
J_path = strcat(path,'Jacobians\Donnell\');
EOM_path = strcat(path,'EOM\Donnell\');

%% Degrees of Freedom

if n_dof == 22
    M = 5; N = 5;
    Mw = 3; Nw = 3;
    M0 = 1; N0 = 1;
else
    M = 3; N = 3;
    Mw = 1; Nw = 1;
    M0 = 1; N0 = 1;
end

n_uv = 2*ceil(M/2)*ceil(N/2);
n_w = ceil(Mw/2)*ceil(Nw/2); 

qt = sym('q%d',[n_dof 1]);
dqt = sym('dtq%d',[n_dof 1]);
q = sym('q%d',[n_dof 1]);
dq = sym('dq%d',[n_dof 1]);
ddq = sym('ddq%d',[n_dof 1]);

%% Symbolic Parameters

Rx = sym('Rx'); Ry = sym('Ry');                    % radii of curvature
t = sym('t');                                      % time
x = sym('x');  y = sym('y');                       % coordinates along the edges
psi = sym('psi'); theta = sym('theta');            % angular coordinates
a = sym('a'); b = sym('b'); h = sym('h');          % curvilinear dimensions of the shell
E = sym('E'); nu = sym('nu'); rho = sym('rho');    % material properties (Young's modulus, Poisson's ratio and density)

% x = Rx * psi and y = Ry * theta
%% Modal amplitudes
% Both time-dependent and a time-independent versions of the amplitudes and
% time derivatives are needed, as the time-dependent amplitudes can be used
% for time-derivatives, but symbolic functions cannot be derived w.r.t.
% these time-dependent ampltiudes as they are symbolic functions
% themselves. For derivatives w.r.t. the amplitudes, they are substituted
% by independent variables.

syms ut11(t) ut13(t) ut31(t) ut33(t) ut15(t) ut51(t) ut35(t) ut53(t) ut55(t) 
syms vt11(t) vt13(t) vt31(t) vt33(t) vt15(t) vt51(t) vt35(t) vt53(t) vt55(t)
syms wt11(t) wt13(t) wt31(t) wt33(t)
syms A11    % CORRECT

syms dtu11(t) dtu13(t) dtu31(t) dtu33(t) dtu15(t) dtu51(t) dtu35(t) dtu53(t) dtu55(t) 
syms dtv11(t) dtv13(t) dtv31(t) dtv33(t) dtv15(t) dtv51(t) dtv35(t) dtv53(t) dtv55(t)
syms dtw11(t) dtw13(t) dtw31(t) dtw33(t)    % CORRECT

U = sym('u%d%d',[M N]);   % time-independent velocities
V = sym('v%d%d',[M N]);
W = sym('w%d%d',[Mw Nw]);

dU = sym('du%d%d',[M N]);   % time-independent velocities
dV = sym('dv%d%d',[M N]);
dW = sym('dw%d%d',[Mw Nw]);

ddU = sym('ddu%d%d',[M N]); % time-independent accelerations
ddV = sym('ddv%d%d',[M N]);
ddW = sym('ddw%d%d',[Mw Nw]);

Ut = {ut11 0 ut13 0 ut15;
         0 0    0 0    0;
      ut31 0 ut33 0 ut35;
         0 0    0 0    0;
      ut51 0 ut53 0 ut55};  % CORRECT
Vt = {vt11 0 vt13 0 vt15;
         0 0   0  0    0;
      vt31 0 vt33  0 vt35;
         0 0   0  0    0;
      vt51 0 vt53 0 vt55};  % CORRECT
Wt = {wt11 0 wt13;
         0 0    0;
      wt31 0 wt33};         % CORRECT
W0 = {A11};

dUt = {dtu11 0 dtu13 0 dtu15;
           0 0     0 0     0;
       dtu31 0 dtu33 0 dtu35;
           0 0     0 0     0;
       dtu51 0 dtu53 0 dtu55};  % CORRECT
dVt = {dtv11 0 dtv13 0 dtv15;
           0 0     0 0     0;
       dtv31 0 dtv33 0 dtv35;
           0 0     0 0     0;
       dtv51 0 dtv53 0 dtv55};  % CORRECT
dWt = {dtw11 0 dtw13;
           0 0     0;
       dtw31 0 dtw33};      	% CORRECT

u = 0; v = 0; w = 0; w0 = 0;
i = 1; j = 1; k = 2*n_w + 1;

%%
for m = 1:2:M
     for n = 1:2:N
         if n <= Nw && m <= Mw
            u = u + Ut{m,n}*cos(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT
            v = v + Vt{m,n}*sin(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
            w = w + Wt{m,n}*sin(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT

            qt(i) = Ut{m,n};
            qt(i+1) = Vt{m,n};
            qt(n_uv+j) = Wt{m,n};    % modes in w-direction are stored at the end of the vector

            dqt(i) = dUt{m,n};
            dqt(i+1) = dVt{m,n};
            dqt(n_uv+j) = dWt{m,n};
            
            q(i) = U(m,n);
            q(i+1) = V(m,n);
            q(n_uv+j) = W(m,n);
            
            dq(i) = dU(m,n);
            dq(i+1) = dV(m,n);
            dq(n_uv+j) = dW(m,n);
            
            ddq(i) = ddU(m,n);
            ddq(i+1) = ddV(m,n);
            ddq(n_uv+j) = ddW(m,n);

            i = i + 2;
            j = j + 1;
        else
            u = u + Ut{m,n}*cos(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT
            v = v + Vt{m,n}*sin(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
                
            qt(k) = Ut{m,n};
            qt(k+1) = Vt{m,n};
                
            dqt(k) = dUt{m,n};
            dqt(k+1) = dVt{m,n};
            
            q(k) = U(m,n);
            q(k+1) = V(m,n);
            
            dq(k) = dU(m,n);
            dq(k+1) = dV(m,n);
                           
            ddq(k) = ddU(m,n);
            ddq(k+1) = ddV(m,n);
            
            k = k + 2;
        end     
    end
end

% % geometric imperfections
% for m = 1:2:M0
%     for n = 1:2:N0
%         w0 = w0 + W0{m,n}*cos(m*pi*x/a)*cos(n*pi*y/b);  % CORRECT
%     end
% end

u = u
v = v
w = w
% w0 = w0

%% Additional terms for Boundary Conditions

uh = 0; vh = 0; u_sk = 0; v_sk = 0; u_ij = 0; v_ij = 0;

for m = 1:2:Mw
    for k = 1:2:Nw
        for s = 1:2:Mw
            u_sk = u_sk + s/(m+s)*Wt{s,k}*sin(k*pi*y/b)*sin((m+s)*pi*x/a);       % CORRECT
        end
    end
    for j = 1:2:N0
        for i = 1:2:M0
            u_ij = u_ij + i/(m+i)*W0{i,j}*sin(j*pi*y/b)*sin((m+i)*pi*x/a);      % CORRECT
        end
    end
    for n = 1:2:Nw
        for k = 1:2:Nw
            for s = 1:2:Mw
                v_sk = v_sk + k/(n+k)*Wt{s,k}*sin(s*pi*x/a)*sin((n+k)*pi*y/b);       % CORRECT
            end
        end
        for j = 1:2:N0
            for i = 1:2:M0
                v_ij = v_ij + j/(n+j)*W0{i,j}*sin(i*pi*x/a)*sin((n+j)*pi*y/b);      % CORRECT
            end
        end
        uh = uh - m*pi/a*(1/2*Wt{m,n}*sin(n*pi*y/b)*u_sk + Wt{m,n}*sin(n*pi*y/b)*u_ij); % CORRECT
        vh = vh - n*pi/b*(1/2*Wt{m,n}*sin(m*pi*x/a)*v_sk + Wt{m,n}*sin(m*pi*x/a)*v_ij);	% CORRECT  
        v_sk = 0; v_ij = 0;
    end
    u_sk = 0; u_ij = 0;
 end
disp('Additional FULL expansion terms for Boundary Conditions added.')

uh = uh
vh = vh

u = u + uh;
v = v + vh;

%% Strain-Displacement (Donnell)

% x and y have to be replaced by psi and theta for differentiation.
u = subs(u, [x y], [Rx*psi Ry*theta]);      % CORRECT
v = subs(v, [x y], [Rx*psi Ry*theta]);      % CORRECT
w = subs(w, [x y], [Rx*psi Ry*theta]);      % CORRECT
w0 = subs(w0, [x y], [Rx*psi Ry*theta]);    % CORRECT

%Donnell (much simpler)
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

% normal forces and bending moments acting on the shell element
N_x = E*h/(1-nu^2)*(eps_x0 + nu*eps_y0);
N_y = E*h/(1-nu^2)*(eps_y0 + nu*eps_x0);
M_x = E*h^3/(12*(1-nu^2))*(k_x + nu*k_y);
M_y = E*h^3/(12*(1-nu^2))*(k_y + nu*k_x);

% force resultants at the edges
N_xa = simplify(subs(N_x,x,a),'Steps',10)
N_x0 = simplify(subs(N_x,x,0),'Steps',10)
N_yb = simplify(subs(N_x,y,b),'Steps',10)
N_y0 = simplify(subs(N_x,y,0),'Steps',10)

% moment resultants at the edges
M_xa = subs(M_x,x,a)
M_x0 = subs(M_x,x,0)
M_yb = subs(M_y,y,b)
M_y0 = subs(M_y,y,0)

%% Strain Energy

U_s = sym('U_s');

% For Donnell's theory the term U_cou that couples bending and stretching
% engery has to be neglected.

Us_fname = sprintf(['Us_ama_%ddof.mat'],n_dof); % filename contains number of dofs 

if(exist(Us_fname,'file') == 2)
    load(Us_fname)       % load Us from file to save calculations
    disp('U_s loaded from file.')
else
    % stretching energy
    t1 = expand(eps_x0^2);                  % expanding saves calculation time
    tic
    t1 = int(t1,x,0,a);
    t1 = int(t1,y,0,b);
    t_x0 = toc

    t2 = expand(eps_y0^2);
    t2 = int(t2,y,0,b);         
    t2 = int(t2,x,0,a);

    t3 = expand(eps_x0*eps_y0);
    tic
    t3 = int(t3,x,0,a);
    t3 = 2*nu*int(t3,y,0,b);
    t_xy0 = toc

    t4 = expand(gam_xy0^2);
    tic
    t4 = int(t4,x,0,a);
    t4 = (1 - nu)/2*int(t4,y,0,b);
    t_gam = toc

    U_str = 1/2*E*h/(1 - nu^2)*(t1 + t2 + t3 + t4);

    % bending energy
    tic
    t5 = expand(k_x^2);
    t5 = int(t5,x,0,a);
    t5 = int(t5,y,0,b);

    t6 = expand(k_y^2);
    t6 = int(t6,x,0,a);
    t6 = int(t6,y,0,b);

    t7 = expand(k_x*k_y);
    t7 = int(t7,x,0,a);
    t7 = 2*nu*int(t7,y,0,b);

    t8 = expand(k_xy^2);
    t8 = int(t8,x,0,a);
    t8 = (1 - nu)/2*int(t8,y,0,b);
    t_ben = toc
 
    U_ben = E*h^3/(24*(1 - nu^2))*(t5 + t6 + t7 + t8);  % FIXED
    
    U_s = U_str + U_ben;
%     U_s = simplify(U_s,'Steps',100);
    disp('U_s calculated.')
    % Save U_s to file
    U_s = subs(U_s,qt,q);
    save(strcat(U_path,Us_fname),'U_s');
    % The simplify-commands have been dropped to avoid errors.
end


%% Kinetic Energy 
% The kinetic energy is not really needed for the equations of motion,
% because the derivation in the Lagrangian approach will simplify it to a
% large degree. It is, however, nice for plotting in post-processing.
% diff(diff(T,dq),t) = rho * h * a * b / 4 * ddq

T_fname = sprintf('T_ama_%ddof.m',n_dof);

if(exist(T_fname,'file') == 2)
    load(T_fname)
    disp('T_s loaded from file.')
else
    tic
    T_s = expand(diff(u,t)^2 + diff(v,t)^2 + diff(w,t)^2);
    T_s = int(T_s,x,0,a);
    T_s = 1/2*rho*h*int(T_s,y,0,b)
    T_s = subs(T_s,[qt;diff(qt,t)],[q;dq]);
    t_T = toc
    save(strcat(T_path,T_fname),'T_s');
    disp('T_s calculated.')
end

%% Damping
% modal damping for simply supported edges
% The complexity of the formulation is again due to the order of dofs.
% The vector of damping coefficients looks like:
% c = [c11 c13 c31 c33 c35 c53 c55]

c = sym('c%d',[3 1]);
% par.zeta*par.rho*par.h*par.a*par.b.*par.eig/2; % numerical values of the damping coefficients
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

F = simplify(a*b/8*F) % CORRECT


%% External Loading

% concentrated normal force
f_N = sym('f_N');       % concentrated normal force
ome_N = sym('ome_N');   % frequency in [rad/s]

W_f = f_N*cos(ome_N*t)*subs(w,[x y],[a/2 b/2]);

W_ext = subs(W_f,qt,q);

%% Lagrange Equations of Motion

dTddq = sym('dtT%d',[n_dof 1]);
dTdq = sym('dqT%d',[n_dof 1]);
dU_s = sym('dUs%d',[n_dof 1]);
Q = sym('Q%d',[n_dof 1]);
L = sym('L%d',[n_dof 1]);


U_s = subs(U_s,qt,q);

for i = 1:n_dof
    dTddq(i) = diff(T_s,dq(i)); % CORRECT
    dTdq(i) = diff(T_s,q(i));   % CORRECT
    dU_s(i) = diff(U_s,q(i));   % CORRECT
    Q(i) = -diff(F,dq(i)) + diff(W_ext,q(i)); % CORRECT
end

dTddq = subs(dTddq,[q;dq],[qt;dqt]);
ddTddqdt = diff(dTddq,t);               % d(dT/dq_i)/dt
ddTddqdt = subs(ddTddqdt,[qt;dqt;diff(qt,t);diff(dqt,t)],[q;dq;dq;ddq]);

%% Lagrange equations

for i = 1:n_dof
    L(i) = ddTddqdt(i) - dTdq(i) + dU_s(i) == Q(i);
end

%% Solve for accelerations

tic
 [d2q1, d2q2, d2q3, ...
  d2q4, d2q5, d2q6, ...
  d2q7, d2q8, d2q9] = solve(L,ddq);
t_sol = toc

d2q = [d2q1; d2q2; d2q3; ...
       d2q4; d2q5; d2q6; ...
       d2q7; d2q8; d2q9];

  
%% Iterative solution

d1q = sym('d1q%d',[n_dof 1]);
d1q(1) = solve(L(1),ddq(1));

for i = 2:n_dof
    tic
    d1q(i) = solve(subs(L(i),ddq(1:i-1),d1q(1:i-1)),ddq(i));
    toc
end

% substitute ddw11 in all other expressions
d1q(1:n_dof-1) = subs(d1q(1:n_dof-1),ddq(n_dof),d1q(n_dof));
d1q = simplify(d1q,'Seconds',300);

for i = 1:n_dof-1
    tic
    d1q(n_dof-i) = subs(d1q(n_dof-i),ddq(n_dof-i+1:n_dof),d1q(n_dof-i+1:n_dof));
    toc
end


%% Equations of Motion, Jacobians and Non-Dimensionalization

T11 = sym('T11'); % period of first eigenfrequency
tau = sym('tau');   % non-dimensionalized time

xi = sym('xi%d',[n_dof 1]);
dxi = sym('dxi%d',[n_dof 1]);

Xi = [xi;dxi];
X = [q;dq];

EOM = [dq;d2q];
EOM_nd = sym('EOM_nd%d',[2*n_dof,1]);
J = sym('J%d',[2*n_dof,2*n_dof]);
J_nd = sym('J_nd%d',[2*n_dof,2*n_dof]);

EOM_nd(1:n_dof) = dxi;
EOM_nd(n_dof+1:end) = T11^2/h*subs(EOM(n_dof+1:end),[q; dq; t],[h*xi; h/T11*dxi; T11*tau]);

tic
for i = 1:2*n_dof
    J(:,i) = (diff(EOM,X(i)));
    J_nd(:,i) = diff(EOM_nd,Xi(i));
end
toc
disp('Jacobians calculated.')

%% Substitutions

% substitute all constants by numerical values
% dU_s = subs(dUs,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11});
% Q = subs(Q, c, par.zeta*par.rho*par.h*par.a*par.b.*par.eig/2);
% Q = subs(Q,{rho a b h},{par.rho par.a par.b par.h});
EOM = subs(EOM, c, par.zeta*par.rho*par.h*par.a*par.b.*par.eig/2);
EOM = subs(EOM,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11});
J = subs(J, c, par.zeta*par.rho*par.h*par.a*par.b.*par.eig/2);
J = subs(J,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11});

EOM_nd = subs(EOM_nd, c, par.zeta*par.rho*par.h*par.a*par.b.*par.eig/2);
EOM_nd = subs(EOM_nd,{rho a b h E Rx Ry nu W0(1,1) T11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T});
J_nd = subs(J_nd, c, par.zeta*par.rho*par.h*par.a*par.b.*par.eig/2);
J_nd = subs(J_nd,{rho a b h E Rx Ry nu W0(1,1) T11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T});

% simplifications
EOM = simplify(EOM,'Seconds',300);
J = simplify(J,'Seconds',480);
EOM_nd = simplify(EOM_nd,'Seconds',300);
J_nd = simplify(J_nd,'Seconds',480);

%% Matlab Functions

J_fname = sprintf(['J_ama_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);
EOM_fname = sprintf(['EOM_ama_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);

J_nd_fname = sprintf(['J_nd_ama_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);
EOM_nd_fname = sprintf(['EOM_nd_ama_%ddof_10ab%d%d_1000h%d_R%d%d.m'],n_dof,10*par.a,10*par.b,1000*par.h,par.Rx,par.Ry);

if (exist(J_fname) == 2)
    J_fun = str2func(J_fname(1:end-2));
    disp('J_fun loaded from file.')
else
    tic
    J_fun = matlabFunction(J,'file','J_fun','Vars',{X [f_N ome_N] t},'File',strcat(J_path,J_fname));
    t_J=toc
    disp('J_nd_fun created with matlabFunction().')
end

if (exist(EOM_fname) == 2)
    EOM_fun = str2func(EOM_fname(1:end-2));
    disp('EOM_nd_fun loaded from file.')
else
    tic
    EOM_fun = matlabFunction(EOM,'file','EOM_fun','Vars',{X, [f_N ome_N], t},'File',strcat(EOM_path,EOM_fname));
    t_EOM=toc
    disp('EOM_fun created with matlabFunction().')
end

if (exist(J_nd_fname) == 2)
    J_nd_fun = str2func(J_nd_filename(1:end-2));
    disp('J_nd_fun loaded from file.')
else
    tic
    J_nd_fun = matlabFunction(J_nd,'file','J_nd_fun','Vars',{Xi [f_N ome_N] tau},'File',strcat(J_path,J_nd_fname));
    t_J=toc
    disp('J_nd_fun created with matlabFunction().')
end

if (exist(EOM_nd_fname) == 2)
    EOM_nd_fun = str2func(EOM_nd_filename(1:end-2));
    disp('EOM_nd_fun loaded from file.')
else
    tic
    EOM_nd_fun = matlabFunction(EOM_nd,'file','EOM_nd_fun','Vars',{Xi, [f_N ome_N], tau},'File',strcat(EOM_path,EOM_nd_fname));
    t_EOM=toc
    disp('EOM_fun created with matlabFunction().')
end

%%

shell.EOM = EOM_fun; 
shell.Jf = J_fun;

shell.Jf_nd = J_nd_fun;
shell.EOM_nd = EOM_nd_fun;
end