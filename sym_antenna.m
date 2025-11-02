function [shell] = sym_antenna(par,n_dof,path)
% This function performs the symbolic computations for the antenna model
% using the ansatz proposed by Jin.

shell = struct;


% path = 'E:\Dokumente\Studium\Master\4th_Semester\Projektarbeit_Strukturdynamik\Matlab\'; % path to project folder
Us_path = strcat(path,'Elastic_Forces\Jin\'); % path to Elastic_Forces folder
Q_path = strcat(path,'External_Forces\Jin\');
T_path = strcat(path,'Kinetic_Energy\Jin\');
J_path = strcat(path,'Jacobians\Jin\');
EOM_path = strcat(path,'EOM\Jin\');

%% Symbolic Parameters

Rx = sym('Rx'); Ry = sym('Ry');                    % radii of curvature
t = sym('t');                                      % time
x = sym('x');  y = sym('y');                       % coordinates along the edges
psi = sym('psi'); theta = sym('theta');            % angular coordinates
a = sym('a'); b = sym('b'); h = sym('h');          % curvilinear dimensions of the shell
E = sym('E'); nu = sym('nu'); rho = sym('rho');    % material properties (Young's modulus, Poisson's ratio and density)

% x = Rx * psi and y = Ry * theta
%% DoFs
% n_dof = 9;
if n_dof == 22
    M = 5; N = 5;
    Mw = 3; Nw = 3;
    M0 = 1; N0 = 1;
else
    M = 10; N = 10;
    Mw = 10; Nw = 10;
    M0 = 1; N0 = 1;
end

% Matlab cannot automatically generate symbolic variables with indices
% including a "0" so they have to assigend manually.
syms u00 v00 w00
syms u10 v10 w10 u01 v01 w01
syms u20 v20 w20 u02 v02 w02
syms u30 v30 w30 u03 v03 w03
syms u40 v40 w40 u04 v04 w04
syms u50 v50 w50 u05 v05 w05
syms u60 v60 w60 u06 v06 w06
syms u70 v70 w70 u07 v07 w07
syms u80 v80 w80 u08 v08 w08
syms u90 v90 w90 u09 v09 w09
syms u100 v100 w100 u010 v010 w010

syms du00 dv00 dw00
syms du10 dv10 dw10 du01 dv01 dw01
syms du20 dv20 dw20 du02 dv02 dw02
syms du30 dv30 dw30 du03 dv03 dw03
syms du40 dv40 dw40 du04 dv04 dw04
syms du50 dv50 dw50 du05 dv05 dw05
syms du60 dv60 dw60 du06 dv06 dw06
syms du70 dv70 dw70 du07 dv07 dw07
syms du80 dv80 dw80 du08 dv08 dw08
syms du90 dv90 dw90 du09 dv09 dw09
syms du100 dv100 dw100 du010 dv010 dw010

syms ddu00 ddv00 ddw00
syms ddu10 ddv10 ddw10 ddu01 ddv01 ddw01
syms ddu20 ddv20 ddw20 ddu02 ddv02 ddw02
syms ddu30 ddv30 ddw30 ddu03 ddv03 ddw03
syms ddu40 ddv40 ddw40 ddu04 ddv04 ddw04
syms ddu50 ddv50 ddw50 ddu05 ddv05 ddw05
syms ddu60 ddv60 ddw60 ddu06 ddv06 ddw06
syms ddu70 ddv70 ddw70 ddu07 ddv07 ddw07
syms ddu80 ddv80 ddw80 ddu08 ddv08 ddw08
syms ddu90 ddv90 ddw90 ddu09 ddv09 ddw09
syms ddu100 ddv100 ddw100 ddu010 ddv010 ddw010

% Modal expansion
U = sym('u%d%d', [10 10]);
V = sym('v%d%d', [10 10]);
W = sym('w%d%d', [10 10]);      % modal amplitudes (time-dependent)
% W0 = sym('A%d%d', [1 1]);     % geometric imperfections

U = [u00 u01 u02 u03 u04 u05 u06 u07 u08 u09 u010;
     u10 U(1,1:10);
     u20 U(2,1:10);
     u30 U(3,1:10);
     u40 U(4,1:10);
     u50 U(5,1:10);
     u60 U(6,1:10);
     u70 U(7,1:10);
     u80 U(8,1:10);
     u90 U(9,1:10);
     u100 U(10,1:10)
];
 
V = [v00 v01 v02 v03 v04 v05 v06 v07 v08 v09 v010;
     v10 V(1,1:10);
     v20 V(2,1:10);
     v30 V(3,1:10);
     v40 V(4,1:10);
     v50 V(5,1:10);
     v60 V(6,1:10);
     v70 V(7,1:10);
     v80 V(8,1:10);
     v90 V(9,1:10);
     v100 V(10,1:10)
];
 
W = [w00 w01 w02 w03 w04 w05 w06 w07 w08 w09 w010;
     w10 W(1,1:10);
     w20 W(2,1:10);
     w30 W(3,1:10);
     w40 W(4,1:10);
     w50 W(5,1:10);
     w60 W(6,1:10);
     w70 W(7,1:10);
     w80 W(8,1:10);
     w90 W(9,1:10);
     w100 W(10,1:10)
];
 
% U = [ u11 u12 ... ]           % the entries with even indices are not
%     [ u21 u22 ... ]             needed
%     [ u31 u32 u33 ]

dU = sym('du%d%d', [10 10]);
dV = sym('dv%d%d', [10 10]);
dW = sym('dw%d%d', [10 10]);    % modal velocities

dU = [du00 du01 du02 du03 du04 du05 du06 du07 du08 du09 du010;
      du10 dU(1,1:10);
      du20 dU(2,1:10);
      du30 dU(3,1:10);
      du40 dU(4,1:10);
      du50 dU(5,1:10);
      du60 dU(6,1:10);
      du70 dU(7,1:10);
      du80 dU(8,1:10);
      du90 dU(9,1:10);
      du100 dU(10,1:10)
];
 
dV = [dv00 dv01 dv02 dv03 dv04 dv05 dv06 dv07 dv08 dv09 dv010;
      dv10 dV(1,1:10);
      dv20 dV(2,1:10);
      dv30 dV(3,1:10);
      dv40 dV(4,1:10);
      dv50 dV(5,1:10);
      dv60 dV(6,1:10);
      dv70 dV(7,1:10);
      dv80 dV(8,1:10);
      dv90 dV(9,1:10);
      dv100 dV(10,1:10)
];
 
dW = [dw00 dw01 dw02 dw03 dw04 dw05 dw06 dw07 dw08 dw09 dw010;
      dw10 dW(1,1:10);
      dw20 dW(2,1:10);
      dw30 dW(3,1:10);
      dw40 dW(4,1:10);
      dw50 dW(5,1:10);
      dw60 dW(6,1:10);
      dw70 dW(7,1:10);
      dw80 dW(8,1:10);
      dw90 dW(9,1:10);
      dw100 dW(10,1:10)
];
  

ddU = sym('ddu%d%d', [10 10]);
ddV = sym('ddv%d%d', [10 10]);
ddW = sym('ddw%d%d', [10 10]);    % modal velocities

ddU = [ddu00 ddu01 ddu02 ddu03 ddu04 ddu05 ddu06 ddu07 ddu08 ddu09 ddu010;
       ddu10 ddU(1,1:10);
       ddu20 ddU(2,1:10);
       ddu30 ddU(3,1:10);
       ddu40 ddU(4,1:10);
       ddu50 ddU(5,1:10);
       ddu60 ddU(6,1:10);
       ddu70 ddU(7,1:10);
       ddu80 ddU(8,1:10);
       ddu90 ddU(9,1:10);
       ddu100 ddU(10,1:10)
];
 
ddV = [ddv00 ddv01 ddv02 ddv03 ddv04 ddv05 ddv06 ddv07 ddv08 ddv09 ddv010;
       ddv10 ddV(1,1:10);
       ddv20 ddV(2,1:10);
       ddv30 ddV(3,1:10);
       ddv40 ddV(4,1:10);
       ddv50 ddV(5,1:10);
       ddv60 ddV(6,1:10);
       ddv70 ddV(7,1:10);
       ddv80 ddV(8,1:10);
       ddv90 ddV(9,1:10);
       ddv100 ddV(10,1:10)
];
 
ddW = [ddw00 ddw01 ddw02 ddw03 ddw04 ddw05 ddw06 ddw07 ddw08 ddw09 ddw010;
       ddw10 ddW(1,1:10);
       ddw20 ddW(2,1:10);
       ddw30 ddW(3,1:10);
       ddw40 ddW(4,1:10);
       ddw50 ddW(5,1:10);
       ddw60 ddW(6,1:10);
       ddw70 ddW(7,1:10);
       ddw80 ddW(8,1:10);
       ddw90 ddW(9,1:10);
       ddw100 ddW(10,1:10)
];

u = 0; v = 0; w = 0; w0 = 0;


% modal extension for simply supported edges
% The formulation might seem unneccessarily complicated. This was done
% to bring the dofs into a nice order for plotting. The vector of modal
% amplitudes will look like: q = [u11 v11 u13 v13 ... u55 v55 w11 ...
% w33]'

n_dof = 2*(M+1)*(N+1) + (Mw+1)*(Nw+1) + 2*2*(N+1) + 2*2*(M+1) + 2*(Nw+1) + 2*(Mw+1);
n_mod = 2*(M+1)*(N+1) + (Mw+1)*(Nw+1);

q = sym('q%d', [n_dof 1]);
dq = sym('dq%d',[n_dof 1]);     % vectors of dofs and their derivatives
ddq = sym('ddq%d', [n_dof 1]);

xi = sym('xi%d', [n_dof 1]);
dxi = sym('dxi%d',[n_dof 1]);   % non-dimensionalized vectors of dofs and their derivatives
    
k = 2*(Mw+1)*(Nw+1) + 1;   % number of modes in w-direction plus 1
j = 1; i = 1;
n_uv = 2*(M+1)*(N+1);

 for m = 0:M
     for n = 0:N
         if n <= Nw && m <= Mw
            u = u + U(m+1,n+1)*cos(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
            v = v + V(m+1,n+1)*cos(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
            w = w + W(m+1,n+1)*cos(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT

            q(i) = U(m+1,n+1);
            q(i+1) = V(m+1,n+1);
            q(n_uv+j) = W(m+1,n+1);    % modes in w-direction are stored at the end of the vector

            dq(i) = dU(m+1,n+1);
            dq(i+1) = dV(m+1,n+1);
            dq(n_uv+j) = dW(m+1,n+1);
            
            ddq(i) = ddU(m+1,n+1);
            ddq(i+1) = ddV(m+1,n+1);
            ddq(n_uv+j) = ddW(m+1,n+1);

            i = i + 2;
            j = j + 1;
        else
            u = u + U(m+1,n+1)*cos(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
            v = v + V(m+1,n+1)*cos(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
                
            q(k) = U(m+1,n+1);
            q(k+1) = V(m+1,n+1);
                
            dq(k) = dU(m+1,n+1);
            dq(k+1) = dV(m+1,n+1);
            
            ddq(k) = ddU(m+1,n+1);
            ddq(k+1) = ddV(m+1,n+1);
                
            k = k + 2;
        end     
    end
end            
disp('Modal extension for movable edges calculated.')

u = u
v = v
w = w

% geometric imperfections
% for m = 1:2:M0
%     for n = 1:2:N0
%         w0 = w0 + W0(m,n)*sin(m*pi*x/a)*sin(n*pi*y/b);  % CORRECT
%     end
% end

q = q
shell.q = q;
shell.dq = dq;

%% Jin's additional functions
syms a_u10 a_u20 a_v10 a_v20 a_w10 a_w20
syms b_u10 b_u20 b_v10 b_v20 b_w10 b_w20
syms da_u10 da_u20 da_v10 da_v20 da_w10 da_w20
syms db_u10 db_u20 db_v10 db_v20 db_w10 db_w20
syms dda_u10 dda_u20 dda_v10 dda_v20 dda_w10 dda_w20
syms ddb_u10 ddb_u20 ddb_v10 ddb_v20 ddb_w10 ddb_w20

a_u = sym('a_u%d%d', [2 N]); a_v = sym('a_v%d%d', [2 N]); a_w = sym('a_w%d%d', [2 Nw]);
b_u = sym('b_u%d%d', [2 M]); b_v = sym('b_v%d%d', [2 M]); b_w = sym('b_w%d%d', [2 Mw]);
a_u = [[a_u10; a_u20] a_u]; a_v = [[a_v10; a_v20] a_v]; a_w = [[a_w10; a_w20] a_w];
b_u = [[b_u10; b_u20] b_u]; b_v = [[b_v10; b_v20] b_v]; b_w = [[b_w10; b_w20] b_w];

da_u = sym('da_u%d%d', [2 N]); da_v = sym('da_v%d%d', [2 N]); da_w = sym('da_w%d%d', [2 Nw]);
db_u = sym('db_u%d%d', [2 M]); db_v = sym('db_v%d%d', [2 M]); db_w = sym('db_w%d%d', [2 Mw]);
da_u = [[da_u10; da_u20] da_u]; da_v = [[da_v10; da_v20] da_v]; da_w = [[da_w10; da_w20] da_w];
db_u = [[db_u10; db_u20] db_u]; db_v = [[db_v10; db_v20] db_v]; db_w = [[db_w10; db_w20] db_w];

dda_u = sym('dda_u%d%d', [2 N]); dda_v = sym('dda_v%d%d', [2 N]); dda_w = sym('dda_w%d%d', [2 Nw]);
ddb_u = sym('ddb_u%d%d', [2 M]); ddb_v = sym('ddb_v%d%d', [2 M]); ddb_w = sym('ddb_w%d%d', [2 Mw]);
dda_u = [[dda_u10; dda_u20] dda_u]; dda_v = [[dda_v10; dda_v20] dda_v]; dda_w = [[dda_w10; dda_w20] dda_w];
ddb_u = [[ddb_u10; ddb_u20] ddb_u]; ddb_v = [[ddb_v10; ddb_v20] ddb_v]; ddb_w = [[ddb_w10; ddb_w20] ddb_w];

xi_x_u = 0; xi_x_v = 0; xi_x_w = 0; xi_y_u = 0; xi_y_v = 0; xi_y_w = 0;
j = 1; k = 1;

n_a = 2*2*(N+1) + 2*2*(M+1);

for n = 0:N
    xi_x_u = xi_x_u + a_u(:,n+1)*cos(n*pi*y/b);
    q(n_mod+j:n_mod+j+1) = a_u(:,n+1);
    dq(n_mod+j:n_mod+j+1) = da_u(:,n+1);
    ddq(n_mod+j:n_mod+j+1) = dda_u(:,n+1);
    j = j + 2;
    xi_x_v = xi_x_v + a_v(:,n+1)*cos(n*pi*y/b);
    q(n_mod+j:n_mod+j+1) = a_v(:,n+1);
    dq(n_mod+j:n_mod+j+1) = da_v(:,n+1);
    ddq(n_mod+j:n_mod+j+1) = dda_v(:,n+1);
    j = j + 2;
    if n <= Nw
        xi_x_w = xi_x_w + a_w(:,n+1)*cos(n*pi*y/b);
        q(n_mod+n_a+k:n_mod+n_a+k+1) = a_w(:,n+1);
        dq(n_mod+n_a+k:n_mod+n_a+k+1) = da_w(:,n+1);
        ddq(n_mod+n_a+k:n_mod+n_a+k+1) = dda_w(:,n+1);
        k = k + 2;
    end
end


for m = 0:M
    xi_y_u = xi_y_u + b_u(:,m+1)*cos(m*pi*x/a);
    q(n_mod+j:n_mod+j+1) = b_u(:,m+1);
    dq(n_mod+j:n_mod+j+1) = db_u(:,m+1);
    ddq(n_mod+j:n_mod+j+1) = ddb_u(:,m+1);
    j = j + 2;
    xi_y_v = xi_y_v + b_v(:,m+1)*cos(m*pi*x/a);
    q(n_mod+j:n_mod+j+1) = b_v(:,m+1);
    dq(n_mod+j:n_mod+j+1) = db_v(:,m+1);
    ddq(n_mod+j:n_mod+j+1) = ddb_v(:,m+1);
    j = j + 2;
    if m <= Mw
        xi_y_w = xi_y_w + b_w(:,m+1)*cos(m*pi*x/a);
        q(n_mod+n_a+k:n_mod+n_a+k+1) = b_w(:,m+1);
        dq(n_mod+n_a+k:n_mod+n_a+k+1) = db_w(:,m+1);
        ddq(n_mod+n_a+k:n_mod+n_a+k+1) = ddb_w(:,m+1);
        k = k + 2;
    end
end

P1x = x*(x/a-1)^2;
P2x = x^2*(x/a-1)/a;
P1y = y*(y/b-1)^2;
P2y = y^2*(y/b-1)/b;

Pxu = [P1x P2x]*xi_x_u; Pxv = [P1x P2x]*xi_x_v; Pxw = [P1x P2x]*xi_x_w;
Pyu = [P1y P2y]*xi_y_u; Pyv = [P1y P2y]*xi_y_v; Pyw = [P1y P2y]*xi_y_w;

u = u + Pxu + Pyu
v = v + Pxv + Pyv
w = w + Pxw + Pyw

shell.q = q;
shell.dq = dq;


%% Strain-Displacement (Donnell)

% x and y have to be replaced by psi and theta for differentiation.
u = subs(u, [x y], [Rx*psi Ry*theta]);      % CORRECT
v = subs(v, [x y], [Rx*psi Ry*theta]);      % CORRECT
w = subs(w, [x y], [Rx*psi Ry*theta]);      % CORRECT
w0 = subs(w0, [x y], [Rx*psi Ry*theta]);    % CORRECT


% Donnell (much simpler than Novozhilov's)
tic
eps_x0 = diff(u,psi)/Rx + w/Rx + 1/2*(diff(w,psi)/Rx)^2 + diff(w,psi)/Rx*diff(w0,psi)/Rx;           % CORRECT x2

eps_y0 = diff(v,theta)/Ry + w/Ry + 1/2*(diff(w,theta)/Ry)^2 + diff(w,theta)/Ry*diff(w0,theta)/Ry;   % CORRECT x2

gam_xy0 = diff(u,theta)/Ry + diff(v,psi)/Rx + diff(w,psi)/Rx*diff(w,theta)/Ry ...
          + diff(w,psi)/Rx*diff(w0,theta)/Ry + diff(w0,psi)/Rx*diff(w,theta)/Ry;                    % CORRECT x2

k_x = -diff(w,psi,2)/(Rx^2);            % CORRECT x2

k_y = -diff(w,theta,2)/(Ry^2);          % CORRECT x2

k_xy = -2*diff(w,psi,theta)/(Rx*Ry);    % CORRECT x2
t_straindisp = toc

% % Jin (linear strain-displacement)
% 
% tic
% eps_x0 = diff(u,psi)/Rx + w/Rx;             % CORRECT x2
% 
% eps_y0 = diff(v,theta)/Ry + w/Ry;           % CORRECT x2
% 
% gam_xy0 = diff(u,theta)/Ry + diff(v,psi)/Rx;                    % CORRECT x2
% 
% k_x = -diff(w,psi,2)/(Rx^2);            % CORRECT x2
% 
% k_y = -diff(w,theta,2)/(Ry^2);          % CORRECT x2
% 
% k_xy = -2*diff(w,psi,theta)/(Rx*Ry);    % CORRECT x2
% t_straindisp = toc

% reverse substitution
u = subs(u, [psi theta], [x/Rx y/Ry]);      % CORRECT
v = subs(v, [psi theta], [x/Rx y/Ry]);      % CORRECT
w = subs(w, [psi theta], [x/Rx y/Ry]);      % CORRECT
w0 = subs(w0, [psi theta], [x/Rx y/Ry]);    % CORRECT

eps_x0  = subs(eps_x0,  [psi theta], [x/Rx y/Ry]);  % CORRECT
eps_y0  = subs(eps_y0,  [psi theta], [x/Rx y/Ry]);  % CORRECT
gam_xy0 = subs(gam_xy0, [psi theta], [x/Rx y/Ry]);  % CORRECT

k_x  = subs(k_x,  [psi theta], [x/Rx y/Ry]);    % CORRECT
k_y  = subs(k_y,  [psi theta], [x/Rx y/Ry]);    % CORRECT
k_xy = subs(k_xy, [psi theta], [x/Rx y/Ry]); 


%% Force and Moment Resultants

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
U_s = sym('Us');

% For Donnell's theory the term U_cou that couples bending and stretching
% engery has to be neglected.

Us_fname = sprintf(['Us_lin_%ddof.mat'],n_dof); % filename contains number of dofs 

if(exist(Us_fname,'file') == 2)
    load(Us_fname)       % load Us from file to save calculations
    disp('Us loaded from file.')
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
%     Us = simplify(Us,'Steps',15);
% The simplify-commands have been dropped to avoid errors.
end

% Save Us to file
save(strcat(Us_path,Us_fname),'U_s');

%% Kinetic Energy (for plotting)
% The kinetic energy is not really needed for the equations of motion,
% because the derivation in the Lagrangian approach will simplify it to a
% large degree. It is, however, nice for plotting in post-processing.
% diff(diff(T,dq),t) = rho * h * a * b / 4 * ddq

T_fname = sprintf(['T_jin_%ddof.mat'],n_dof);

du = subs(u,q,dq); dv = subs(v,q,dq); dw = subs(w,q,dq);

if(exist(T_fname,'file') == 2)
    load(T_fname)
    disp('T_s loaded from file.')
else
    tic
    T_s = expand(du^2 + dv^2 + dw^2);
    T_s = int(T_s,x,0,a);
    T_s = 1/2*rho*h*int(T_s,y,0,b);
    t_T = toc
    disp('T_s calculated.')
    save(strcat(T_path,T_fname),'T_s')
end

%% Damping #TODO: series expansion

syms c1 c2 c3
F = 0; F_a1 = 0; F_a2 = 0; F_a3 = 0; F_b1 = 0; F_b2 = 0; F_b3 = 0;

for m = 1:M+1
    for n = 1:N+1
        if m <= Mw+1 && n <= Nw+1
            if m == Mw+1
                for i = 1:2
                    F_a1 = F_a1 + da_u(i,n)^2 + da_v(i,n)^2 + da_w(i,n)^2;
                end
            end
            F = F + c1*(dU(m,n)^2 + dV(m,n)^2 + dW(m,n)^2);
        elseif n > Nw+1 && m <= 4 && n <= 4 && not(n == m && n == 4)
            if m == 4
                for i = 1:2
                    F_a2 = F_a2 + da_u(i,n)^2 + da_v(i,n)^2;
                end
            end
            F = F + c2*(dU(m,n)^2 + dV(m,n)^2);
        else
            if m == M+1 &&  n >= 4
                for i = 1:2
                    F_a3 = F_a3 + da_u(i,n)^2 + da_v(i,n)^2;
                end
            end
            F = F + c3*(dU(m,n)^2 + dV(m,n)^2);
        end
    end
    if m <= Mw+1
           for i = 1:2
               F_b1 = F_b1 + db_u(i,m)^2 + db_v(i,m)^2 + db_w(i,m)^2;
           end
    elseif m > Mw+1 && m <= 3
           for i = 1:2
               F_b2 = F_b2 + db_u(i,m)^2 + db_v(i,m)^2;
           end
    else
           for i = 1:2
               F_b3 = F_b3 + db_u(i,m)^2 + db_v(i,m)^2;
           end
    end
    
end

F = F + c1*(F_a1 + F_b1)+ c2*(F_a2 + F_b2) + c3*(F_a3 + F_b3);


% F = simplify(a*b/8*F)

%% External Loading

% concentrated normal force
f_N = sym('f_N');       % concentrated normal force
ome_N = sym('ome_N');   % frequency in [rad/s]

W_f = f_N*cos(ome_N*t)*subs(w,[x y],[a/2 b/2]);

W_ext = W_f;

%% Artificical Springs

syms k_x0_u k_x0_v k_x0_w K_x0_w
syms k_x1_u k_x1_v k_x1_w K_x1_w
syms k_y0_u k_y0_v k_y0_w K_y0_w
syms k_y1_u k_y1_v k_y1_w K_y1_w

% k = num2cell(par.k);
% [k_x0_u,k_x0_v,k_x0_w,K_x0_w,...
%     k_x1_u,k_x1_v,k_x1_w,K_x1_w,...
%     k_y0_u,k_y0_v,k_y0_w,K_y0_w,...
%     k_y1_u,k_y1_v,k_y1_w,K_y1_w] = deal(k{:});

k_sp = [k_x0_u,k_x0_v,k_x0_w,K_x0_w,...
        k_x1_u,k_x1_v,k_x1_w,K_x1_w,...
        k_y0_u,k_y0_v,k_y0_w,K_y0_w,...
        k_y1_u,k_y1_v,k_y1_w,K_y1_w];
    
U_sp = 1/2*int(subs(k_x0_u*u^2 + k_x0_v*v^2 + k_x0_w*w^2 + K_x0_w*diff(w,x),x,0)... 
        + subs(k_x1_u*u^2 + k_x1_v*v^2 + k_x1_w*w^2 + K_x1_w*diff(w,x),x,a),y,0,b) ...
        + 1/2*int(subs(k_y0_u*u^2 + k_y0_v*v^2 + k_y0_w*w^2 + K_y0_w*diff(w,y),y,0) ...
        + subs(k_y1_u*u^2 + k_y1_v*v^2 + k_y1_w*w^2 + K_y1_w*diff(w,y),y,b),x,0,a);

%% Derivatives of T with respect to the generalized coordinates

dT = sym('dT%d',[n_dof,1]);
dTdq = sym('dTdq%d',[n_dof,1]);

for i=1:n_dof
    dT(i) = diff(T_s,dq(i));
    dTdq(i) = diff(T_s,q(i));
end

dTdt = subs(dT,dq,ddq);
% There are no quadratic terms in the generalized velocities in dT so
% taking the derivative w.r.t. time corresponds to replacing the modal
% velocities by the accelerations.

%% Lagrange Equations

EOM = sym('EOM%d',[n_dof 1]);
L = sym('L%d',[n_dof 1]);
J = sym('J_nd_%d%d',[2*n_dof,2*n_dof]);
Q = sym('Q%d',[n_dof 1]);
dU_s = sym('dU%d',[n_dof 1]);

tic
for i = 1:n_dof
    Q(i) = -diff(F,dq(i)) + diff(W_ext,q(i));
    dU_s(i) = diff(U_s+U_sp,q(i));
    L(i) = dTdt(i) - dTdq(i) + dU_s(i) == Q(i);
end
toc

%% Solve L for the accelerations
% The output from the "solve"-command cannot be assigned to dynamically
% generated variables but have to assigend manually.

% tic
% [d2q1, d2q2, d2q3, d2q4, d2q5, d2q6, d2q7, d2q8, d2q9,...
%  d2q10, d2q11, d2q12, d2q13, d2q14, d2q15, d2q16, d2q17, d2q18,...
%  d2q19, d2q20, d2q21, d2q22, d2q23, d2q24, d2q25, d2q26, d2q27,...
%  d2q28, d2q29, d2q30, d2q31, d2q32, d2q33, d2q34, d2q35, d2q36] = solve(L,ddq);
% toc

% M=N=3
  [d2q1, d2q2, d2q3, d2q4, d2q5, d2q6, d2q7, d2q8, d2q9, ...
   d2q10, d2q11, d2q12, d2q13, d2q14, d2q15, d2q16, d2q17, d2q18, ...
   d2q19, d2q20, d2q21, d2q22, d2q23, d2q24, d2q25, d2q26, d2q27, ...
   d2q28, d2q29, d2q30, d2q31, d2q32, d2q33, d2q34, d2q35, d2q36, ...
   d2q37, d2q38, d2q39, d2q40, d2q41, d2q42, d2q43, d2q44, d2q45, ...
   d2q46, d2q47, d2q48, d2q49, d2q50, d2q51, d2q52, d2q53, d2q54, ...
   d2q55, d2q56, d2q57, d2q58, d2q59, d2q60, d2q61, d2q62, d2q63, ...
   d2q64, d2q65, d2q66, d2q67, d2q68, d2q69, d2q70, d2q71, d2q72, ...
   d2q73, d2q74, d2q75, d2q76, d2q77, d2q78, d2q79, d2q80, d2q81, ...
   d2q82, d2q83, d2q84, d2q85, d2q86, d2q87, d2q88, d2q89, d2q90,...
   d2q91, d2q92, d2q93, d2q94, d2q95, d2q96] = solve(L,ddq);

%% Save the righ-hand side as variable.

% d2q = [d2q1; d2q2; d2q3; d2q4; d2q5; d2q6; d2q7; d2q8; d2q9;...
%       d2q10; d2q11; d2q12; d2q13; d2q14; d2q15; d2q16; d2q17; d2q18;...
%       d2q19; d2q20; d2q21; d2q22; d2q23; d2q24; d2q25; d2q26; d2q27;...
%       d2q28; d2q29; d2q30; d2q31; d2q32; d2q33; d2q34; d2q35; d2q36];


d2q = [d2q1; d2q2; d2q3; d2q4; d2q5; d2q6; d2q7; d2q8; d2q9; ...
   d2q10; d2q11; d2q12; d2q13; d2q14; d2q15; d2q16; d2q17; d2q18; ...
   d2q19; d2q20; d2q21; d2q22; d2q23; d2q24; d2q25; d2q26; d2q27; ...
   d2q28; d2q29; d2q30; d2q31; d2q32; d2q33; d2q34; d2q35; d2q36; ...
   d2q37; d2q38; d2q39; d2q40; d2q41; d2q42; d2q43; d2q44; d2q45; ...
   d2q46; d2q47; d2q48; d2q49; d2q50; d2q51; d2q52; d2q53; d2q54; ...
   d2q55; d2q56; d2q57; d2q58; d2q59; d2q60; d2q61; d2q62; d2q63; ...
   d2q64; d2q65; d2q66; d2q67; d2q68; d2q69; d2q70; d2q71; d2q72; ...
   d2q73; d2q74; d2q75; d2q76; d2q77; d2q78; d2q79; d2q80; d2q81; ...
   d2q82; d2q83; d2q84; d2q85; d2q86; d2q87; d2q88; d2q89; d2q90; ...
   d2q91; d2q92; d2q93; d2q94; d2q95; d2q96];

save(strcat(EOM_path,sprintf('EOM_lin_%d.mat',n_dof)),'d2q')

%% Compute the Jacobian
X = [q;dq];
rhs = [dq;d2q];

% Jacobian
for i = 1:2*n_dof
    tic
    J(:,i) = diff(rhs,X(i));
    toc
end

%% Non-Dimensionalization

T11 = sym('T11'); % period of first eigenfrequency
tau = sym('tau');   % non-dimensionalized time

xi = sym('xi%d',[n_dof,1]);
dxi = sym('dxi%d',[n_dof,1]);
EOM_nd = sym('EOM_nd%d',[2*n_dof,1]);
J_nd = sym('J_nd_%d%d',[2*n_dof,2*n_dof]);
Xi = [xi;dxi];


EOM_nd(1:n_dof) = dxi;
EOM_nd(n_dof+1:end) = T11^2/h*subs(d2q,[q; dq; t],[h*xi; h/T11*dxi; T11*tau]);
EOM_nd = subs(EOM_nd,{h T11},{par.h par.T11});
disp('EOM_nd calculated.')

tic
for i = 1:2*n_dof
%     J(:,i) = (diff(EOM,X(i)));
    J_nd(:,i) = diff(EOM_nd,Xi(i));
end
toc
disp('Jacobians calculated.')


% Jacobian
disp('Non-dimensionalization complete!')

shell.J_nd = J_nd;

%% Substitute symbolic variables by numerical values
% substitute all constants by numerical values

% dU_s = subs(dU_s,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11});
% Q = subs(Q, {c1 c2 c3}, {par.c});
% Q = subs(Q,{rho a b h W0(1,1)},{par.rho par.a par.b par.h par.A11});
EOM = subs(rhs,{c1 c2 c3}, {par.c});
EOM = simplify(subs(EOM,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11}),'Seconds',600);
J = subs(J,{c1 c2 c3}, {par.c});
J = simplify(subs(J,{rho a b h E Rx Ry nu W0(1,1)},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11}),'Seconds',1200);
EOM_nd = subs(EOM_nd,{c1 c2 c3}, {par.c});
%%
EOM_nd = simplify(subs(EOM_nd,{rho a b h E Rx Ry nu W0(1,1) T11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T11}),'Seconds',600);
J_nd = subs(J_nd,{c1 c2 c3}, {par.c});
J_nd = simplify(subs(J_nd,{rho a b h E Rx Ry nu W0(1,1) T11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T11}),'Seconds',1200);
disp('Substitutions complete!')

%% Matlab functions

% The filenames are quite long as they contain most parameters that might
% be subject to change.

% dUs_fname = sprintf(['dUs_jin_%ddof_10a%d_R%d.m'],n_dof,10*par.a,par.Rx);
% Q_fname = sprintf(['Q_jin_%ddof_10a%d_R%d%d.m'],n_dof,10*par.a,par.Rx);
EOM_fname = sprintf(['EOM_lin_%ddof_10a%d_R%ds.m'],n_dof,10*par.a,par.Rx);
J_fname = sprintf(['J_lin_%ddof_10a%d_R%ds.m'],n_dof,10*par.a,par.Rx);

% if (exist(Q_fname) == 2)
%     Q_fun = str2func(Q_fname(1:end-2));
%     disp('Q_fun loaded from file.')
% else
%     tic
%     Q_fun = matlabFunction(Q,'file','Q_fun','Vars',{dq [f_N ome_N] t},'File',strcat(Q_path,Q_fname));
%     t_Q=toc
%     disp('Q_fun created with matlabFunction().')
% end
% 
% if (exist(dUs_fname) == 2)
%     dUs_fun = str2func(dUs_fname(1:end-2));
%     disp('dUs_fun loaded from file.')
% else
%     tic
%     dUs_fun = matlabFunction(dU_s,'file','dU_fun','Vars',{q k_sp t},'File',strcat(Us_path,dUs_fname));
%     t_dUs=toc
%     disp('dUs_fun created with matlabFunction().')
% end
%%
if (exist(EOM_fname) == 2)
    EOM_fun = str2func(EOM_fname(1:end-2));
    disp('EOM_fun loaded from file.')
else
    tic
    EOM_fun = matlabFunction(EOM,'file','EOM_fun','Vars',{X k_sp [f_N ome_N] t},'File',strcat(EOM_path,EOM_fname));
    t_EOM=toc
    disp('EOM_fun created with matlabFunction().')
end
%%
if (exist(J_fname) == 2)
    J_fun = str2func(J_fname(1:end-2));
    disp('J_fun loaded from file.')
else
    tic
    J_fun = matlabFunction(J,'file','J_fun','Vars',{X k_sp t},'File',strcat(J_path,J_fname));
    t_J=toc
    disp('J_fun created with matlabFunction().')
end

%%

EOM_fname = sprintf(['EOM_jin_nd_%ddof_10a%d_1000h%d_R%ds.m'],n_dof,10*par.a,par.Rx);
if (exist(EOM_fname) == 2)
    EOM_fun = str2func(EOM_fname(1:end-2));
    disp('EOM_fun loaded from file.')
else
    tic
    EOM_fun = matlabFunction(EOM,'file','EOM_fun','Vars',{X k_sp [f_N ome_N] t},'File',strcat(EOM_path,EOM_fname));
    t_EOM=toc
    disp('EOM_fun created with matlabFunction().')
end

%% Matlab functions for non-dimensionalized EOM

% % substitute all constants by numerical values
% J_nd = subs(J_nd, c, c_val);
% J_nd = subs(J_nd,{rho a b h E Rx Ry nu W0(1,1) T_11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T});
% EOM_nd = subs(EOM_nd,{rho a b h E Rx Ry nu W0(1,1) T_11},{par.rho par.a par.b par.h par.E par.Rx par.Ry par.nu par.A11 par.T});
% EOM_nd = subs(EOM_nd, c, c_val);

% names for function files
J_nd_fname = sprintf(['J_jin_nd_%ddof_10a%d%d_R%d.m'],n_dof,10*par.a,par.Rx);
EOM_nd_fname = sprintf(['EOM_jin_nd_%ddof_10a%d_R%d.m'],n_dof,10*par.a,par.Rx);

%%
if (exist(J_nd_fname) == 2)
    J_nd_fun = str2func(J_nd_fname(1:end-2));
    disp('J_nd_fun loaded from file.')
else
    tic
    J_nd_fun = matlabFunction(J_nd,'file','J_nd_fun','Vars',{Xi k_sp tau},'File',strcat(J_path,J_nd_fname));
    t_J=toc
    disp('J_nd_fun created with matlabFunction().')
end

%%

if (exist(EOM_nd_fname) == 2)
    EOM_nd_fun = str2func(EOM_nd_fname(1:end-2));
    disp('EOM_nd_fun loaded from file.')
else
    tic
    EOM_nd_fun = matlabFunction(EOM_nd,'file','EOM_nd_fun','Vars',{Xi, [f_N ome_N], k_sp, tau},'File',strcat(EOM_path,EOM_nd_fname));
    t_EOM=toc
    disp('EOM_fun created with matlabFunction().')
end

% %% Return function handles
% 
% shell.Q = Q_fun;
% shell.dU = dUs_fun;
% shell.Jf = J_fun;
% shell.dU_nd = dUs_nd_fun;
% shell.Q_nd = Q_nd_fun;
% shell.Jf_nd = J_nd_fun;
% shell.EOM_nd = EOM_nd_fun;

%%
disp('Symbolic Computations Completed!')
(h*rho*((a^3*b*da_u10^2)/105 + (a^3*b*da_u11^2)/210 + (a^3*b*da_v10^2)/105 + (a^3*b*da_v11^2)/210 + (a^3*b*da_w10^2)/105 + (a^3*b*da_w11^2)/210 + (a^3*b*da_u20^2)/105 + (a^3*b*da_u21^2)/210 + (a^3*b*da_v20^2)/105 + (a^3*b*da_v21^2)/210 + (a^3*b*da_w20^2)/105 + (a^3*b*da_w21^2)/210 + (a*b^3*db_u10^2)/105 + (a*b^3*db_u11^2)/210 + (a*b^3*db_v10^2)/105 + (a*b^3*db_v11^2)/210 + (a*b^3*db_w10^2)/105 + (a*b^3*db_w11^2)/210 + (a*b^3*db_u20^2)/105 + (a*b^3*db_u21^2)/210 + (a*b^3*db_v20^2)/105 + (a*b^3*db_v21^2)/210 + (a*b^3*db_w20^2)/105 + (a*b^3*db_w21^2)/210 + a*b*du00^2 + (a*b*du01^2)/2 + (a*b*du10^2)/2 + (a*b*du11^2)/4 + a*b*dv00^2 + (a*b*dv01^2)/2 + (a*b*dv10^2)/2 + (a*b*dv11^2)/4 + a*b*dw00^2 + (a*b*dw01^2)/2 + (a*b*dw10^2)/2 + (a*b*dw11^2)/4 - (a^3*b*da_u10*da_u20)/70 - (a^3*b*da_u11*da_u21)/140 - (a^3*b*da_v10*da_v20)/70 - (a^3*b*da_v11*da_v21)/140 - (a^3*b*da_w10*da_w20)/70 - (a^3*b*da_w11*da_w21)/140 - (a*b^3*db_u10*db_u20)/70 - (a*b^3*db_u11*db_u21)/140 - (a*b^3*db_v10*db_v20)/70 - (a*b^3*db_v11*db_v21)/140 - (a*b^3*db_w10*db_w20)/70 - (a*b^3*db_w11*db_w21)/140 + (a^2*b*da_u10*du00)/6 + (a^2*b*da_u11*du01)/12 - (a^2*b*da_u20*du00)/6 - (a^2*b*da_u21*du01)/12 + (a*b^2*db_u10*du00)/6 + (a^2*b*da_v10*dv00)/6 + (a^2*b*da_v11*dv01)/12 - (a*b^2*db_u20*du00)/6 + (a*b^2*db_u11*du10)/12 - (a^2*b*da_v20*dv00)/6 - (a^2*b*da_v21*dv01)/12 - (a*b^2*db_u21*du10)/12 + (a*b^2*db_v10*dv00)/6 + (a^2*b*da_w10*dw00)/6 + (a^2*b*da_w11*dw01)/12 - (a*b^2*db_v20*dv00)/6 + (a*b^2*db_v11*dv10)/12 - (a^2*b*da_w20*dw00)/6 - (a^2*b*da_w21*dw01)/12 - (a*b^2*db_v21*dv10)/12 + (a*b^2*db_w10*dw00)/6 - (a*b^2*db_w20*dw00)/6 + (a*b^2*db_w11*dw10)/12 - (a*b^2*db_w21*dw10)/12 + (a^2*b^2*da_u10*db_u10)/72 + (a^2*b^2*da_v10*db_v10)/72 + (a^2*b^2*da_w10*db_w10)/72 - (a^2*b^2*da_u10*db_u20)/72 - (a^2*b^2*da_u20*db_u10)/72 - (a^2*b^2*da_v10*db_v20)/72 - (a^2*b^2*da_v20*db_v10)/72 - (a^2*b^2*da_w10*db_w20)/72 - (a^2*b^2*da_w20*db_w10)/72 + (a^2*b^2*da_u20*db_u20)/72 + (a^2*b^2*da_v20*db_v20)/72 + (a^2*b^2*da_w20*db_w20)/72 - (a^2*b^2*da_u10*db_u11)/(6*pi^2) - (a^2*b^2*da_u11*db_u10)/(6*pi^2) + (2*a^2*b^2*da_u10*db_u11)/pi^4 + (2*a^2*b^2*da_u11*db_u10)/pi^4 - (a^2*b^2*da_v10*db_v11)/(6*pi^2) - (a^2*b^2*da_v11*db_v10)/(6*pi^2) + (2*a^2*b^2*da_u11*db_u11)/pi^4 + (2*a^2*b^2*da_v10*db_v11)/pi^4 + (2*a^2*b^2*da_v11*db_v10)/pi^4 - (a^2*b^2*da_w10*db_w11)/(6*pi^2) - (a^2*b^2*da_w11*db_w10)/(6*pi^2) - (48*a^2*b^2*da_u11*db_u11)/pi^6 + (2*a^2*b^2*da_v11*db_v11)/pi^4 + (2*a^2*b^2*da_w10*db_w11)/pi^4 + (2*a^2*b^2*da_w11*db_w10)/pi^4 + (288*a^2*b^2*da_u11*db_u11)/pi^8 - (48*a^2*b^2*da_v11*db_v11)/pi^6 + (2*a^2*b^2*da_w11*db_w11)/pi^4 + (288*a^2*b^2*da_v11*db_v11)/pi^8 - (48*a^2*b^2*da_w11*db_w11)/pi^6 + (a^2*b^2*da_u10*db_u21)/(6*pi^2) - (a^2*b^2*da_u11*db_u20)/(6*pi^2) - (a^2*b^2*da_u20*db_u11)/(6*pi^2) + (a^2*b^2*da_u21*db_u10)/(6*pi^2) + (288*a^2*b^2*da_w11*db_w11)/pi^8 - (2*a^2*b^2*da_u10*db_u21)/pi^4 + (2*a^2*b^2*da_u11*db_u20)/pi^4 + (a^2*b^2*da_v10*db_v21)/(6*pi^2) - (a^2*b^2*da_v11*db_v20)/(6*pi^2) + (2*a^2*b^2*da_u20*db_u11)/pi^4 - (2*a^2*b^2*da_u21*db_u10)/pi^4 - (a^2*b^2*da_v20*db_v11)/(6*pi^2) + (a^2*b^2*da_v21*db_v10)/(6*pi^2) + (2*a^2*b^2*da_u11*db_u21)/pi^4 + (2*a^2*b^2*da_u21*db_u11)/pi^4 - (2*a^2*b^2*da_v10*db_v21)/pi^4 + (2*a^2*b^2*da_v11*db_v20)/pi^4 + (a^2*b^2*da_w10*db_w21)/(6*pi^2) - (a^2*b^2*da_w11*db_w20)/(6*pi^2) + (2*a^2*b^2*da_v20*db_v11)/pi^4 - (2*a^2*b^2*da_v21*db_v10)/pi^4 - (a^2*b^2*da_w20*db_w11)/(6*pi^2) + (a^2*b^2*da_w21*db_w10)/(6*pi^2) - (48*a^2*b^2*da_u11*db_u21)/pi^6 + (2*a^2*b^2*da_v11*db_v21)/pi^4 - (48*a^2*b^2*da_u21*db_u11)/pi^6 + (2*a^2*b^2*da_v21*db_v11)/pi^4 - (2*a^2*b^2*da_w10*db_w21)/pi^4 + (2*a^2*b^2*da_w11*db_w20)/pi^4 + (2*a^2*b^2*da_w20*db_w11)/pi^4 - (2*a^2*b^2*da_w21*db_w10)/pi^4 + (288*a^2*b^2*da_u11*db_u21)/pi^8 - (48*a^2*b^2*da_v11*db_v21)/pi^6 + (2*a^2*b^2*da_w11*db_w21)/pi^4 + (288*a^2*b^2*da_u21*db_u11)/pi^8 - (48*a^2*b^2*da_v21*db_v11)/pi^6 + (2*a^2*b^2*da_w21*db_w11)/pi^4 + (288*a^2*b^2*da_v11*db_v21)/pi^8 - (48*a^2*b^2*da_w11*db_w21)/pi^6 + (288*a^2*b^2*da_v21*db_v11)/pi^8 - (48*a^2*b^2*da_w21*db_w11)/pi^6 + (a^2*b^2*da_u20*db_u21)/(6*pi^2) + (a^2*b^2*da_u21*db_u20)/(6*pi^2) + (288*a^2*b^2*da_w11*db_w21)/pi^8 + (288*a^2*b^2*da_w21*db_w11)/pi^8 - (2*a^2*b^2*da_u20*db_u21)/pi^4 - (2*a^2*b^2*da_u21*db_u20)/pi^4 + (a^2*b^2*da_v20*db_v21)/(6*pi^2) + (a^2*b^2*da_v21*db_v20)/(6*pi^2) + (2*a^2*b^2*da_u21*db_u21)/pi^4 - (2*a^2*b^2*da_v20*db_v21)/pi^4 - (2*a^2*b^2*da_v21*db_v20)/pi^4 + (a^2*b^2*da_w20*db_w21)/(6*pi^2) + (a^2*b^2*da_w21*db_w20)/(6*pi^2) - (48*a^2*b^2*da_u21*db_u21)/pi^6 + (2*a^2*b^2*da_v21*db_v21)/pi^4 - (2*a^2*b^2*da_w20*db_w21)/pi^4 - (2*a^2*b^2*da_w21*db_w20)/pi^4 + (288*a^2*b^2*da_u21*db_u21)/pi^8 - (48*a^2*b^2*da_v21*db_v21)/pi^6 + (2*a^2*b^2*da_w21*db_w21)/pi^4 + (288*a^2*b^2*da_v21*db_v21)/pi^8 - (48*a^2*b^2*da_w21*db_w21)/pi^6 + (288*a^2*b^2*da_w21*db_w21)/pi^8 - (2*a^2*b*da_u10*du10)/pi^2 + (24*a^2*b*da_u10*du10)/pi^4 - (a^2*b*da_u11*du11)/pi^2 + (12*a^2*b*da_u11*du11)/pi^4 - (2*a^2*b*da_u20*du10)/pi^2 + (24*a^2*b*da_u20*du10)/pi^4 - (a^2*b*da_u21*du11)/pi^2 + (12*a^2*b*da_u21*du11)/pi^4 - (2*a*b^2*db_u10*du01)/pi^2 + (24*a*b^2*db_u10*du01)/pi^4 - (2*a*b^2*db_u20*du01)/pi^2 - (2*a^2*b*da_v10*dv10)/pi^2 - (a*b^2*db_u11*du11)/pi^2 + (24*a*b^2*db_u20*du01)/pi^4 + (24*a^2*b*da_v10*dv10)/pi^4 - (a^2*b*da_v11*dv11)/pi^2 + (12*a*b^2*db_u11*du11)/pi^4 + (12*a^2*b*da_v11*dv11)/pi^4 - (2*a^2*b*da_v20*dv10)/pi^2 - (a*b^2*db_u21*du11)/pi^2 + (24*a^2*b*da_v20*dv10)/pi^4 - (a^2*b*da_v21*dv11)/pi^2 + (12*a*b^2*db_u21*du11)/pi^4 + (12*a^2*b*da_v21*dv11)/pi^4 - (2*a*b^2*db_v10*dv01)/pi^2 + (24*a*b^2*db_v10*dv01)/pi^4 - (2*a*b^2*db_v20*dv01)/pi^2 - (2*a^2*b*da_w10*dw10)/pi^2 - (a*b^2*db_v11*dv11)/pi^2 + (24*a*b^2*db_v20*dv01)/pi^4 + (24*a^2*b*da_w10*dw10)/pi^4 - (a^2*b*da_w11*dw11)/pi^2 + (12*a*b^2*db_v11*dv11)/pi^4 + (12*a^2*b*da_w11*dw11)/pi^4 - (2*a^2*b*da_w20*dw10)/pi^2 - (a*b^2*db_v21*dv11)/pi^2 + (24*a^2*b*da_w20*dw10)/pi^4 - (a^2*b*da_w21*dw11)/pi^2 + (12*a*b^2*db_v21*dv11)/pi^4 + (12*a^2*b*da_w21*dw11)/pi^4 - (2*a*b^2*db_w10*dw01)/pi^2 + (24*a*b^2*db_w10*dw01)/pi^4 - (2*a*b^2*db_w20*dw01)/pi^2 - (a*b^2*db_w11*dw11)/pi^2 + (24*a*b^2*db_w20*dw01)/pi^4 + (12*a*b^2*db_w11*dw11)/pi^4 - (a*b^2*db_w21*dw11)/pi^2 + (12*a*b^2*db_w21*dw11)/pi^4))/2
