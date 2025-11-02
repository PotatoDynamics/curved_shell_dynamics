%% Declaration of Variables
% This script should demonstrate the non-orthogonality of the
% shape-functions used by Amabili in his 2005 paper. Here, the modes are
% expanded for simplicity into only one modal amplitude in x-, y- and z-direction.

syms u11(t) v11(t) w11(t)
syms u(x,y,t) v(x,y,t) w(x,y,t)
syms uh(x,y,t) vh(x,y,t)
syms a b

%% Modal expansion
% The expansions correspond to equations 16 and 25 in the case M = N = Mw =
% Nw = 1; geometric imperfections are neglected.

u = u11*cos(pi*x/a)*sin(pi*y/b);
uh = -pi/a*1/2*w11*sin(pi*y/b)*1/2*w11*sin(pi*y/b)*sin(2*pi*x/a); % additional terms

v = v11*sin(pi*x/a)*cos(pi*y/b);
vh = -pi/b*1/2*w11*sin(pi*x/a)*1/2*w11*sin(pi*x/a)*sin(2*pi*y/b); % additional terms

w = w11*sin(pi*x/a)*sin(pi*y/b);

u = u + uh;
v = v + vh;

%% Velocities and kinetic energy
% The kinetic energy is calculated and the modal mass is neglected.

syms du dv dw

du = diff(u,t);
dv = diff(v,t);
dw = diff(w,t);

T = simplify(int(int(du^2 + dv^2 + dw^2,x,0,a),y,0,b),'Steps',100);


%% Equations of Motion
% In order to evaluate the derivatives w.r.t. to the generalized
% coordinates and velocities for the Lagrange equations, the time-dependent
% symbolic variables have to replaced by time-independent variants.
% d/dt dT/ddq - dT/dq + ... = Q

syms u1 v1 w1
syms du1 dv1 dw1
syms ddu1 ddv1 ddw1

T_s = subs(T,{u11 v11 w11 diff(u11,t) diff(v11,t) diff(w11,t)},{u1 v1 w1 du1 dv1 dw1});

dTddq = [diff(T_s,du1); diff(T_s,dv1); diff(T_s,dw1)];  % derivatives w.r.t. generalized coordinates
dTdq = [diff(T_s,u1); diff(T_s,v1); diff(T_s,w1)];

dTddq = subs(dTddq,{u1 v1 w1 du1 dv1 dw1},{u11 v11 w11 du11 dv11 dw11}); % substitute by time-dependent variables again

dTddqdt = diff(dTddq,t);
dTddqdt = subs(dTddqdt,{diff(du11,t) diff(dv11,t) diff(dw11,t) diff(u11,t) diff(v11,t) diff(w11,t)},{ddu1 ddv1 ddw1 du1 dv1 dw1});
dTddqdt = subs(dTddqdt,{u11 v11 w11 du11 dv11 dw11}, {u1 v1 w1 du1 dv1 dw1});

% The complete equations of motion containing elastic energy and external
% forces have to solved for the accelerations.