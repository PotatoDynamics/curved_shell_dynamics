%% Parameters of the Antenna
% Here is a set of parameters for the antenna taken from the JSV paper.
par = struct;
par.rho = 7815.42;
par.E = 201.12e9;
par.nu = 0.31;
par.Rx = 10;
par.Ry = 0.019;
par.a = 0.5;
par.b = 0.019;
par.h = 0.00011;
zeta = 0.004;
par.c = zeta*par.rho*par.h*par.a*par.b/2*[20 100 250];
par.A11 = 0.0;
par.T11 = 1/20;
