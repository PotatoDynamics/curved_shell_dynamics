% This script is only for plotting shape functions.
%% Shape Functions

a = sym('a'); b = sym('b');
x = sym('x');  y = sym('y');

M = 3; N = 3;
Mw = 1; Nw = 1;
M0 = 1; N0 = 1;

U = sym('u%d%d', [M N]);
V = sym('v%d%d', [M N]);
W = sym('w%d%d', [Mw Nw]);      % modal amplitudes (time-dependent)
W0 = sym('A%d%d', [M0 N0]);

u = 0; v = 0; w = 0;

for m = 1:2:M
     for n = 1:2:N
         if n <= Nw && m <= Mw
            u = u + U(m,n)*cos(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT
            v = v + V(m,n)*sin(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
            w = w + W(m,n)*sin(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT

        else
            u = u + U(m,n)*cos(m*pi*x/a)*sin(n*pi*y/b);     % CORRECT
            v = v + V(m,n)*sin(m*pi*x/a)*cos(n*pi*y/b);     % CORRECT
        end     
    end
end

%% Dimensions and Resolution

a = 0.1;
b = 0.1;

x_grid = [0:0.005:0.1];
y_grid = [0:0.005:0.1];

%% Plots

for i = 1:4
    figure(i)
    q_u = zeros(22,1);
    q_u(1+(i-1)*2) = 1;
    surf(x_grid,y_grid,u_plot_22(q_u,x_grid',y_grid))
    xlabel('x')
    ylabel('y')
    
    q_w = zeros(22,1);
    q_w(18+i) = 1;
    figure(i+4)
    surf(x_grid,y_grid,w_plot_22(q_w,x_grid',y_grid))
    xlabel('x')
    ylabel('y') 
end

 