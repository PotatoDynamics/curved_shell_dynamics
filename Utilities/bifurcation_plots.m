% This script reads time series from the bifurcation experiment and
% generates time series, phase portrait and spectrum plots.
%% Parameters
n = 223;    % number of file to be read
n_dof = 9;
wn = 952.26;

A = importdata(sprintf('Bifurcation_%d.csv',n));
disp(sprintf('Solution for f = %dN loaded.',1.75*n));
T = A(1,:)/wn;
Y = A(2:19,:);
Y(10:18,:) = 1/(2*pi)*Y(10:18,:);

%% Time Series
q = 9;  % generalized coordinate to be plotted; 9: w_1,1

figure(1)
plot(T(end-1010:end)-T(end-1010)+1,Y(q,end-1010:end)); hold on
% plot(T(end-1000:end),cos(2*pi*0.8*wn*T(end-1000:end)))


%% Stroboscopic Map (pedestrian version, inaccurate)

F = cos(2*pi*0.8*T*wn);
[pks,ind] = findpeaks(F);
% If the sampling frequency in simulation is not a multiple of the excitation
% frequency, the peaks identified by 'findpeaks' do not exactly match the
% real locations as there is not sample at the true maximum.

x = Y(q,ind(1:end));
dx = Y(q+n_dof,ind(1:end));

figure
plot(x,dx,'k.','MarkerSize',4)
title('Stroboscopic map for $f_N=1379\,N$.\\Sampling Frequency $F_S=0.8\,\omega_{1,1}$.','interpreter','latex')
xlabel('$w_{1,1}/h$','interpreter','latex')
ylabel('$\dot{w}_{1,1}/h\omega_{1,1}$','interpreter','latex')
% figure(2)
% plot(T(end-1000:end)-T(end-1000)+1,Y(1,end-1000:end))
% 
% figure(3)
% plot(T(end-1000:end)-T(end-1000)+1,Y(2,end-1000:end))

%% Phase Portrait

figure(2)
plot(Y(q,end-2000:end),Y(q+n_dof,end-2000:end))
title('Phase Portrait for $f_N=1344N$.','interpreter','latex')
xlabel('$w_{1,1}/h$','interpreter','latex')
ylabel('$\dot{w}_{1,1}/h\omega_{1,1}$','interpreter','latex')

%% Fourier Transform/PSD


t = T(end-20000:end);
y = Y(q,end-20000:end)';

% FFT
L = floor(length(t)/2);
Y_f = fft(y,L,1);
Fs = L/(t(L)-t(1));
P2 = 1/L*abs(Y_f);
P1 = P2(1:L/2+1,:);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(1:(L/2))/L;

figure(4)
plot(f/1000,(P1(2:end)))
xlim([0 5])

% PSD
L = floor(length(t)/2);
Y_f = fft(y,L,1);
Y_f = Y_f(1:L/2+1);
Fs = L/(t(L)-t(1));
P1 = 1/(Fs*L)*abs(Y_f).^2;
P1(2:end-1) = 2*P1(2:end-1);
f = 0:Fs/L:Fs/2;

figure(5)
plot(f/1000,10*log10(P1))
xlabel('f/[kHz]')
ylabel('PSD')
xlim([0 5])

