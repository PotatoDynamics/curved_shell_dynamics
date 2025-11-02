% This script reads time series from the free response and generates a
% backbone curve. The time series is split into batches of 2000 samples,
% for each batch the average oscillation amplitude and frequency are
% determined by FFT.

%% Relevant parameters
ome11 = 952.26; % first linear eigenfrequency
batch = 2000;   % batch size
n = 25;         % number of batches per file (50000 samples/file)
N_up = 216;     % determines which files are loaded
N_low = 0;
N = N_up + N_low;

I_tot = zeros(N*n,1);   % indices for frequencies
MAX_tot = zeros(N*n,1); % maximum deflection
MIN_tot = zeros(N*n,1); % minimum deflection

os = 10; % number of zeros for zero-padding; over-sampling

k = 50;

%% Read data from file and evaluate spectrum and average amplitude

for j = 1:N_up
    
    fname = sprintf('Free_Ama_%i.csv',j);
    A = importdata(fname);
    T = A(1,:)'/ome11;
    Y = A(2:end,:)';
    
    n = length(T)/batch;
    
    P1 = zeros(n*os*batch/2,1);
    I = zeros(n,1);
    minA = zeros(n,1);
    maxA = zeros(n,1);
    M = zeros(n,1);
    
    for i = 1:n
        t = T(1+(i-1)*batch:i*batch);
        y = Y(1+(i-1)*batch:i*batch,9);
        minA(i) = mean(mink(y,k));
        maxA(i) = mean(maxk(y,k));
        L = floor(length(t));
        Y_f = fft(y,10*L,1);
        Fs = L/(t(L)-t(1));
        P2 = abs(Y_f/L);
        P1(1+(i-1)*os*batch/2:i*os*batch/2) = P2(1:os*L/2);
        P1(1+(i-1)*os*batch/2:10+(i-1)*os*batch/2) = 0;
        P1(2+(i-1)*os*batch/2:i*os*batch/2-1) = 2*P1(2+((i-1)*os*batch/2):i*os*batch/2-1);
        [M(i),I(i)] = max(P1(2+((i-1)*os*batch/2):i*os*batch/2-1));
    end
    
    f = Fs*(1:(os*L/2))/(os*L);
    save(strcat(fname(1:end-3),'mat'),'minA','maxA','P1','I')
    I_tot(1+(j-1)*n:j*n) = I;
    MAX_tot(1+(j-1)*n:j*n) = maxA;
    MIN_tot(1+(j-1)*n:j*n) = minA;
      
end

%% Plots

figure
plot(f(I_tot)/ome11,abs(MIN_tot))
hold on
% plot(BB_min(:,1),BB_min(:,2))

figure
plot(f(I_tot)/ome11,abs(MAX_tot))
hold on 
% plot(BB_max(:,1),BB_max(:,2))
hold on 
% plot(BB_max_free(:,1),BB_max_free(:,2),'ro')

%% Spectrum without oversampling/zero-padding (low resolution)

f = Fs*(1:(L/2))/L;
P1 = zeros(n*batch/2,1);
I = zeros(n,1);
minA = zeros(n,1);
maxA = zeros(n,1);
M = zeros(n,1);

for i = 1:n
    t = T(1+(i-1)*batch:i*batch);
    y = Y(1+(i-1)*batch:i*batch,9);
    minA(i) = min(y);
    maxA(i) = max(y);
    L = floor(length(t));
    Y_f = fft(y,L,1);
    Fs = L/(t(L)-t(1));
    P2 = abs(Y_f/L);
    P1(1+(i-1)*batch/2:i*batch/2) = P2(1:L/2);
    P1(2+(i-1)*batch/2:i*batch/2-1) = 2*P1(2+((i-1)*batch/2):i*batch/2-1);
    [M(i),I(i)] = max(P1(2+((i-1)*batch/2):i*batch/2-1));
end

%% 
figure
plot(f(I_tot),abs(MIN_tot));

%% Only one batch for experimentation

t = T(1+(1-1)*batch:1*batch);
y = Y(1+(1-1)*batch:1*batch,9);
L = floor(length(t));
Y_f = fft(y,10*L,1);
Fs = L/(t(L)-t(1));
f = Fs*(0:(5*L))/(10*L);
P2 = abs(Y_f/L);
P1 = P2(1:10*L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
[M,I] = max(P1);

plot(f,P1)

%%
figure
plot(f,P1(2+((1-1)*batch/2):1*batch/2+1))
figure
plot(f(I),minA)

