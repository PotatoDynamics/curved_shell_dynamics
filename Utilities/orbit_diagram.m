%% Settings
% This script reads time series from csv-files and generates a bifurcation
% diagram by extracting the unique response amplitudes.
% Each file to be read contains the last 1000 cycles of the corresponding
% force step. At an excitation frequency of 0.8*omega_11 and a sampling
% frequency of 50*omega_11 this corresponds to fsize samples.

N = 800;    % index of last file to be read
n1 = 800;    % last file in the batch
n0 = 751;     % index of first file
dF = 1.75;  % step size of forcing amplitude
dT = 1250;  % number of ca
fs = 50;        % samples per cycle of omega_11
fsize = 62501;  % number of samples in each file
n_dof = 19;     % number of degrees of freedom plus 1 for the time vector
ome11 = 952.26;

%% Make vector of force steps

F = linspace(1.75,1400,800);
amplitudes = cell(2,N);

%% Memory Allocation

A = zeros(n_dof,(n1-n0+1)*fsize);

%% Import Data from file

for i = n0:n1
    A(:,1+(i-n0)*fsize:(i-n0+1)*fsize) = importdata(sprintf('BifRes_%i.csv',i));
end

%% Find unique amplitudes

for i = n0:n1
    [pks,idx] = findpeaks(A(10,1+(i-n0)*fsize:(i-n0+1)*fsize));
    tol = 0.005;
    pks_uni = uniquetol(pks,tol);
    amplitudes{1,i} = F(i);
    amplitudes{2,i} = round(pks_uni,3);
end

%% Save Results

save('Bif_IC_750800_f175_df1400.mat','tol','amplitudes');

%% Plot Results

figure(3)
for i = 1:n1
    plot(amplitudes{1,i},amplitudes{2,i},'k.','MarkerSize',0.5);
    hold on;
    xlabel('$F_{N}$ in N','interpreter','latex');
    ylabel('$\hat{w}_{1,1}/h$','interpreter','latex');
%     xlim([0 1400])
%     ylim([0 3]);
end

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%%
% This section plots the bifurcation diagram, when loaded from the result
% files.

figure(4)
for i = 1:800
    plot(Amplitudes{1,i},Amplitudes{2,i},'k.','MarkerSize',0.5);
    hold on;
    xlabel('$F_{N}$ in N','interpreter','latex');
    ylabel('$\hat{w}_{1,1}/h$','interpreter','latex');
%     xlim([0 1400])
%     ylim([0 3]);
end

%% Plotting
% This part was done by Richard Martin. It makes a pixel graphic for the 
% bifurcation diagram

disp('Plotting...')
larry = nan(1,N);

for i=1:N
    larry(1,i)=length(Amplitudes{2,i});
end

maxlarry = max(larry);
uniquepeaksmat = nan(maxlarry,N);

for i=1:N
    nmissingfields = maxlarry-length(Amplitudes{2,i});
    uniquepeaksmat(:,i) = [Amplitudes{2,i}'; nan(nmissingfields,1)];
end

xmin=0;
xmax=max(uniquepeaksmat(:));

xpeaks=xmin:10^(-3):xmax;
uniquepeaksimg = zeros(length(xpeaks),N);

for i=1:N
    for ii=1:length(xpeaks)
        if ~isempty(find(amplitudes{2,i}==round(xpeaks(ii),3),1))
            uniquepeaksimg(ii,i)=1;
        end
    end
end

imagesc(F,xpeaks,-uniquepeaksimg);colormap gray;axis xy;
xlabel('$F_{N}$','interpreter','latex')
ylabel('$\hat{w}_{1,1}$','interpreter','latex')

%% Vector form
% This section turns the cell array containing the force amplitudes and 
% the unique amplitudes into a vector, so that it can be plotted properly.

num_points = sum(larry,2);
f_vec = nan(num_points,1);
a_vec = nan(num_points,1);
ind = 1;

for i = 1:N
    n = larry(1,i);
    f_vec(ind:ind+n-1) = Amplitudes{1,i}*ones(n,1);
    a_vec(ind:ind+n-1) = Amplitudes{2,i};
    ind = ind + n;
end

A = [f_vec a_vec];
