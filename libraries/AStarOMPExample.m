% script to demonstrate use of AStarOMP function.
% creates a random sparse vector x and a random observation matrix Phi.
% computes random observation of x as y = Phi*x
% reconstructs x from y via AStarOMP.
%
% Nazim Burak Karahanoglu, Sabancý University, 09.2012
%   email: karahanoglu at sabanciuniv.edu


%clc
% data properties
K = 45;     % sparsity
N = 256;    % # all dimensions 
M = 64;    % # observations

% create a random sparse vector
display('Creating random sparse test vector x...');
x = zeros(N,1);
x_ind= randi(N,K,1);   %K non zero entries
% ensure K-sparseness by cancelling repetations of the same index
done = 0;
while done == 0
    done = 1;
    x_ind = sort(x_ind);
    for i = 1:K-1
        if x_ind(i) == x_ind(i+1)
            x_ind(i+1) = randi(N,1,1);
            done = 0;
        end;
    end;
end;
x(x_ind) = randn(K,1);          % standard Gaussian entries

% random observation matrix
display('Creating random observation matrix Phi...');
Phi = randn(M,N);
% normalize observation matrix.
for k=1:N
    Phi(:,k) = Phi(:,k)/norm(Phi(:,k));
end;

% observed vector
y = Phi*x;

%A*OMP options
%options, may be omitted
options.I = 3;  % I>=1
options.B = 2;  % B>=1
options.P = 200;   % P>=1
options.beta = 1.2;
options.AuxMode = 'AdapMul'; %possible: 'Add','Adap', 'Mul', 'AdapMul'
options.alpha = 0.97;
options.Control = 'Off';  %possible: 'On'/'Off'
options.Display = 'Off';  %possible: 'On'/'Off'
options.Kmax = M/2;
options.eps = 1e-7;  %possible: 'On'/'Off'
% run A*OMP with residue based termination
%%
display('Running A*OMP...');
[xhat,iter] = AStarOMP(y, Phi, N, options);
% % Alternative: run A*OMP with termination based on sparsity level
% options.Kmax = K;
% options.eps = 0; 
% tic
% xhat = AStarOMP(y, Phi, N,options);
% toc
%% display results 1
xhat_ind = find(xhat);
display('RECONSTRUCTION RESULT:')
if (isequal(x_ind, xhat_ind))
    display(' x is exactly reconstructed.');
else
    noIdentifiedComp = sum(ismember(xhat_ind, x_ind));
    display(' Exact reconstruction failed:') 
    display(['   ' num2str(noIdentifiedComp) ' entires of x were correctly reconstructed.']); 
    display(['   ' num2str(K-noIdentifiedComp) ' entires of x were missed.']); 
end

% plot signals 
MSE = mse(x - xhat);
display([' Mean Squared Error: ' num2str(MSE)]);
% figure; plot(1:N,x,'--b'); hold on; plot(1:N,xhat,'-.r'); plot(1:N,(x-xhat),':g');
% legend('Test signal','Reconstructed signal','Reconstruction Error');

%% display results 1
xhat_ind = find(xhat);
display('RECONSTRUCTION RESULT:')
if (isequal(x_ind, xhat_ind))
    display(' x is exactly reconstructed.');
else
    noIdentifiedComp = sum(ismember(xhat_ind, x_ind));
    display(' Exact reconstruction failed:') 
    display(['   ' num2str(noIdentifiedComp) ' entires of x were correctly reconstructed.']); 
    display(['   ' num2str(K-noIdentifiedComp) ' entires of x were missed.']); 
end

% plot signals 
MSE = mse(x - xhat);
display([' Mean Squared Error: ' num2str(MSE)]);
figure; plot(1:N,x,'--b'); hold on; plot(1:N,xhat,'-.r'); plot(1:N,(x-xhat),':g');
legend('Test signal','Reconstructed signal','Reconstruction Error');