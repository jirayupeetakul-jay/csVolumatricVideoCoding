function [ empiricaliterwatch_sigma, xhat ] = Demo_1( y,A,Eta,Etader,par )

% This Demo is for noiseless setting, Eta is SoftThresholding function, Etader is
% offered, A is random matrix of which every element has iid N(0,1/n)
% distribution. You can also pretend Etader is not provided and type "Null"
% to run AMP. The sparsity of x0 is k=100. 
%
% Input Arguments:
%  y:      observables
%  A:      the function handle that represents matrix A. A(x,1) means A*x,
%          and A(x,2) means A'*x
%  Eta:    the function handle which is a generic denoise function,
%          xhat_temp=Eta(temp,sigma);
%  Etader: the derivative function of Eta, if you can't provide this
%          funciton, then please input "Null"
%  niter:  the maximum number of iterations
%  par:    indicates whether we need all the estimates in the whole process
%          or only the final estimate we obtain, 1 means that we need all 
%          the estimates, and 0 means we
%          need the final estimate only.
%
% Output Parameters:
%  empiricaliterwatch_sigma: the standard deviation of the estimate in 
%          each iteration in empirical setting
%  xhat:   the estimate of xhat in each iteration (when par=1) 
%          or the final estimate (when par=0)


if nargin<1
    
    %n=1500;
    N=5000;
    A=@SNM;
    Eta=@SoftThreshold;
    Etader=@SoftThresholdDer;
    %Etader='Null'; 
    par=cell(2,1);
    par{1}=0; par{2}='Auto';
    
    k=100; % sparsity
    x0=zeros(N,1);
    pick=randperm(N);  
    x0(pick(1:k))=random('norm',10,1,[k,1]); % generate x0 with element~N(10,1)
    y=A(x0,1);

end

[ empiricaliterwatch_sigma, xhat ] = GenericAMP( y,A,Eta,Etader,par );

plot(empiricaliterwatch_sigma);
xlabel('iteration');
ylabel('sigma');

%disp(xhat);

end

function [newv]=SNM(v,par)   

% Standard normal matrix
    
n=1500;
N=5000;
rng('default');
%M=random('norm',0,1/sqrt(n),n,N); % elements~N(0,1/n)
M=random('norm',0,1,n,N);

if par==1
    newv=M*v;
else
    newv=M'*v;
end

end

function [v_new]=SoftThreshold(vv,tt)

% Soft thresholding function

v_new=(abs(vv)> tt).*(abs(vv)-tt).*sign(vv);

end


function [v_der]=SoftThresholdDer(vv,tt)

% Derivative function of Soft thresholding

v_der=abs(vv)> tt;

end