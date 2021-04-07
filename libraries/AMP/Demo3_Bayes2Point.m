function [empiricaliterwatch_sigma,sigma0_niter,xhat]=Demo3_Bayes2Point(y,A,Eta,Etader,par,sigma_noise,rho)
% [empiricaliterwatch_sigma,sigma0_niter,xhat]=demo3_Bayes2Point(y,A,Eta,Etader,par,sigma_noise,rho)
% Demo function of Generic_AMP.m and GenericSE.m, using Bayes2Point
% setting, the matrix A corresponding to the function handle A is a random
% Gaussian matrix with mean 0 and variance 1/n.
% 
% Input Arguments:
%  y:      observables
%  A:      the function handle that represents matrix A. A(x,1) means A*x,
%          and A(x,2) means A'*x
%  Eta:    the function handle which is a generic denoise function,
%          xhat_temp=Eta(temp,sigma);
%  Etader: the derivative function of Eta, if you cannot provide this
%          function, then please input 'Null' 
%  par:    par indicates whether we need all the estimates in the whole 
%          process or only the final estimate we obtain, 1 means that we need 
%          all the estimates, and 0 means we need the final estimate only.
%  sigma_noise: the standard deviation of the noise
%  rho:    the sparsity of the problem
%
% Output Arguments:
%  empiricaliterwatch_sigma: the standard deviation of the estimate in 
%                            each iteration in empirical setting
%  sigma0_niter: the standard deviation in each iteration in theoretical setting
%  xhat:   the estimate of xhat in each iteration (when par=1) or the
%          final estimate (when par=0)

if nargin<1
    n=1500;
    N=5000;
    A=@SNM;
    Eta=@Bayes2Point_Eta;
    Etader=@Bayes2Point_eta_prime;
    par=cell(2,1);
    par{1}=0; par{2}=10;
    niter=10;
    
    sigma_noise=sqrt(0.1);
    
    Px.option=1;
    Px.dist='2Point';
    mu=10;
    Px.param=mu;
        
    k=100; 
    rho=k/n;

    x0=zeros(N,1);
    pick=randperm(N);  
    x0(pick(1:k))=Px.param; 
    y=A(x0,1)+random('normal',0,sigma_noise,n,1);
    delta=n/N;

end

[empiricaliterwatch_sigma,xhat]=GenericAMP(y,A,Eta,Etader,par);
sigma0_niter=GenericSE(delta,rho,Eta,Px,sigma_noise,niter);

plot(empiricaliterwatch_sigma,'b-');
hold on;
plot(sigma0_niter,'r-');
hold off;
xlabel('iteration');
ylabel('sigma');
legend('Empirical','Theoretical');

end


function [newv]=SNM(v,par)   
n=1500;
N=5000;
rng('default');
M=random('normal',0,1/sqrt(n),n,N);
% M=random('normal',0,1,n,N);
if par<1.5
    newv=M*v;
else
    newv=M'*v;
end
end

function eta=Bayes2Point_Eta(v,sigma)
mu=10; epsilon=100/5000;
sigma_sq=sigma^2;
v_exp_mu=-(v-mu).^2/sigma_sq/2;
v_exp_0=-v.^2/sigma_sq/2;
numerator=mu*epsilon*exp(v_exp_mu);
denominator=(1-epsilon)*exp(v_exp_0)+epsilon*exp(v_exp_mu);
eta=numerator./denominator;
end

function eta_prime=Bayes2Point_eta_prime(v,sigma)
mu=10; epsilon=100/5000;
sigma_sq=sigma^2;
v_exp_mu=-(v-mu).^2/sigma_sq/2;
v_exp_0=-v.^2/sigma_sq/2;
denominator=(1-epsilon)*exp(v_exp_0)+epsilon*exp(v_exp_mu);
prime_numerator=mu^2*epsilon*(1-epsilon)*exp(v_exp_0).*exp(v_exp_mu);
prime_denominator=sigma_sq*(denominator.^2);
eta_prime=prime_numerator./prime_denominator;
end