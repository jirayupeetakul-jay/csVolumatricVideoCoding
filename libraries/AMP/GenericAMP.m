function [ empiricaliterwatch_sigma, xhat ] = GenericAMP( y,A,Eta,Etader, par )

% Inputs:
%   y :          observations
%   A :          a function handle that represents matrix A, A(x,1) means
%                A*x, A(x,2) means A'*x 
%   Eta :        a function handle which is a generic denoiser, 
%                xhat=Eta(temp,sigma) 
%   Etader :     a function handle which is the derivative function of the 
%                denoise function Eta, if you can't provide this derivative
%                function, please input "Null"
%   niter :      the maximum number of iterations
%   par:         a cell with two elements, the first denotes whether we need 
%                all the estimates in whole process or just the final estimate 
%                we obtain, "1" means need, "0" means not; the second
%                denotes how many iteration you want AMP to run, if you
%                input a positive integer t, then AMP runs t times
%                iterations for you, if you input the string 'Auto', then
%                AMP will try to runs 100 times iterations and stop when
%                the ratio of l2 norm of (x_(t+1)-x_(t)) and l2 norm of
%                x_(t) less then 0.01


% Related functions: Eta_der_Estimate

% check: @A,@Eta,@Etader, sigma is the estimated std instead of variance.
% sigma_w seems useless in this function


n=length(y);
lengthN=A(zeros(n,1),2);
N=length(lengthN);

pick=randperm(N);
DtmIndex=pick(1:5);
DtmNorms=zeros(5,1);
I=eye(N);
for i=1:5
    DtmNorms(i)=norm( A(I(:,DtmIndex(i)),1) ).^2;
end

if (sum(DtmNorms>1.1)+sum(DtmNorms<0.9))>0 % We need to normalize A matrix
    disp('It is necessary to normalize the A matrix, please wait...');
    tempA=zeros(n,N);
    tempA_ra=zeros(n,N);
    normalize_time_total=0;
    for j=1:N
        t0=cputime;
        tempA(:,j)=A(I(:,j),1);  
        normalize_time=(cputime-t0)/60;
        normalize_time_total=normalize_time_total+normalize_time;
        
        if j/(N/100)==fix(j/(N/100))
            normalize_time_remain=normalize_time_total*(N-j)/j;
            percent=j/(N/100);
            disp(['Normalizing has been through ' num2str(percent) '%.' 10 'The estimated remaining time for Normalzing is ' num2str(normalize_time_remain) ' minutes.']);
        end
        
    end
    
    for j=1:N                              %remove average
        tempA_ra(:,j)=tempA(:,j)-mean(tempA(:,j));
    end
    
    colnormA=(sqrt(sum(tempA_ra.^2,1)))';
    ind=find(colnormA==0);
    colnormA(ind)=(sqrt(sum(abs(tempA(:,ind)).^2,1)))';
    
    disp('Normalizing ends, Iteration starting...');
else
    disp('It is not necessary to normalize the A matrix, Iteration starting...');
    colnormA=ones(N,1);
end


% Denote normalized A matrix as AA, then, when we calculate AA*v, we need 
% to do A*(v./colnormA); when we calculate AA'*v, we need to do (A'*v)./colnormA.

%%%%%

par1=logical(par{1});
par2=par{2};
if ischar(par2)
    niter=100;
else
    niter=par2;
end

empiricaliterwatch_sigma=zeros(niter+1,1);
xall=zeros(N,niter+1);

mx=zeros(N,1);
mz=y-A(mx./colnormA,1);

for iter=1:niter
    disp(['iteration = ' num2str(iter)]);
    temp_z=A(mz,2)./colnormA+mx;
    sigma_hat= norm(mz)/sqrt(n);
    mx=Eta(temp_z,sigma_hat);
    
    if strcmpi(Etader,'Null')
        mz=y-A(mx./colnormA,1)+mz*Eta_der_Estimate(temp_z,sigma_hat,Eta )*N/n;
    else
        mz=y-A(mx./colnormA,1)+mz*sum(Etader(temp_z,sigma_hat))/n;
    end

    empiricaliterwatch_sigma(iter)=sigma_hat;
    xall(:,iter+1)=mx./colnormA;    
    
    if niter==100 && abs(empiricaliterwatch_sigma(iter)-norm(mz)/sqrt(n))<0.001
        break;
    end

end

empiricaliterwatch_sigma(iter+1)=norm(mz)/sqrt(n);

if niter==100
    empiricaliterwatch_sigma=empiricaliterwatch_sigma(1:(iter+1));
    xall=xall(:,1:(iter+1));
end

if iter==100
    fprintf('Iteration reaches the maximum (100) times,\nthe algorithm does not converge within 100 iterations.\n')
end

if par1
    xhat=xall;
else
    xhat=xall(:,end);
end

end






