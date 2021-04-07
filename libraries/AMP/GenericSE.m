function sigma0_niter=GenericSE(delta,rho,eta,Px,sigma_noise,niter)
% GenericSE State evolution for AMP in generic case 
% (State evolution after normalizing A)
%   sigma0_niter=GenericSE_R(delta,rho,eta,Px,sigma_noise,niter)
%   The goal of this function is to characterize the state evolution of AMP,
%   and it outputs sigma_0, sigma_1, ... , sigma_niter.
%  
%   x_0i ~ (1-epsilon)*delta_0+epsilon*G
%** Notice: If the distribution G is a discrete distribution, and takes on 
%   more than 500 different values, say, 1000 different values, then use
%   the command  set(0,'RecursionLimit',1000) to change the recursion
%   limit. However, be aware that exceeding your available stack space can 
%   crash MATLAB and/or your computer.
%
%
%  Input Arguments:
%   delta:  the undersampling factor, delta=n/N
%   rho:    the relative signal sparsity
%   eta:    the function (handle) eta that we use in GenericAMP
%   Px:     the struct that characterizes the distribution of G, also the distribution of x_0i
%   sigma_noise: the standard deviation of the random noise
%   niter:  the maximum number of iterations
%
%   Px has four fields: Px.option, Px.dist, Px.param, Px.pdf
%   Px.option indicates the option to input the distribution of G:
%   Note: x_0i ~ (1-epsilon)*delta_0+epsilon*G, we input the distribution G. 
%   Px.option=1 to input as a string and parameters like 'Gaussian', ( mu , sigma^2 );
%   Px.option=2 to input the distibution as 'PointMasses' sum ( w_i * delta_(mu_i) );
%   Px.option=3 to input the p.d.f. (function handle) of G. 
%** (IMPORTANT NOTICE: If the p.d.f. equals 0 on some interval, you should 
%   input the p.d.f. as a piecewise function. DO NOT OMIT THIS PART.)
%
%   For Px.option=1, distribution parameters:
%   'Gaussian'   [mean,variance]
%   'Uniform'    [a,b] (support) -- Uniform distribution (continuous)
%   '2Point'     [mu]  x_0i ~ (1-epsilon)*delta_0+epsilon*delta_mu
%                      ( x_0i takes on 0 with probability (1-epsilon), and 
%                       takes on mu with probability epsilon. )
%   '3Point'     [mu]  x_0i ~ (1-epsilon)*delta_0+epsilon/2*delta_mu+epsilon/2*delta_-mu
%                      ( x_0i takes on 0 with probability (1-epsilon), 
%                       takes on mu with probability epsilon/2, and takes 
%                       on -mu with probability epsilon/2. )
%   'Signs'      []    G ~ 0.5*delta_-1+0.5*delta_1
%                      ( G takes on -1 and 1, both with probability 0.5 )
%                      No parameter is required.
%   'Power'      [k]   G takes on 1, 1/2, 1/3, ... , 1/k with probability
%                      1/N, and takes on 0 with probability 1-k/N. 
%   'Beta'       [shape_alpha,shape_beta]
%   'Binomial'   [n,p]  where n denotes the number of trials, p denotes the
%                       success probability in each trial.
%   'Chi-Square' [df]
%   'TypeIGumbel'[location,scale]
%   'Exponential'[rate]  rate=1/mean
%   'Gamma'      [shape,rate]  rate=1/scale
%   'F'          [df1,df2]
%   'Geometric'  [p]  where p is the success probability in each trial. 
%                IMPORTANT: support is {0,1,2,...}. Check help of geopdf
%                for the Geometric p.m.f. in MATLAB.
%   'GEV'        [location,scale,shape]
%   'Lognormal'  [mu,sigma]
%   't'          [df]
%   'Weibull'    [lambda,k]
%   'Unid'       [a,b] (support is the set {a,a+1,...,b}) -- Uniform distribution (discrete)
%   'Rayleigh'   [sigma]
%   'Generalized Pareto'   [location,scale,shape]
%   'Poisson'    [lambda]
%   'Noncentral Chi-Square' [df,mu]
%   'Noncentral t'          [df,mu]
%   'Noncentral F'          [df_1,df_2,lambda]
%   'Negative Binomial'     [r,p]
%   'Hypergeometric'        [A,B,N]  choose k elements from A, and choose
%                                      N-k elements from B
%   Please specify the p.m.f. (option 2) or p.d.f. (option 3) of the distribution other than those listed above. 
%   Example:  Px.option=1;  Px.dist='Gaussian';  Px.param=[10;2];
%
%   For option=2, Px.param=[X_w,X_x]; (Px.dist is not required.)
%   where X_x is the column vector consisting of the values that the
%   discrete random variable takes on, and X_w is the corresponding weight
%   vector (column vector).
%   Example:  Px.option=2;   X_x=[1;2;3;4;5;6;7;8];
%             X_w=[0.01;0.003;0.2;0.6;0.0001;0.002;0.1;0.0849];
%             Px.param=[X_w,X_x];
%
%   For option=3, input Px.pdf as a function handle.
%   Example:  Px.option=3;  Px.pdf=@normpdf;
%** (IMPORTANT NOTICE: If the p.d.f. equals 0 on some interval, you should 
%   input the p.d.f. as a piecewise function. DO NOT OMIT THIS PART.)
%
% 
%  Output Arguments:
%   sigma0_niter:  [sigma_0; sigma_1; ... ; sigma_niter]

fprintf('\nTHEORETICAL ANALYSIS ...\n\n');

if nargin<4
    error('GenericSE: Not Enough Input Arguments.');
end

% Default input parameters
if nargin<5
    sigma_noise=0;
end
if nargin<6
    niter=100;
end

epsilon=rho*delta;
sigma0_niter=zeros(niter+1,1);
    
switch Px.option
    case 1
        if strcmpi(Px.dist,'Gaussian') || strcmpi(Px.dist,'Normal')
            X_mu=Px.param(1);
            X_sigma=sqrt(Px.param(2));
            E_xi2=(X_mu^2+X_sigma^2)*epsilon;
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_Gaussian=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*(normpdf(x,X_mu,X_sigma)).*normpdf(z);
                Expectation_1=integral2(fun_integral_Gaussian,-inf,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end            
        elseif strcmpi(Px.dist,'Uniform') || strcmpi(Px.dist,'Unif') || strcmpi(Px.dist,'UniformContinuous') || strcmpi(Px.dist,'UnifC')
            X_lower=Px.param(1);
            X_upper=Px.param(2);
            E_xi2=epsilon*(((X_lower+X_upper)^2)/4+((X_upper-X_lower)^2)/12);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_Uniform=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*(normpdf(z)).*(unifpdf(x,X_lower,X_upper));
                Expectation_1=integral2(fun_integral_Uniform,X_lower,X_upper,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'2Point') || strcmpi(Px.dist,'2Points')
            X_mu=Px.param;
            E_xi2=epsilon*(X_mu^2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_2Point=@(z) ((1-epsilon)*(eta(sigma_temp*z,sigma_temp)).^2+epsilon*(eta(X_mu+sigma_temp*z,sigma_temp)-X_mu).^2).*normpdf(z);
                % Expectation=integral(fun_integral_2Point,-inf,inf);
                Expectation=integral(fun_integral_2Point,-30,30);
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'3Point') || strcmpi(Px.dist,'3Points')
            X_mu=Px.param;
            E_xi2=epsilon*(X_mu^2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_3Point=@(z) ((1-epsilon)*(eta(sigma_temp*z,sigma_temp)).^2+epsilon/2*((eta(sigma_temp*z+X_mu,sigma_temp)-X_mu).^2+(eta(sigma_temp*z-X_mu,sigma_temp)+X_mu).^2)).*normpdf(z);
                Expectation=integral(fun_integral_3Point,-inf,inf);
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Signs') || strcmpi(Px.dist,'Sign')
            X_x=[-1;1];
            X_w=[0.5;0.5];
            E_xi2=epsilon;
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            l=2;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                summation=@(z) 0;
                for iterl=1:l
                    X_w_temp=X_w(iterl);
                    X_x_temp=X_x(iterl);
                    m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                    m2=@(z) m1(z).*m1(z);
                    sum_temp=@(z) X_w_temp*m2(z);
                    summation=@(z) summation(z)+sum_temp(z);
                end
                fun_integral_Signs=@(z) summation(z).*normpdf(z);
                Expectation_1=integral(fun_integral_Signs,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Power')
            X_k=Px.param;
            X_prob_0=1-X_k/N;
            X_w=(ones(X_k+1,1))./N;
            X_w(1)=X_prob_0;
            X_x=[0;(1./(1:k))'];
            E_xi2=(sum((X_x.^2).*X_w))*epsilon;
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            l=length(X_x);
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                summation=@(z) 0;
                for iterl=1:l
                    X_w_temp=X_w(iterl);
                    X_x_temp=X_x(iterl);
                    m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                    m2=@(z) m1(z).*m1(z);
                    sum_temp=@(z) X_w_temp*m2(z);
                    summation=@(z) summation(z)+sum_temp(z);
                end
                fun_integral_Power=@(z) summation(z).*normpdf(z);
                Expectation=integral(fun_integral_Power,-inf,inf);
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Beta')
            X_a=Px.param(1);
            X_b=Px.param(2);
            E_xi2=epsilon*((X_a/(X_a+X_b))^2+X_a*X_b/((X_a+X_b)^2)/(X_a+X_b+1));
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_Beta=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*betapdf(x,X_a,X_b).*normpdf(z);
                Expectation_1=integral2(fun_integral_Beta,0,1,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Bino') || strcmpi(Px.dist,'Bin') || strcmpi(Px.dist,'Binomial')
            X_n=Px.param(1);
            X_p=Px.param(2);
            X_w=binopdf(0:X_n,X_n,X_p);
            X_x=0:X_n;
            E_xi2=(sum((X_x.^2).*X_w))*epsilon;
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            l=length(X_x);
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                summation=@(z) 0;
                for iterl=1:l
                    X_w_temp=X_w(iterl);
                    X_x_temp=X_x(iterl);
                    m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                    m2=@(z) m1(z).*m1(z);
                    sum_temp=@(z) X_w_temp*m2(z);
                    summation=@(z) summation(z)+sum_temp(z);
                end
                fun_integral_Binomial=@(z) summation(z).*normpdf(z);
                Expectation_1=integral(fun_integral_Binomial,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Chi2') || strcmpi(Px.dist,'Chi-Square') || strcmpi(Px.dist,'Chi-Squared') || strcmpi(Px.dist,'ChiSquare') || strcmpi(Px.dist,'ChiSquared')
            X_df=Px.param;
            E_xi2=epsilon*X_df*(X_df+2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_ChiSquare=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*chi2pdf(x,X_df).*normpdf(z);
                Expectation_1=integral2(fun_integral_ChiSquare,0,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'EV') || strcmpi(Px.dist,'ExtremeValue') || strcmpi(Px.dist,'Type1Extreme') || strcmpi(Px.dist,'TypeIExtreme') || strcmpi(Px.dist,'Gumbel') || strcmpi(Px.dist,'Type1Gumbel') || strcmpi(Px.dist,'TypeIGumbel')
            X_mu=Px.param(1);
            X_beta=Px.param(2);
            euler_gamma=0.577215664901533;
            E_xi2=epsilon*((pi^2)*X_beta^2/6+(X_mu+X_beta*euler_gamma)^2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_TypeIGumbel=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*evpdf(x,X_mu,X_beta).*normpdf(z);
                Expectation_1=integral2(fun_integral_TypeIGumbel,-inf,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Exp') || strcmpi(Px.dist,'Exponential')
            X_lambda=Px.param;
            X_mu=1./X_lambda;
            E_xi2=epsilon*2*(X_mu^2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_Exp=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*exppdf(x,X_mu).*normpdf(z);
                Expectation_1=integral2(fun_integral_Exp,0,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Gamma') || strcmpi(Px.dist,'Gam')
            X_alpha=Px.param(1);
            X_beta=Px.param(2);
            X_k=X_alpha;
            X_theta=1./X_beta;
            E_xi2=epsilon*X_alpha*(X_alpha+1)/(X_beta^2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_Gamma=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*gampdf(x,X_k,X_theta).*normpdf(z);
                Expectation_1=integral2(fun_integral_Gamma,0,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'F') || strcmpi(Px.dist,'F-distribution') || strcmpi(Px.dist,'Fdistribution') || strcmpi(Px.dist,'F-dist') || strcmpi(Px.dist,'Fdist')
            X_df1=Px.param(1);
            X_df2=Px.param(2);
            assert(X_df2>4,'ERROR: df2 should be greater than 4 to have second moment.');
            E_xi2=epsilon*(1/((1-2/X_df2)^2)+2*(X_df2^2)*(X_df1+X_df2-2)/X_df1/((X_df2-2)^2)/(X_df2-4));
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_F=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*fpdf(x,X_df1,X_df2).*normpdf(z);
                Expectation_1=integral2(fun_integral_F,0,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'GEV') || strcmpi(Px.dist,'GeneralizedExtremeValue') || strcmpi(Px.dist,'Generalized EV')
            X_location=Px.param(1);
            X_scale=Px.param(2);
            X_shape=Px.param(3);
            assert(X_shape<0.5,'ERROR: shape parameter xi should be less than 1/2 to have finite second moment.');
            euler_gamma=0.577215664901533;
            if abs(X_shape)<eps
                E_xi2=((X_location+X_scale*euler_gamma)^2+(X_scale^2)*(pi^2)/6)*epsilon;
            else
                E_xi2=(X_location^2+(2*X_scale*(X_shape*X_location-X_scale)*gamma(1-X_shape)+X_scale^2-2*X_location*X_scale*X_shape+gamma(1-2*X_shape)*(X_scale^2))/(X_shape^2))*epsilon;
            end
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_GEV=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*gevpdf(x,X_shape,X_scale,X_location).*normpdf(z);
                Expectation_1=integral2(fun_integral_GEV,-inf,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Lognormal') || strcmpi(Px.dist,'Log-Normal') || strcmpi(Px.dist,'logn') || strcmpi(Px.dist,'Galton')
            X_mu=Px.param(1);
            X_sigma=Px.param(2);
            E_xi2=epsilon*exp(2*(X_mu+X_sigma^2));
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_LogNormal=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*lognpdf(x,X_mu,X_sigma).*normpdf(z);
                Expectation_1=integral2(fun_integral_LogNormal,0,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'t-distribution') || strcmpi(Px.dist,'t') || strcmpi(Px.dist,'tDistribution') || strcmpi(Px.dist,'t-dist') || strcmpi(Px.dist,'tdist')
            X_df=Px.param;
            assert(X_df>2,'ERROR: df should be greater than 2 to have finite second moment.');
            E_xi2=epsilon/(1-2/X_df);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_t=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*tpdf(x,X_df).*normpdf(z);
                Expectation_1=integral2(fun_integral_t,-inf,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'wbl') || strcmpi(Px.dist,'Weibull') || strcmpi(Px.dist,'Weib')
            X_lambda=Px.param(1);
            X_k=Px.param(2);
            E_xi2=epsilon*(X_lambda^2)*gamma(1+2/X_k);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_Weibull=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*wblpdf(x,X_lambda,X_k).*normpdf(z);
                Expectation_1=integral2(fun_integral_Weibull,0,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'UniformDiscrete') || strcmpi(Px.dist,'UniformD') || strcmpi(Px.dist,'unid')
            X_a=Px.param(1);
            X_b=Px.param(2);
            X_n=X_b-X_a+1;
            X_w=ones(X_n,1)./X_n;
            X_x=(X_a:X_b)';
            E_xi2=(sum((X_x.^2).*X_w))*epsilon;
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            l=length(X_x);
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                summation=@(z) 0;
                for iterl=1:l
                    X_w_temp=X_w(iterl);
                    X_x_temp=X_x(iterl);
                    m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                    m2=@(z) m1(z).*m1(z);
                    sum_temp=@(z) X_w_temp*m2(z);
                    summation=@(z) summation(z)+sum_temp(z);
                end
                fun_integral_UniformDiscrete=@(z) summation(z).*normpdf(z);
                Expectation_1=integral(fun_integral_UniformDiscrete,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Rayleigh') || strcmpi(Px.dist,'Rayl')
            X_sigma=Px.param;
            E_xi2=epsilon*2*(X_sigma^2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_Rayleigh=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*raylpdf(x,X_sigma).*normpdf(z);
                Expectation_1=integral2(fun_integral_Rayleigh,0,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Generalized Pareto') || strcmpi(Px.dist,'GeneralizedPareto') || strcmpi(Px.dist,'gp')
            X_location=Px.param(1);
            X_scale=Px.param(2);
            X_shape=Px.param(3);
            assert(X_shape<0.5,'ERROR: shape parameter xi should be less than 1/2 to have second moment.');
            E_xi2=epsilon*((X_location+X_scale/(1-X_shape))^2+(X_scale^2)/((1-X_shape)^2)/(1-2*X_shape));
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_GeneralizedPareto=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*gppdf(x,X_shape,X_scale,X_location).*normpdf(z);
                Expectation_1=integral2(fun_integral_GeneralizedPareto,-inf,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Poisson') || strcmpi(Px.dist,'Poi')
            X_lambda=Px.param;
            X_x_lower=max(0,floor(X_lambda-40*sqrt(X_lambda)));
            X_x_upper=ceil(X_lambda+40*sqrt(X_lambda));
            X_x=(X_x_lower:1:X_x_upper)';
            l=length(X_x);
            X_w=poisspdf(X_x,X_lambda);
            E_xi2=epsilon*X_lambda*(1+X_lambda);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration); 
                summation=@(z) 0;
                for iterl=1:l
                    X_w_temp=X_w(iterl);
                    X_x_temp=X_x(iterl);
                    m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                    m2=@(z) m1(z).*m1(z);
                    sum_temp=@(z) X_w_temp*m2(z);
                    summation=@(z) summation(z)+sum_temp(z);
                end
                fun_integral_Poisson=@(z) summation(z).*normpdf(z);
                Expectation_1=integral(fun_integral_Poisson,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'Noncentral Chi-Square') || strcmpi(Px.dist,'Non-central Chi-Square') || strcmpi(Px.dist,'ncx2')
            X_df=Px.param(1);
            X_mu=Px.param(2);
            E_xi2=epsilon*((X_df+X_mu)^2+2*(X_df+2*X_mu));
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_NoncentralChiSquare=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*ncx2pdf(x,X_df,X_mu).*normpdf(z);
                Expectation_1=integral2(fun_integral_NoncentralChiSquare,0,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end                
        elseif strcmpi(Px.dist,'Noncentral t') || strcmpi(Px.dist,'Non-central t') || strcmpi(Px.dist,'nct')
            X_df=Px.param(1);
            X_mu=Px.param(2);
            assert(X_df>2,'ERROR: df should be greater than 2 to have second moment.');
            E_xi2=epsilon*(1+X_mu^2)/(1-2/X_df);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_NonCentralt=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*nctpdf(x,X_df,X_mu).*normpdf(z);
                Expectation_1=integral2(fun_integral_NonCentralt,-inf,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end               
        elseif strcmpi(Px.dist,'Noncentral F') || strcmpi(Px.dist,'Non-central F') || strcmpi(Px.dist,'ncf')
            X_df1=Px.param(1);
            X_df2=Px.param(2);
            X_lambda=Px.param(3);
            assert(X_df2>4,'ERROR: df2 should be greater than 4 to have second moment.');
            E_xi2=epsilon*((X_df2/X_df1)^2)*(X_lambda^2+2*X_lambda*X_df1+4*X_lambda+X_df1^2+2*X_df1)/(X_df2-2)/(X_df2-4);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                fun_integral_NonCentralF=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*ncfpdf(x,X_df1,X_df2,X_lambda).*normpdf(z);
                Expectation_1=integral2(fun_integral_NonCentralF,-inf,inf,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end       
        elseif strcmpi(Px.dist,'Geometric') || strcmpi(Px.dist,'Geo')
            X_prob=Px.param;
            X_x=(0:floor(-30/log(1-X_prob)))';
            X_w=geopdf(X_x,X_prob);
            E_xi2=epsilon*(1-X_prob)*(2-X_prob)/(X_prob^2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            l=length(X_x);
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                summation=@(z) 0;
                for iterl=1:l
                    X_w_temp=X_w(iterl);
                    X_x_temp=X_x(iterl);
                    m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                    m2=@(z) m1(z).*m1(z);
                    sum_temp=@(z) X_w_temp*m2(z);
                    summation=@(z) summation(z)+sum_temp(z);
                end
                fun_integral_Geometric=@(z) summation(z).*normpdf(z);
                Expectation_1=integral(fun_integral_Geometric,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end
        elseif strcmpi(Px.dist,'NB') || strcmpi(Px.dist,'Negative Binomial') || strcmpi(Px.dist,'NegativeBinomial') || strcmpi(Px.dist,'nbin')
            X_r=Px.param(1);
            X_p=Px.param(2);
            X_x=(0:(ceil(X_r*(1-X_p)/X_p+10*sqrt(X_r*(1-X_p)/X_p^2))))';
            X_w=nbinpdf(X_x,X_r,X_p);
            E_xi2=epsilon*X_p*X_r*(1+X_p*X_r)/((1-X_p)^2);
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            l=length(X_x);
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                summation=@(z) 0;
                for iterl=1:l
                    X_w_temp=X_w(iterl);
                    X_x_temp=X_x(iterl);
                    m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                    m2=@(z) m1(z).*m1(z);
                    sum_temp=@(z) X_w_temp*m2(z);
                    summation=@(z) summation(z)+sum_temp(z);
                end
                fun_integral_NB=@(z) summation(z).*normpdf(z);
                Expectation_1=integral(fun_integral_NB,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end           
        elseif strcmpi(Px.dist,'Hypergeometric') || strcmpi(Px.dist,'HyperGeo') || strcmpi(Px.dist,'hyge')
            X_A=Px.param(1);
            X_B=Px.param(2);
            X_N=Px.param(3);
            X_M=X_A+X_B;
            X_K=X_A;
            X_x=(max(0,X_N-X_B):1:min(X_N,X_A))';
            X_w=hygepdf(X_x,X_M,X_K,X_N);
            E_xi2=(sum((X_x.^2).*X_w))*epsilon;
            sigma0=sqrt(sigma_noise^2+E_xi2/delta);
            sigma0_niter(1)=sigma0;
            l=length(X_x);
            for iteration=1:niter
                disp(['iteration = ' num2str(iteration)]);
                sigma_temp=sigma0_niter(iteration);
                summation=@(z) 0;
                for iterl=1:l
                    X_w_temp=X_w(iterl);
                    X_x_temp=X_x(iterl);
                    m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                    m2=@(z) m1(z).*m1(z);
                    sum_temp=@(z) X_w_temp*m2(z);
                    summation=@(z) summation(z)+sum_temp(z);
                end
                fun_integral_Hypergeometric=@(z) summation(z).*normpdf(z);
                Expectation_1=integral(fun_integral_Hypergeometric,-inf,inf);
                fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
                Expectation_2=integral(fun_integral_0,-inf,inf);
                Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
                sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
            end            
        else
            error('Please specify the distribution of x_0i using p.m.f. (option 2) or p.d.f. (option 3).');
        end
        
    case 2
        X_w=Px.param(:,1);
        X_x=Px.param(:,2);
        E_xi2=(sum((X_x.^2).*X_w))*epsilon;
        sigma0=sqrt(sigma_noise^2+E_xi2/delta);
        sigma0_niter(1)=sigma0;
        l=length(X_x);
        for iteration=1:niter
            disp(['iteration = ' num2str(iteration)]);
            sigma_temp=sigma0_niter(iteration);
            summation=@(z) 0;
            for iterl=1:l
                X_w_temp=X_w(iterl);
                X_x_temp=X_x(iterl);
                m1=@(z) eta(X_x_temp+sigma_temp*z,sigma_temp)-X_x_temp;
                m2=@(z) m1(z).*m1(z);
                sum_temp=@(z) X_w_temp*m2(z);
                summation=@(z) summation(z)+sum_temp(z);
            end
            fun_integral_Pointmasses=@(z) summation(z).*normpdf(z);
            Expectation_1=integral(fun_integral_Pointmasses,-inf,inf);
            fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
            Expectation_2=integral(fun_integral_0,-inf,inf);
            Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
            sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
        end       
        
    case 3
        X_pdf=Px.pdf;
        E_xi2_fh=@(x) (x.^2).*X_pdf(x);
        E_xi2=(integral(E_xi2_fh,-inf,inf))*epsilon;
        sigma0=sqrt(sigma_noise^2+E_xi2/delta);
        sigma0_niter(1)=sigma0;
        for iteration=1:niter
            disp(['iteration = ' num2str(iteration)]);
            sigma_temp=sigma0_niter(iteration);
            fun_integral_pdf=@(x,z) ((eta(x+sigma_temp*z,sigma_temp)-x).^2).*(X_pdf(x)).*normpdf(z);
            Expectation_1=integral2(fun_integral_pdf,-inf,inf,-inf,inf);
            fun_integral_0=@(z) ((eta(sigma_temp*z,sigma_temp)).^2).*normpdf(z);
            Expectation_2=integral(fun_integral_0,-inf,inf);
            Expectation=(1-epsilon)*Expectation_2+epsilon*Expectation_1;
            sigma0_niter(iteration+1)=sqrt(Expectation/delta+sigma_noise^2);
        end
        
end

end