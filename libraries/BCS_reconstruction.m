function [reconstructed_image] = BCS_reconstruction(y, phi, theta, image_reconstruction_algorithm, sub_pixels)
%___MINIMIZATION___
    %___SOLVE THE LP (L1-NORM SOLUTION)___
    switch image_reconstruction_algorithm
        case 'l1_eq_pd' %Stable
            MIN_ENERGY          = theta'*y;
            opt_results         = l1eq_pd(MIN_ENERGY, theta, [], y, 1e-3); % L1-magic toolbox
            inv_transform       = ifwht(opt_results);
            reconstructed_image = reshape(inv_transform, sub_pixels, sub_pixels); %raster
        case 'l1qc_logbarrier' %Stable
            MIN_ENERGY          = theta'*y;
            sigma               = 0.005;
            K                   = size(MIN_ENERGY,2);
            ethetalon           = sigma*sqrt(K)*sqrt(1 + 2*sqrt(2)/sqrt(K));
            opt_results         = l1qc_logbarrier(MIN_ENERGY, theta, [], y, ethetalon, 1e-3);
            inv_transform       = ifwht(opt_results);
            reconstructed_image = reshape(inv_transform, sub_pixels, sub_pixels); %raster
        case 'tvdantzig_logbarrier' %Possible with some matrices
            MIN_ENERGY          = theta'*y;
            ethetalon           = 5e-3;
            opt_results         = tvdantzig_logbarrier(MIN_ENERGY, theta, [], y, ethetalon, 1e-3, 5, 1e-8, 1500);
            inv_transform       = ifwht(opt_results);
            reconstructed_image = reshape(inv_transform, sub_pixels, sub_pixels); %raster
        case 'twist2' %Stable
            MIN_ENERGY          = theta'*y;
            %There is a initialize function inside use y instead
            tau                 = 0.001;    
            opt_results         = TwIST(y,theta,tau,'Lambda', 1e-3);
            inv_transform       = ifwht(opt_results);
            reconstructed_image = reshape(inv_transform, sub_pixels, sub_pixels); %raster
        case 'omp' %Not even close
            MIN_ENERGY = theta'*y;
            %%Solves the problem Ab=x s.t. |x|_0<S
            %%Input: measurement matrix A, 
%                    observation b 
%                    sparsity S, this parameter is the most problem because
%                    in general we do not know.
%                    we can cheating by assigning equal to original dimension.
            %%Output x obtained through orthogonal matching pursuit
            opt_results         = OMP(phi,y,length(y));
            inv_transform       = ifwht(opt_results);
            reconstructed_image = reshape(inv_transform, sub_pixels, sub_pixels); %raster
        case 'matching_pursuit' %Not even close
            MIN_ENERGY = theta'*y;
            %%Solves the problem Ab=x s.t. |x|_0<S
            %%Input: measurement matrix A, 
%                    observation b 
%                    sparsity S, this parameter is the most problem because
%                    in general we do not know.
%                    we can cheating by assigning equal to original dimension.
            %%Output x obtained through orthogonal matching pursuit
            opt_results = matching_pursuit(phi,y,30);
        case 'nesta' %Not even close
            MIN_ENERGY               = theta'*y;
            opts                     = [];
            opts.Verbose             = 0;
            opts.tolvar              = 1e-8;
            b                        = y;     % the data
            A                        = orth(phi')';
            At                       = [];
            delta                    = 0;
            muf                      = 1e-8;
            [opt_results,niter,resid,outData] = NESTA(A,At,b,muf,delta,opts);
            inv_transform            = ifwht(opt_results);
            reconstructed_image      = reshape(inv_transform, sub_pixels, sub_pixels); %raster
    end
end