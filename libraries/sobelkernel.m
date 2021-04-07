function [kernel, s, d] = sobelkernel(size, varargin)
%SOBELKERNEL returns sobel kernel
%   KERNEL=SOBELKERNEL(SIZE) Returns Sobel filter of predefined size.
%
%   KERNEL=SOBELKERNEL(SIZE, 'NORMALISED') The Sobel matrix should be
%   normalised if proper derivative estimator is required.
%
%   [KERNEL, S, D]=SOBELKERNEL(SIZE, 'NORMALISED') Returns also smoothing S
%   and derivative D components individually. The kernel is then a
%   multiplication of these two vectors: S'*D. This can be usefull for
%   convolution acceleration via separable kernels as illustrated at:
%		http://blogs.mathworks.com/steve/2006/10/04/separable-convolution/
%   
%   Example
%   -------
%       sobelkernel(4)
%
%   See also IMFILTER, FSPECIAL.
%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.1 $  $Date: 2013/02/13 16:58:01 $
%   For method description see:
%       http://stackoverflow.com/questions/9567882/sobel-filter-kernel-of-large-size
% Parameter checking.
if ~isempty(varargin) && strcmp(varargin(1), 'normalise')
    normalisation = 1/8;
else
    normalisation = 1;
end
% The dafault 3x3 Sobel kernel.
s = normalisation * [1 2 1];
d = [1 0 -1];
% Convolve the default 3x3 kernel to the desired size.
for i=1:size-3
    s = normalisation * conv([1 2 1], s);
    d = conv([1 2 1], d);
end
% Output the kernel.
kernel = s'*d;
end