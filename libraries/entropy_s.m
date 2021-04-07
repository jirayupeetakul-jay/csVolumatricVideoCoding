function H = entropy_s(I)
% Assume I in the range 0..1
p = hist(I(:), linspace(0,1,256));  % create histogram
p(p==0) = []; % remove zero entries that would cause log2 to return NaN
p = p/numel(I); % normalize histogram to unity
H = -sum(p.*log2(p)); % entropy definition

