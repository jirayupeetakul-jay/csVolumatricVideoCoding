function H = Measurement_Entropy(y, total_pixels)

% assumed that y is vector.
y = y(:);

nbin = unique(y);
frequency = hist(y,nbin);

p = frequency/sum(frequency);
p(p==0) = 1;
H = sum(-p.*log2(p));

% total pixel is used to calculate the real entropy
H = H*numel(y)/total_pixels;