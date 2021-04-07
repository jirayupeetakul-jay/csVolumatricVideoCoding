
function [y_q rate] = SDPC_Coding(y, quantizer_bitdepth, num_rows,num_cols, Phi)

y_max = max(y(:));  y_min = min(y(:));

q = (y_max - y_min)/2^quantizer_bitdepth;

DC_block = ones(256,1)*128;

DC_Measure = Phi*DC_block;

[yq y_q y_index] = SDPC_Encode(y,q,DC_Measure);

total_pixels = num_rows*num_cols;

rate = Measurement_Entropy(yq(:),total_pixels);

end

