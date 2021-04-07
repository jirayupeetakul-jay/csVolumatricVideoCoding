function adc_out = adc(R, totalbits, X)
%
% adc_out = adc(R, totalbits, X)
% R: quantization input range, i.e., [min_in, max_in];
% totalbits: total number of bits to digitize the input signal
% X: analog input vector
% adc_out: ADC outputs
%
sz = size(X);
if sz(1) > sz(2)
  X = X.';
end
q_lev = 2^totalbits;
midPnt = q_lev/2;    % center point
R_max = R(2);
R_min = R(1);
step = (R_max - R_min)/q_lev;
offset = 0.5*step;
if (length(step) > 1 | step <= 0),
  error('Quantization range = [min_value, max_value],  max_value > min_value')
end
% min_max clamping
X(find(X>=R(2))) = R_max;
X(find(X<=R(1))) = R_min;
% quantization
adc_out = (round((X-R_min)./step))*step;
adc_out = adc_out + R_min;
adc_out(find(adc_out > R_max-offset)) = R_max-step;
deltaR = diff(R);
if deltaR > 0
  
  x_min = R(1) - deltaR/10;
  x_max = R(2) + deltaR/10;
else
  error('The upper limit must greater than the lower limit\n')
  
end
% subplot(2,1,1)
% plot((1:length(X)),[X; adc_out]);grid;title('ADC In/Out signals');
% axis([1 length(X) x_min x_max])
% legend('Inputs','Outputs','Location','northwest')
% subplot(2,1,2);grid; 
% plot(X-adc_out);grid;
% axis([1 length(X) x_min x_max])
% title('Qauntization Errors')
return