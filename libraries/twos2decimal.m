function [decimal] = twos2decimal(x,bits)
if numel(x) > 1
decimal = arrayfun(@twos2decimal,x,ones(size(x))*bits);
return;
end
if isnan(x)
decimal = nan;
return;
end
if bitget(x,bits) == 1,
decimal = (bitxor(x, 2^bits-1) + 1) * -1;
else
decimal = x;
end