function decimal= mybin2real(binarystring)
% Converts an input string of 1s and 0s to a number.
Ind = strfind(binarystring, '.');
L = length(binarystring);
if isempty(Ind)
    Ind = L+1;
end
Num1 = binarystring(1:Ind-1);
LN1 = length(Num1);
Num2 = binarystring(Ind+1:end);
LN2 = length(Num1);
dec1=0;
for ii = 1 : LN1
    dec1 = dec1 + str2double(Num1(LN1-ii+1)) * 2^(ii-1);
end
dec2=0;
for ii = 1 : length(Num2) 
    dec2 = dec2 + str2double(Num2(ii)) * 2^-(ii);
end
decimal=dec1+dec2;