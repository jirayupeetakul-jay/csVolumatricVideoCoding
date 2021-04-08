function varargout = Huff(xC)
% Encoder: [data,bits] = Huff(xC);
% Decoder: xC = Huff(data);
if isstruct(xC{1})
    Encoder = 0; Decoder = 1;
else
    Encoder = 1; Decoder = 0;
end

if Encoder
   bits = 0;
   len = length(xC);
   data = cell(len,1);
   for ii = 1:len
       data{ii} = mat2huff(xC{ii});
       bits = bits + bytes(data{ii})*8;
   end
   varargout{1} = data;
   varargout{2} = bits;
end

if Decoder
   data = xC;
   len = length(data); 
   for ii = 1:len
       xC{ii} = huff2mat(data{ii});
   end
   varargout{1} = xC;
end

end
       
       
       