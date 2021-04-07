function [outSTR,Sizevec,evalSTR] = run_length_com(vec)
% function [outSTR,Sizevec,evalSTR] = compress(vec)
%
% Run length compression of a vector or matrix by creating
% a string of Matlab evaluable code! Could thus also be 
% called 'data2mcode' or something like that.
%
% This function replaces streams of same numbers with 
% Matlab code and returns a string which can be evaluated!
% Really funny way to save memory or transform i.e. pictures
% to Matlab executable code. I used it to make a Simulink
% block icon without the need of a file on the disk.
% (See the Simulink model in the ddfread.m on MatlabCentral)
%
% vec      a vector, matrix of any dimension
% outSTR   the string which describes the data in vector format
% Sizevec  the size of the input data 
% evalSTR  the string which can be evaluated 
%
% Example:
%          x = ones(5)
%          x([4 5 7]) = 0.5
%          [cS,Sv,eS] = compress(x)
%          y=eval(eS)
%          y-x
%                                  (c) 2002, B. Kaeferstein
%                                    berthold@kaefersten.de
Sizevec = size(vec);
vec = vec(:)';
diffvec = diff(vec);
[res,idx] = find(diffvec~=0);
idx=[0,idx,size(vec,2)];
compressSTR = '[';
outSTR      = '[';
for repDATA = 1 : size(idx,2)-1
    
    idxBegin = idx (repDATA)  +1;
    idxEnd   = idx (repDATA+1);  
    
    repBegin = vec (idxBegin);
    repEnd   = vec (idxEnd);
         
    test = [repBegin:repEnd];
    if idxBegin == idxEnd %compress only if values are different!
        compressSTR = [compressSTR,num2str(repEnd)];
        outSTR      = [outSTR,num2str(repEnd)];
    else    
        compressSTR = [compressSTR,'[',num2str(1),':',num2str(idxEnd-idxBegin+1),']*0+',num2str(repEnd)];
        outSTR      = [outSTR,'[',num2str(1),':',num2str(idxEnd-idxBegin+1),']*0+',num2str(repEnd)];    
    end
    
    if repDATA  ~= size(idx,2)-1
        compressSTR = [compressSTR,','];
        outSTR      = [outSTR,','];
    end
    
    if mod(repDATA,15)==0
        compressSTR= [compressSTR,''',...',char(10),'''']; % start sometimes a new line
    end
        
end
compressSTR = [compressSTR,']'];
     outSTR = [outSTR,']'];
evalSTR     =['reshape(eval([''',compressSTR,''']),[',num2str(Sizevec),'])'];