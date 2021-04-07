
function [Store_Byte,Counts,Table]=Arith_Code(Data)
% Complelet Program for coding stream of symbols into stream of bits
% by uisng Airhtmetic coding Algoirthm.....
% input to the function :- array of "Data" 
%              by Mohammed M. Siddeq
%                  Date 2010/7/3
Table=0; New_Data=0; 
%Example for input :- 
%Data=[3 0 -1 129 9 -255 255 0 0 -3];  .......
% Note\ This program accept negrative, posative and zero data
%' Compute Table of Symbols ...'
Table(1)=Data(1);
S_AC=size(Data);
    for jAC=1:S_AC(2)
        S_2AC=size(Table);Flag=0;
        for kAC=1:S_2AC(2)
            if (Table(kAC)==Data(jAC))
                Flag=1;
            end;
        end;
        if (Flag==0) Table(S_2AC(2)+1)=Data(jAC); end;
    end;
  
 %' Compute the probability of the symbols...'   
 
 %%-----------------------------------------
  S_2AC=size(Table); 
  Counts(1:S_2AC(2))=0;
   for kAC=1:S_AC(2)
        iAC=1;
        while (Table(iAC)~=Data(kAC))
            iAC=iAC+1;
        end;
          New_Data(kAC)=iAC;
          Counts(iAC)=Counts(iAC)+1;
    end;
 %%%-------------------------------------
 code=0;
    %' Apply Arithmetic coding .... Now' 
    code = arithenco(New_Data,Counts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    %%% Converting bits to ASCII and save it into a file....
     %'converts bits to bytes... Now'
    D_Bits=char(code+48); % convert bits to characters '1' and '0' 
    
    N_Bits=8; %%% Choose number of bits for convertion 
    i_B=1; Size_AC=size(code);
    D_Bits(Size_AC(2)+1:Size_AC(2)+N_Bits)='0';
    i_B=1;k_Loc=1;
    while (i_B<=Size_AC(2))
     e8bits(1:N_Bits)=D_Bits(i_B:i_B+N_Bits-1);
     Store_Byte(k_Loc)=bin2dec(e8bits);
     k_Loc=k_Loc+1;
     i_B=i_B+N_Bits; 
    end;
       
    
       
       
end       
    % Exapmle :- %%%%%%%% Decode with Arithmetic Decoding %%%%%%%%
    %dseq = arithdeco(code,Counts,length(Data));
      % Note \ this instruction used for loading saved data from the a file...
      % Using \ S=load('c:\mm.dat','-mat');   
