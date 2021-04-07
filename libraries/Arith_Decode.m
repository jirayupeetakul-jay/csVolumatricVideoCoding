function [Data]=Arith_Decode(Store_Byte,Counts,Table)
%Compelet program for  Arithmetic Decoding : read  compressed file to be
% converted into    
%%%  input :- read infomation from compressed file:-
%%%     1- Store_Byte
%%%     2- Counts
%%%     3- Table
       
       Counts=double(Counts);   
      
      Store_Byte=double(Store_Byte);
      %' Compute data size'
      length_Data=0;
      Size_Counts=size(Counts);
      for iAD=1:Size_Counts(2)
      length_Data=length_Data+Counts(iAD);
      end;
      %'Read Data from the array...'
      Size_Data=size(Store_Byte);
      
      Loc=1; code=0;
      for iAD=1:Size_Data(2)
          e8bits=dec2bin(Store_Byte(iAD),8);
          code(Loc:Loc+7)=e8bits(1:8);
          Loc=Loc+8;
      end;
      code=code-48;
      %'Apply Arithmetic Decode..'
      New_Data = arithdeco(code,Counts,length_Data);
      
      Data=0;
      %'Return original data'  
      for iAD=1:length_Data
       Data(iAD)=Table(New_Data(iAD));
      end;