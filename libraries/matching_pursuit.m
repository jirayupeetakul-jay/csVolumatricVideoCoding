function x_sol=matching_pursuit(A,b,S)

%%Solves the problem Ab=x s.t. |x|_0<S
%%Input: measurement matrix A, observation b and sparsity S
%%Output x obtained through matching pursuit



iter=0;
residue=b;
exclusion_set=zeros(size(b,1),S);
x_sol=zeros(size(A,2),1);
while iter<S
    
  iter=iter+1;
  projection=A'*residue;
  [~,m]=max(abs(projection));
  x_sol(m)=projection(m);
  
  exclusion_set(:,iter)=A(:,m);
  
  residue=residue-A(:,m).*x_sol(m);
  
  A(A==exclusion_set(:,iter))=0;
  
end
end