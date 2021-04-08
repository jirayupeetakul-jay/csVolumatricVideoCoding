function x_sol=OMP(A_orig,b,S)

%%Solves the problem Ab=x s.t. |x|_0<S
%%Input: measurement matrix A, observation b and sparsity S
%%Output x obtained through orthogonal matching pursuit

for k=1:size(A_orig,2)
    
 A(:,k)=A_orig(:,k)/norm(A_orig(:,k),2);
end



iter=0;
residue=b;

x_sol=zeros(size(A_orig,2),1);
while iter<S
  %fprintf('iter=%d\n',iter);
  iter=iter+1;
  projection=A'*residue;
  [~,m]=max(abs(projection));

  n=m+iter-1;
  inculsion_set(:,iter)=A_orig(:,n);
  
  
  lambda_s=pinv(inculsion_set'*inculsion_set)*(inculsion_set'*b);
  
  
  
  max_set(iter)=n;
  
  x_sol(max_set)=lambda_s;
  
  residue=b-inculsion_set*lambda_s;
  
  A(:,m)=[];
  
end
end