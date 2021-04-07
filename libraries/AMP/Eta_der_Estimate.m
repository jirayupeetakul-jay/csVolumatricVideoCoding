function [ ed ] = Eta_der_Estimate( v,sigma,Eta )

p=length(v);
sigma_zz=1e-3;
zz=randn(p,1);
ed=sum( zz.*( Eta(v+zz*sigma_zz,sqrt(sigma^2+sigma_zz^2))-Eta(v,sigma) ))/sigma_zz/p;

end

