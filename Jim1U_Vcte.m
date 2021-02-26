function Jim1_U_tilde = Jim1U_Vcte(syspar,T,gip1,gim1,tau,etai,ki,Xio,Oo,Oxi)
alpha = syspar.alpha;
Jim1_U = [0 0 0;  
          etai*alpha/tau 0 0; 
          0 0 gim1/ki];
Jim1_U_tilde = T*Jim1_U/T;
end