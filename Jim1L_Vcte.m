function Jim1_L_tilde = Jim1L_Vcte(syspar,T,gip1,gim1,tau,etai,ki,Xio,Oo,Oxi)
alpha = syspar.alpha;
Jim1_L = [0 0 0; 
          0 0 0;
          0 0 gim1/ki];
Jim1_L_tilde = T*Jim1_L/T;
end