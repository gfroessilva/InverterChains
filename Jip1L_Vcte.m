function Jip1_L_tilde = Jip1L_Vcte(syspar,T,gip1,gim1,tau,etai,ki,Xio,Oo,Oxi)
alpha = syspar.alpha;
Jip1_L = [0 0 0;
          0 0 0;
          0 0 gip1/ki];
Jip1_L_tilde = T*Jip1_L/T;
end