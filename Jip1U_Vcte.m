function Jip1_U_tilde = Jip1U_Vcte(syspar,T,gip1,gim1,tau,etai,ki,Xio,Oo,Oxi)
alpha = syspar.alpha;
Jip1_U = [0 0 0;  
          etai*alpha/tau 0 0; 
          0 0 gip1/ki];
Jip1_U_tilde = T*Jip1_U/T; 
end