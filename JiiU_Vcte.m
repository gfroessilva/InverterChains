function Jii_U_tilde = JiiU_Vcte(syspar,T,gip1,gim1,tau,etai,ki,Xio,Oo,Oxi)
alpha = syspar.alpha;
Jii_U = [0 1 0;
    0 -Oo/tau Oxi/tau;
    0 -Xio/ki -(gim1 + gip1)/ki];

Jii_U_tilde = T*Jii_U/T;
end