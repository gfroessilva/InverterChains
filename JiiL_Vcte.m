function Jii_L_tilde = JiiL_Vcte(syspar,T,gip1,gim1,tau,etai,ki,Xio,Oo,Oxi)
alpha = syspar.alpha;
Jii_L = [0 1 0;
    -etai*(2*alpha)/tau -Oo/tau Oxi/tau;
    -0 -Xio/ki -(gim1 + gip1)/ki];

Jii_L_tilde = T*Jii_L/T;
end