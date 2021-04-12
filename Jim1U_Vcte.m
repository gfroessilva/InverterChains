function Jim1_U_tilde = Jim1U_Vcte(syspar,T,comm,tau,etai,ki,Xio,Xid,Oo,Oxi,ctrl_opt)
gim1 = comm(1);
gip1 = comm(2);
try
    lim1 = comm(3);
    lip1 = comm(4);
catch
end
alpha = syspar.alpha;

switch ctrl_opt
    case 1
        Jim1_U = [0 0 0;
            etai*alpha/tau 0 0;
            0 0 gim1/ki];
    case 2
        Jim1_U = [0 0 0;
            etai*alpha/tau 0 0;
            Xid*alpha/ki 0 gim1/ki];
    case 21
        Jim1_U = [0 0 0;
            etai*alpha/tau 0 0;
            Xid*alpha/ki 0 gim1/ki];        
    case 22
        Jim1_U = [0 0 0;
            0 0 0;
            0 0 gip1/ki];        
    case 23
        Jim1_U = [0 0 0;
            etai*alpha/tau 0 0;
            (lim1)/ki 0 gim1/ki];
    case 24
        Jim1_U = [0 0 0;
            etai*alpha/tau 0 0;
            2*(lim1*alpha)/ki 0 gim1/ki];
end

Jim1_U_tilde = T*Jim1_U/T;
end