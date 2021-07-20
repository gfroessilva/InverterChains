function Jim1_L_tilde = Jim1L_Vcte(syspar,T,comm,tau,etai,ki,Xio,Xid,Oo,Oxi,ctrl_opt)
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
        Jim1_L = [0 0 0;
            0 0 0;
            0 0 gim1/ki];
    case 2
        Jim1_L = [0 0 0;
            0 0 0;
            0 0 gim1/ki];        
    case 21
        Jim1_L = [0 0 0;
            0 0 0;
            0 0 gim1/ki];        
    case 22
        Jim1_L = [0 0 0;
            0 0 0;
            0 0 gip1/ki];        
    case 23
        Jim1_L = [0 0 0;
            0 0 0;
            lim1/ki 0 gim1/ki];     
    case 24
        Jim1_L = [0 0 0;
            0 0 0;
            0 0 gip1/ki];        
end

Jim1_L_tilde = T*Jim1_L/T;
end