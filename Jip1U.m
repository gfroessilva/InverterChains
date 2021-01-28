function Jip1_U_tilde = Jip1U(Gamma,syspar,comm,ctrlgain,ctrlxtra,Vbar)
tau = ctrlgain(1);
np = ctrlgain(2);
nq = ctrlgain(3);
kxi = ctrlgain(4);
kzeta = ctrlgain(5);
l = ctrlxtra(1);
m = ctrlxtra(2);
o = ctrlxtra(3);
x = ctrlxtra(4);
y = ctrlxtra(5);
gim1 = comm(1);
gip1 = comm(2);
him1 = comm(3);
hip1 = comm(4);
Bij = syspar.Bij;
Bii = syspar.Bii;

Jip1_U = [0                  0 0               0                 0;
          -np*Bij*Vbar^2/tau 0 0               0                 0;
          0                  0 np*Bij*Vbar/tau 0                 0;
          0                  0 0               gip1/kxi          0;
          0                  0 0               0        hip1/kzeta];
Jip1_U_tilde = Gamma*Jip1_U/Gamma; 
end