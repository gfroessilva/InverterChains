function Jii_L_tilde = JiiL(Gamma,syspar,comm,ctrlgain,ctrlxtra,Vbar)
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

Jii_L = [0               1      0                              0         0;
    -np/tau*2*Bij*Vbar^2 -1/tau 0                              y/tau     0;
    0                    0      -nq*(2*Bij-2*Bii)*Vbar/tau-o/tau 0   m/tau;
    0                    -x/kxi 0                      -(gip1+gim1)/kxi  0;
    0                    0      -l/kzeta             0 -(hip1+him1)/kzeta];
Jii_L_tilde = Gamma*Jii_L/Gamma;
end