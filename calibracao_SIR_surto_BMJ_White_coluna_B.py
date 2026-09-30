import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import curve_fit
N=763.0
t=np.arange(14,dtype=float)
B=np.array([3,8,26,76,225,298,258,253,189,128,68,29,14,4],dtype=float)
def rhs(y,t,beta,gamma):
    S,I,R=y
    return (-beta*S*I/N,beta*S*I/N-gamma*I,gamma*I)
def Icurve(tt,beta,gamma):
    return odeint(rhs,[760.,3.,0.],tt,args=(beta,gamma))[:,1]
popt,_=curve_fit(Icurve,t,B,p0=(1.5,.3),bounds=(0,np.inf),maxfev=20000)
beta,gamma=popt
print('beta =',beta,'day^-1')
print('gamma =',gamma,'day^-1')
print('R0 =',beta/gamma)
print('1/gamma =',1/gamma,'days')
print('RSS =',np.sum((Icurve(t,beta,gamma)-B)**2))
print('RMSE =',np.sqrt(np.mean((Icurve(t,beta,gamma)-B)**2)))
tfine=np.linspace(0,13,500)
plt.figure(figsize=(8.2,5.2))
plt.plot(tfine,Icurve(tfine,beta,gamma),linewidth=2,label='SIR ajustado')
plt.plot(t,B,'o',markersize=5,label='Dados: B (confined to bed)')
plt.xlabel('Dias desde 22 de janeiro de 1978')
plt.ylabel('Numero de individuos em I(t)')
plt.title('Ajuste SIR usando somente a coluna B')
plt.grid(alpha=.25); plt.legend(); plt.tight_layout()
plt.savefig('calibracao_SIR_surto_BMJ_White_coluna_B.png',dpi=180)
plt.show()
