import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import curve_fit

N = 763.0
t = np.arange(14, dtype=float)
B = np.array([3,8,26,76,225,298,258,253,189,128,68,29,14,4], dtype=float)
C = np.array([0,0,0,0,9,17,105,162,176,166,150,85,47,20], dtype=float)
Iobs = B + C

I0 = Iobs[0]
S0 = N - I0
R0 = 0.0

def sir_rhs(y, t, beta, gamma):
    S, I, R = y
    return (-beta*S*I/N, beta*S*I/N-gamma*I, gamma*I)

def infected_curve(t_eval, beta, gamma):
    sol = odeint(sir_rhs, (S0,I0,R0), t_eval, args=(beta,gamma))
    return sol[:,1]

popt, pcov = curve_fit(infected_curve, t, Iobs,
                       p0=(1.5,0.3),
                       bounds=(0.0,np.inf), maxfev=20000)
beta_hat, gamma_hat = popt

print("beta =", beta_hat, "day^-1")
print("gamma =", gamma_hat, "day^-1")
print("R0 =", beta_hat/gamma_hat)
print("1/gamma =", 1/gamma_hat, "days")
print("RSS =", np.sum((infected_curve(t,beta_hat,gamma_hat)-Iobs)**2))

tfine=np.linspace(0,13,500)
plt.figure(figsize=(8.2,5.2))
plt.plot(tfine,infected_curve(tfine,beta_hat,gamma_hat),
         linewidth=2,label="SIR ajustado")
plt.plot(t,Iobs,"o",markersize=5,label=r"Dados: $B+C$")
plt.xlabel("Dias desde 22 de janeiro de 1978")
plt.ylabel(r"Número de indivíduos em $I(t)$")
plt.title("Ajuste SIR aos dados do surto de influenza de 1978")
plt.grid(alpha=.25); plt.legend(); plt.tight_layout()
plt.savefig("calibracao_SIR_surto_BMJ_White.png",dpi=180)
plt.show()
