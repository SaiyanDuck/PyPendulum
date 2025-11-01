"""
Feito por: Thales serafim santore

Computação 1 - POLI UFRJ 25.2

23/09/2025

PyPendulum versão Alpha 0.01
"""

# ===================== ROADMAP Ver 1.0 =========================
r"""
+---------------------+
|  Início do Programa |
+---------------------+
           |
           v
+-----------------------------------+      
|         Exibir Menu Principal     |      +-------------------------+
| 1. Pêndulo Simples                |----->| MÓDULO PÊNDULO SIMPLES  |
| 2. Pêndulo com Resistência (TODO) |      +-------------------------+
| 3. Pêndulo Forçado (TODO)         |                 |
| 4. Pêndulo Duplo (TODO)           |                 |
| S. Sair                           |                 v
+-----------------------------------+      +--------------------------------+
           | (Se 'S')                      | Solicitar Parâmetros Iniciais  |
           v                               | (L, m, g, theta0, t_final)     |
+---------------------+                    +--------------------------------+
|   Fim do Programa   |                                 |
+---------------------+                                 v
                                             +--------------------------------+
                                             |   Loop de Simulação (passo Δt) |
                                             |  - Calcular nova vel. angular  |
                                             |  - Calcular nova pos. angular  |
                                             |  - Calcular Energias (Ek, Ep)  |
                                             |  - Armazenar resultados        |
                                             +--------------------------------+
                                                         |
                                                         v
                                             +--------------------------------+
                                             |       Gerar Saídas Gráficas    |
                                             |  - Plotar Animação do Pêndulo  |
                                             |  - Plotar Gráfico θ vs. t       |
                                             |  - Plotar Gráfico Energias vs. t|
                                             +--------------------------------+
                                                         |
                                                         v
                                             +--------------------------------+
                                             |  Retornar ao Menu Principal    |
                                             +--------------------------------+
                                                         |
                                                         ----------------------
                                                                  | (volta para o menu)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45,solve_ivp
#import matplotlib.animation as


pyPendulumVer = "Alpha 0.01"
g = 9.8         #m/s^2     
l = 1           #m
tf1 = 10        #s
m= 5            #kg

D2phiDt = lambda t,S :[S[1],-g/l*np.sin(S[0])]

def penduloSimples(S0:list,tf,l):
    t_eval = np.linspace(0, tf, 1000) 
    sol = solve_ivp(D2phiDt,t_span=[0,tf],y0 = S0, t_eval = t_eval, method='RK45',rtol=1e-6, atol=1e-9)
    print(sol)
    tempo = sol.t
    theta_valores = sol.y[0] # Primeira linha da solução é theta
    omega_valores = sol.y[1] # Segunda linha é omega
    kinEnergy = 1/2*m*(l*omega_valores)**2
    PotEnergy = m*g*l*(1-np.cos(theta_valores))

    return np.array([tempo,theta_valores,omega_valores,kinEnergy,PotEnergy])



#============ Plot do espaço de fase ==============
theta_min, theta_max = -(2+1/2) * np.pi, (2+1/2) * np.pi
omega_min, omega_max = -8, 8

theta_grid = np.linspace(theta_min, theta_max, 30) # 30 pontos no eixo theta
omega_grid = np.linspace(omega_min, omega_max, 30) # 30 pontos no eixo omega
THETA, OMEGA = np.meshgrid(theta_grid, omega_grid)


dTHETA_dt = np.zeros_like(THETA)
dOMEGA_dt = np.zeros_like(OMEGA)
for j in range(len(theta_grid)):
    for i in range(len(omega_grid)):
        # Pega o ponto (theta, omega) atual da grade
        ponto_atual = [THETA[i, j], OMEGA[i, j]]
        
        # Calcula o vetor de fluxo nesse ponto usando sua função EDO
        fluxo = D2phiDt(0, ponto_atual) # O tempo 't=0' não importa aqui
        
        # Guarda as componentes do vetor
        dTHETA_dt[i, j] = fluxo[0]
        dOMEGA_dt[i, j] = fluxo[1]

# ===========plot do retrato de fase=============

def plotRetratoFasePenduloSimples(S:list = None, S0:list = None):
    theta_min, theta_max = -(2+1/2) * np.pi, (2+1/2) * np.pi            #limites inferiores e superiores do eixo horizontal (theta)
    omega_min, omega_max = -10, 10                                        #limites inferiores e superiores do eixo vertical   (omega)

    theta_grid = np.linspace(theta_min, theta_max, 30)                  # 30 pontos no eixo horizontal  (theta)
    omega_grid = np.linspace(omega_min, omega_max, 30)                  # 30 pontos no eixo vertical    (omega)
    THETA, OMEGA = np.meshgrid(theta_grid, omega_grid)


    dTHETA_dt = np.zeros_like(THETA)
    dOMEGA_dt = np.zeros_like(OMEGA)
    for i in range(len(theta_grid)):
        for j in range(len(omega_grid)):
            ponto_atual = [THETA[j, i], OMEGA[j, i]]            # Pega o ponto (theta, omega) atual da grade
            fluxo = D2phiDt(0, ponto_atual)                     # Calcula o vetor de fluxo nesse ponto usando sua função EDO. O tempo 't=0' não importa aqui
            dTHETA_dt[j, i] = fluxo[0]
            dOMEGA_dt[j, i] = fluxo[1]
    
    
    #================= Plot do grafico======================
    fig, ax = plt.subplots(figsize=(12, 8))

    # CAMADA 1: Plotar o retrato de fase (pano de fundo)
    ax.streamplot(THETA, OMEGA, dTHETA_dt, dOMEGA_dt, density=1.2, linewidth=0.8, color='gray')    # Usando uma cor mais suave (cinza) e densidade menor para não poluir o gráfico

    if S is not None and S0 is not None:
        Sol = penduloSimples(S0,tf1,l)
        # CAMADA 2: Plotar a trajetória específica por cima
        ax.plot(Sol[1], Sol[2], color='blue', linewidth=2.5, label=f'Trajetória (θ₀={S0[1]:.2f} rad)')

        #CAMADA 3 (Bônus): Marcar o ponto inicial da trajetória
        #ax.plot(Sol[1], Sol[2], 'go', markersize=12, label='Ponto Inicial')

    # Configurações finais do gráfico
    ax.set_title('Retrato de Fase com Trajetória Específica')
    ax.set_xlabel('Ângulo (θ) [rad]')
    ax.set_ylabel('Velocidade Angular (ω) [rad/s]')
    ax.grid(True)
    ax.legend()
    ax.set_xlim(theta_min, theta_max)
    ax.set_ylim(omega_min, omega_max)

    plt.show()


teste1 = [np.pi/3,0]
teste2 = [3*np.pi/2,0]
teste3 = [np.pi,0]
teste4 = [-8,6.4]

i = 1
while i <= 4:

    testeplot = penduloSimples(eval(f"teste{i}"),tf1,l)

    plt.title('Teste %i para os $v_{ini}$ de $\Theta = %.2f$ e $\omega = %.2f$' %(i,eval(f'teste{i}[0]'),eval(f'teste{i}[1]')))
    plt.plot(testeplot[0],testeplot[1])
    plt.grid(True)
    plt.show()

    plt.clf()
    plt.title('Espaço de fase do teste %i para os $v_{ini}$ de $\Theta = %.2f$ e $\omega = %.2f$' %(i,eval(f'teste{i}[0]'),eval(f'teste{i}[1]')))
    plt.plot(testeplot[1],testeplot[2])
    plt.grid(True)
    plt.show()

    plotRetratoFasePenduloSimples(  [testeplot[1],testeplot[2]], eval(f"teste{i}"))

    plt.clf()
    plt.title('Grafico da energia do teste %i para os $v_{ini}$ de $\Theta = %.2f$ e $\omega = %.2f$' %(i,eval(f'teste{i}[0]'),eval(f'teste{i}[1]')))
    plt.plot(testeplot[0],testeplot[3],color = 'blue')
    plt.plot(testeplot[0],testeplot[4],color = 'red')
    plt.plot(testeplot[0],testeplot[3]+testeplot[4],color ='yellow')
    plt.legend(['K','U','Etot'])
    plt.grid(True)
    plt.show()
    i+=1

