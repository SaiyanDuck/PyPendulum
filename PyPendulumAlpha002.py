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
import matplotlib.animation as animation


pyPendulumVer = "Alpha 0.02"
g = 9.8         #m/s^2     
l = 1           #m
tf1 = 10        #s
m = 5           #kg

D2phiDt = lambda t,S :[S[1],-g/l*np.sin(S[0])]

def penduloSimples(S0:list,tf,l):
    t_eval = np.linspace(0, tf, 1000) 
    sol = solve_ivp(D2phiDt,t_span=[0,tf],y0 = S0, t_eval = t_eval, method='RK45',rtol=1e-6, atol=1e-9)
    print(sol)
    tempo = sol.t
    theta_valores = sol.y[0]                    # Primeira coluna da solução é theta
    omega_valores = sol.y[1]                    # Segunda coluna é omega

    #Calculo das Energias
    kinEnergy = 1/2*m*(l*omega_valores)**2      
    PotEnergy = m*g*l*(1-np.cos(theta_valores))

    return np.array([tempo,theta_valores,omega_valores,kinEnergy,PotEnergy])

# ===========Plot do retrato de fase=============

def plotRetratoFasePenduloSimples(S:list = None, S0:list = None):
    theta_min, theta_max = -(2+1/2) * np.pi, (2+1/2) * np.pi              # limites inferiores e superiores do eixo horizontal (theta)
    omega_min, omega_max = -10, 10                                        # limites inferiores e superiores do eixo vertical   (omega)

    theta_grid = np.linspace(theta_min, theta_max, 30)                    # 30 pontos no eixo horizontal  (theta)
    omega_grid = np.linspace(omega_min, omega_max, 30)                    # 30 pontos no eixo vertical    (omega)
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
        ax.plot(Sol[1][1], Sol[2][1], 'go', markersize=12, label='Ponto Inicial')

    # Configurações finais do gráfico
    ax.set_title('Retrato de Fase com Trajetória Específica')
    ax.set_xlabel('Ângulo (θ) [rad]')
    ax.set_ylabel('Velocidade Angular (ω) [rad/s]')
    ax.grid(True)
    ax.legend()
    ax.set_xlim(theta_min, theta_max)
    ax.set_ylim(omega_min, omega_max)

    plt.show()




#======================Plot Da animação=============================
def PenduloSimplesAnimado(S,L):
    #Converter os ângulos (rad) para coordenadas cartesianas (x, y)
    theta_valores = S[1]
    tempo = S[0]
    x_valores = L * np.sin(theta_valores)
    y_valores = -L * np.cos(theta_valores)

    #Plot do grafico inicial
    fig, ax = plt.subplots(figsize=(8, 8))
    limite = L * 1.2
    ax.set_xlim(-limite, limite)
    ax.set_ylim(-limite, limite)
    ax.set_aspect('equal') # Garante que o aspecto seja 1:1
    ax.grid(True)
    ax.set_title('Animação do Pêndulo Simples')

    # Cria os elementos gráficos que serão animados.
    # Começamos com eles "vazios".
    # A vírgula depois de 'line' é um truque para desempacotar a lista de um elemento que o plot retorna.
    ponto_fixo, = ax.plot(0, 0, 'ko', markersize=5) # Ponto de pivô
    line, = ax.plot([], [], 'o-', lw=3, color='royalblue', markersize=10) # A haste e a massa
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes) # Texto para mostrar o tempo

    # Função de inicialização: desenha o fundo de cada quadro.
    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    # Função de animação: é chamada sequencialmente para cada quadro.
    def animate(i):
        # O argumento 'i' é o número do quadro, que usamos como índice nos nossos arrays de dados.
        x_atual = x_valores[i]
        y_atual = y_valores[i]

        # Atualiza a posição da linha do pêndulo. Ela vai de (0,0) até (x,y).
        line.set_data([0, x_atual], [0, y_atual])

        # Atualiza o texto do tempo
        time_text.set_text(f'Tempo: {tempo[i]:.1f}s')

        # Retorna uma tupla de todos os artistas que foram modificados.
        # Isso é necessário para a otimização 'blit=True'.
        return line, time_text

    # Cria a animação!
    # O FuncAnimation mantém uma referência à animação, então precisamos guardá-lo em uma variável.
    ani = animation.FuncAnimation(
        fig,
        animate,
        frames=len(tempo), # Número de quadros é o total de pontos de tempo
        init_func=init,
        interval=20,  # Atraso entre os quadros em milissegundos (1000/50 = 20 para 50 FPS)
        blit=True       # Otimização: redesenha apenas as partes que mudaram
    )

    # Exibe a animação
    plt.show()


    # --- 3. BLOCO DE EXECUÇÃO PRINCIPAL (INTERATIVO) ---

    if __name__ == "__main__":
        print("--- Simulador Interativo de Pêndulo Simples ---")
        print("Vamos configurar a simulação para baixas amplitudes (θ <= 45°).\n")

        # Coletando e validando o ângulo inicial
        while True:
            try:
                theta0_graus = float(input("Digite o ângulo inicial em graus (ex: 30): "))
                if abs(theta0_graus) <= 45:
                    break # Sai do loop se o ângulo for válido
                else:
                    print("Por favor, insira um ângulo menor ou igual a 45 graus.")
            except ValueError:
                print("Entrada inválida. Por favor, digite um número.")

        # Coletando outros parâmetros
        try:
            l_pendulo = float(input("Digite o comprimento do pêndulo em metros (ex: 2): "))
            t_final = float(input("Digite o tempo total de simulação em segundos (ex: 10): "))
        except ValueError:
            print("Entrada inválida. Usando valores padrão (L=2.0m, t=10s).")
            l_pendulo = 2.0
            t_final = 10.0

        # Convertendo ângulo para radianos e montando as condições iniciais
        theta0_rad = np.radians(theta0_graus)
        omega0 = 0.0
        s0 = [theta0_rad, omega0]

        # Rodando a simulação
        print("\nCalculando o movimento...")
        tempo, theta_valores = D2phiDt(s0, t_final, l_pendulo)
        print("Simulação concluída!")

        # Gerando e exibindo a animação
        print("Gerando a animação...")
        D2phiDt(tempo, theta_valores, l_pendulo)



teste1 = [np.pi/3,0]
teste2 = [3*np.pi/2,0]
teste3 = [np.pi,0]
teste4 = [-8,6.4]

teste11 = penduloSimples(teste1,tf1,l)

print(PenduloSimplesAnimado(teste1,l))


"""
i = 1
while i <= 4:

    testeplot = penduloSimples(eval(f"teste{i}"),tf1,l)

    #Plot da função theta(t)        (posição em relação ao tempo)
    plt.title('Teste %i para os $v_{ini}$ de $\Theta = %.2f$ e $\omega = %.2f$' %(i,eval(f'teste{i}[0]'),eval(f'teste{i}[1]')))
    plt.plot(testeplot[0],testeplot[1])
    plt.grid(True)
    plt.ylabel('Ângulo (θ) [rad]')
    plt.xlabel('Tempo (t) [seg]')
    plt.show()

    #Plot da função omega(theta)    (Espaço de fase)
    plt.clf()
    plt.title('Espaço de fase do teste %i para os $v_{ini}$ de $\theta = %.2f$ e $\omega = %.2f$' %(i,eval(f'teste{i}[0]'),eval(f'teste{i}[1]')))
    plt.plot(testeplot[1],testeplot[2])
    plt.grid(True)
    plt.ylabel('Velocidade Angular (ω) [rad/s]')
    plt.xlabel('Ângulo (θ) [rad]')
    plt.show()

    plotRetratoFasePenduloSimples(  [testeplot[1],testeplot[2]], eval(f"teste{i}"))

    #Plot da função De enregia
    plt.clf()               
    plt.title('Grafico da energia do teste %i para os $v_{ini}$ de $\Theta = %.2f$ e $\omega = %.2f$' %(i,eval(f'teste{i}[0]'),eval(f'teste{i}[1]')))
    plt.plot(testeplot[0],testeplot[3],color = 'blue')
    plt.plot(testeplot[0],testeplot[4],color = 'red')
    plt.plot(testeplot[0],testeplot[3]+testeplot[4],color ='yellow')
    plt.ylabel('Energia (Joule)')
    plt.xlabel('Ângulo (θ) [rad]')
    plt.legend(['K','U','Etot'])
    plt.grid(True)
    plt.show()
    i+=1
"""