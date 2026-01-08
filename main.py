"""
Feito por: Thales serafim santore
Computação 1 - POLI UFRJ 25.2
PyPendulum - Versão 0.91 (Alguns fixes no codigo Final)

Descrição:
    Simulador interativo de pêndulos (Simples, Forçado/Amortecido e Duplo).
    Utiliza integração numérica (Runge-Kutta 45) para resolver as equações de movimento.
    Oferece visualização dupla (Espaço Físico + Espaço de Fase) e análise de sensibilidade (Caos).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp
import sys
import os
import json
import datetime







# ================= CONFIGURAÇÕES GLOBAIS =================

G = 9.81  # Aceleração da gravidade (m/s^2)
ARQUIVO_CONFIG = "pypendulum_config.json"





# ================= 1. UTILITÁRIOS E GERENCIAMENTO =================


def carregar_estado():
    """
    Carrega as últimas configurações usadas pelo usuário de um arquivo JSON.
         
    Retorna:
        dict: Dicionário contendo parâmetros como 'L', 'massa', 'tf', etc.
              Retorna valores padrão se o arquivo não existir.
    """
    # Valores padrão atualizados para o caso caótico solicitado
    estado_padrao = {
        "L": 1.0, "massa": 1.0, "tf": 10.0,
        "theta0": 45.0, "omega0": 0.0,
        "amortecimento": 0.5, "amp_forca": 1.2, "freq_forca": 0.66,
        "theta1": 120.0, "omega1": 0.0, "theta2": 30.0, "omega2": 0.0, 
        "ultimo_acesso": ""
    }
    if os.path.exists(ARQUIVO_CONFIG):
        try:
            with open(ARQUIVO_CONFIG, 'r') as f:
                dados = json.load(f)
            # Garante compatibilidade mesclando chaves novas com antigas
            for k, v in estado_padrao.items():
                if k not in dados: 
                    dados[k] = v
            return dados
        except: 
            return estado_padrao
    return estado_padrao


def salvar_estado(estado):
    """
    Salva as configurações atuais no arquivo JSON ao sair do programa.
    
    Args:
        estado (dict): O dicionário com as configurações atuais da simulação.
        
    Retorna:
        None.
    """
    try:
        estado["ultimo_acesso"] = str(datetime.datetime.now())
        with open(ARQUIVO_CONFIG, 'w') as f: 
            json.dump(estado, f, indent=4)
    except: 
        pass


def salvar_dados_dat(cabecalho, dados_matriz, nome_arquivo="pendulo_dados.dat"):
    """
    Exporta os resultados numéricos para um arquivo de texto (.dat).
    
    Args:
        cabecalho (str): String contendo os nomes das colunas.
        dados_matriz (list): Lista de arrays (tempo, posições, energias, etc).
        nome_arquivo (str): Caminho para salvar o arquivo. Padrão: "pendulo_dados.dat".
        
    Retorna:
        None.
    """
    try:
        with open(nome_arquivo, 'w') as f:
            f.write(cabecalho + "\n")
            linhas, colunas = len(dados_matriz[0]), len(dados_matriz)
            for i in range(linhas):
                # Formata cada linha com 4 casas decimais
                f.write("\t".join([f"{dados_matriz[j][i]:.4f}" for j in range(colunas)]) + "\n")
        print(f">> Exportado com sucesso: {nome_arquivo}")
    except Exception as e: print(f"Erro na exportação: {e}")


def normalizar_angulo(rad):
    """
    Converte qualquer ângulo para o intervalo [-pi, pi].
    
    Args:
        rad (float): O ângulo de entrada em radianos.
        
    Returns:
        float: O ângulo equivalente normalizado entre -PI e +PI.

    OBS: Essencial para evitar que o gráfico de fase cresça indefinidamente (ex: 1000 graus).
    """
    return (rad + np.pi) % (2 * np.pi) - np.pi







# ================= 2. FÍSICA E CÁLCULO NUMÉRICO =================


#     ========= 2.1 EDO's ==========


def edo_generica_simples(t, S, L, b, F, wf):
    """
    Define as Equações Diferenciais Ordinárias (EDO) para o pêndulo simples.
    
    Args:
        t (float): Instante de tempo atual.
        S (list): Vetor de estado [theta, omega] (posição e velocidade angular).
        L (float): Comprimento do fio (m).
        b (float): Coeficiente de amortecimento.
        F (float): Amplitude da força externa de excitação.
        wf (float): Frequência angular da força externa.
        
    Retorna:
        list: Derivadas do estado [d_theta/dt, d_omega/dt].
    """
    theta, omega = S
    dtheta = omega
    domega = F*np.cos(wf*t) - b*omega - (G/L)*np.sin(theta)
    return [dtheta, domega]


def edo_pendulo_duplo(t, S, L1, L2, m1, m2):
    """
    Define as EDOs acopladas do Pêndulo Duplo (Lagrangiana).
    
    Args:
        t (float): Instante de tempo atual.
        S (list): Vetor de estado [theta1, omega1, theta2, omega2].
        L1, L2 (float): Comprimentos das hastes 1 e 2.
        m1, m2 (float): Massas dos pêndulos 1 e 2.
        
    Returns:
        list: Derivadas [omega1, acel_angular1, omega2, acel_angular2].
    """
    theta1, omega1, theta2, omega2 = S
    delta = theta1 - theta2
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    den2 = (L2 / L1) * den1
    
    domega1 = (m2 * L1 * omega1**2 * np.sin(delta) * np.cos(delta) +
           m2 * G * np.sin(theta2) * np.cos(delta) +
           m2 * L2 * omega2**2 * np.sin(delta) -
           (m1 + m2) * G * np.sin(theta1)) / den1
           
    domega2 = (-m2 * L2 * omega2**2 * np.sin(delta) * np.cos(delta) +
           (m1 + m2) * G * np.sin(theta1) * np.cos(delta) -
           (m1 + m2) * L1 * omega1**2 * np.sin(delta) -
           (m1 + m2) * G * np.sin(theta2)) / den2
    return [omega1, domega1, omega2, domega2]




#     ========= 2.2 Energia ==========


def calcular_energias(tipo, sol_y, params):
    """
    Calcula a Energia Cinética (K), Potencial (U) e Mecânica (E) para cada instante.
    
    Args:
        tipo (int): 1 ou 2 para Simples/Forçado, 3 para Duplo.
        sol_y (numpy.ndarray): Matriz com a solução das EDOs (linhas=variáveis, colunas=tempo).
        params (list): Lista dos parâmetros físicos usados na simulação (L, m, etc).
        
    Returns:
        tuple: Três arrays numpy (K, U, E_total).
    """
    # ...
    if tipo in [1, 2]: # Simples
        theta, omega = sol_y[0], sol_y[1]
        L, m = params[2], params[3]
        K = 0.5 * m * (L * omega)**2
        U = m * G * L * (1 - np.cos(theta))
        return K, U, K+U
        
    elif tipo == 3: # Duplo
        theta1, omega1, theta2, omega2 = sol_y
        L, m = params[4], params[5]
        # Velocidade escalar ao quadrado (Lei dos cossenos para massa 2)
        v1_sq = (L*omega1)**2
        v2_sq = (L*omega1)**2 + (L*omega2)**2 + 2*L*L*omega1*omega2*np.cos(theta1-theta2)
        
        K = 0.5 * m * v1_sq + 0.5 * m * v2_sq
        y1 = -L * np.cos(theta1)
        y2 = y1 - L * np.cos(theta2)
        U = m*G*y1 + m*G*y2 - (-(m+m)*G*(L+L)) # Ajuste de referência zero
        return K, U, K+U


#     ====== 2.3 Resolver EDO's =======


def resolver_sistema(tipo, params, tf):
    """
    Configura e executa a integração numérica (Runge-Kutta 45).

    Args:
        tipo (int)      : 1 (Simples), 2 (Forçado/Amortecido) ou 3 (Duplo).
        params (list)   : Lista contendo condições iniciais e constantes físicas na ordem correta.
        tf (float)      : Tempo final da simulação em segundos.

    Mapeamento da lista 'params':
        Tipo 1 (Simples): [theta0, omega0, L, m]
        Tipo 2 (Forçado): [theta0, omega0, L, m, b, A, wf]
        Tipo 3 (Duplo)  : [theta1, omega1, theta2, omega2, L, m]

    Retorna (Lista de Arrays):
        Simples/Forçado: [tempo, theta, omega, K, U, E_total]
        Duplo          : [tempo, theta1, omega1, theta2, omega2, K, U, E_total]
    """
    fps = 30
    t_eval = np.linspace(0, tf, int(round(tf * fps)))
    
    if tipo == 1:       #Pendulo simples
        args = (params[2], 0, 0, 0)
        y0 = params[0:2]
    elif tipo == 2:     #pendulo forçado amortecido
        args = (params[2], params[4], params[5], params[6]) 
        y0 = params[0:2]
    elif tipo == 3:     #Pendulo duplo
        args = (params[4], params[4], params[5], params[5])
        y0 = params[0:4]
        sol = solve_ivp(edo_pendulo_duplo, [0, tf], y0, args=args, t_eval=t_eval, method='RK45', rtol=1e-8)
        K, U, E = calcular_energias(3, sol.y, params)
        return [sol.t, *sol.y, K, U, E]

    sol = solve_ivp(edo_generica_simples, [0, tf], y0, args=args, t_eval=t_eval, method='RK45')
    K, U, E = calcular_energias(tipo, sol.y, params)
    return [sol.t, sol.y[0], sol.y[1], K, U, E]             






# ================= 3. ANIMAÇÃO (DUAL VIEW) =================

# ===== 3.1 Motor de animação =====

def animar_dual(tipo, dados, config):
    """
    Configura a figura e inicia o loop da animação principal (Física + Fase).
    
    Args:
        tipo (int): Tipo de simulação (1, 2 ou 3).
        dados (list): Output da função 'resolver_sistema'.
        config (dict): Configurações globais (para ler comprimento L, etc).
        
    Returns:
        None.
    """
    L = config["L"]
    t_arr = dados[0]
    dt = int((t_arr[1]-t_arr[0])*1000) # Intervalo em ms
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    
    # --- CONFIG AX1: FÍSICA ---
    lim = 2.2*L if tipo==3 else 1.2*L
    ax1.set_xlim(-lim, lim); ax1.set_ylim(-lim, lim); ax1.set_aspect('equal')
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.set_title("Simulação Física (Real)", fontweight='bold')
    ax1.set_xlabel("x (m)"); ax1.set_ylabel("y (m)")
    ax1.plot(0, 0, 'k+', ms=10, markeredgewidth=2) # Pivô
    txt = ax1.text(0.05, 0.95, '', transform=ax1.transAxes, va='top', bbox=dict(facecolor='white', alpha=0.9))

    # --- CONFIG AX2: ESPAÇO DE FASE ---
    ax2.set_title("Espaço de Fase (θ vs ω)", fontweight='bold')
    ax2.set_xlabel("Ângulo θ (rad)"); ax2.set_ylabel("Velocidade ω (rad/s)")
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.set_xlim(-np.pi - 0.2, np.pi + 0.2) # Eixo X travado em [-PI, +PI]
    
    # Inicialização dos objetos gráficos (linhas vazias)
    if tipo == 3: # Pêndulo Duplo
        theta1, omega1, theta2, omega2 = dados[1], dados[2], dados[3], dados[4]
        x1, y1 = L*np.sin(theta1), -L*np.cos(theta1)
        x2, y2 = x1 + L*np.sin(theta2), y1 - L*np.cos(theta2)
        
        w_max = max(np.max(np.abs(omega1)), np.max(np.abs(omega2))) + 2
        ax2.set_ylim(-w_max, w_max)

        l, = ax1.plot([], [], 'k-', lw=2); pts, = ax1.plot([], [], 'o', ms=8, color='firebrick')      
        tr, = ax1.plot([], [], 'b-', lw=1, alpha=0.4)                                                 
        ph_tr1, = ax2.plot([], [], 'b-', lw=0.5, alpha=0.6, label='Haste 1')
        ph_pt1, = ax2.plot([], [], 'bo', ms=4)                                                  
        ph_tr2, = ax2.plot([], [], 'r-', lw=0.5, alpha=0.6, label='Haste 2')
        ph_pt2, = ax2.plot([], [], 'ro', ms=4)
        ax2.legend(loc='upper right')
        
        #As list vazias são buffers de memoria da trajetoria do pendulo para desenhar, ele para de salvar quando chega a 100 frames
        ani = animation.FuncAnimation(fig, update_dual_duplo, frames=len(t_arr), 
                                      fargs=(l, pts, tr, txt, ph_tr1, ph_pt1, ph_tr2, ph_pt2, x1, y1, x2, y2, theta1, omega1, theta2, omega2, t_arr, [], [], [], [], [], []), 
                                      interval=dt, blit=True)       
    else: # Simples/Forçado
        th, w = dados[1], dados[2]
        x, y = L*np.sin(th), -L*np.cos(th)         
        w_max = np.max(np.abs(w)) + 2
        ax2.set_ylim(-w_max, w_max)
        
        l, = ax1.plot([], [], '-', lw=2, color='#333')
        pt, = ax1.plot([], [], 'o', ms=15, color='navy')
        
        # OBS: Retrato de fase de fundo (streamplot) só é válido para sistemas autônomos (sem dependência explícita do tempo)
        is_autonomous = (tipo == 1) or (tipo == 2 and config["amp_forca"] == 0)
        
        if is_autonomous:
            Y, X = np.mgrid[-w_max:w_max:20j, -np.pi:np.pi:20j]
            if tipo == 1: U, V = Y, -(G/L)*np.sin(X)
            else: b = config["amortecimento"]; U, V = Y, -b*Y - (G/L)*np.sin(X)
            
            ax2.streamplot(X, Y, U, V, color='silver', density=0.8, linewidth=0.6, arrowsize=1)
            
        ph_l, = ax2.plot([], [], 'r-', lw=2, alpha=0.8); ph_p, = ax2.plot([], [], 'go', ms=8)
        
        ani = animation.FuncAnimation(fig, update_dual_simples, frames=len(t_arr), 
                                      fargs=(l, pt, txt, ph_l, ph_p, x, y, th, w, t_arr, [], []), 
                                      interval=dt, blit=True)

    plt.tight_layout()
    plt.show()
    
    # Opção de salvar GIF (requer pillow)
    if input("Salvar animação? (s/n): ").lower() == 's':
        try: 
            nome_arq = f"anim_{datetime.datetime.now().strftime('%H%M')}.gif"
            print(">> Gerando GIF (Aguarde,isso pode demorar)")
            ani.save(nome_arq, writer='pillow', fps=30)
            print(f">> Salvo: {nome_arq}")
        except Exception as e: print(f"Erro ao salvar: {e}")


def animar_sensibilidade(todos_dados, L):
    """
    Gerencia a animação de múltiplas trajetórias para visualizar o caos.
    
    Args:
        todos_dados (list): Lista de resultados de múltiplas simulações com CI próximas.
        L (float): Comprimento da haste para escala do gráfico.
        
    Returns:
        None.
    """
    tempo = todos_dados[0][0]
    dt = int((tempo[1]-tempo[0])*1000)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    lim = 2.2 * L
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect('equal'); ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_title("Efeito Borboleta: Divergência de Trajetórias", fontweight='bold')
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    
    ax.plot(0, 0, 'k+', ms=10)
    time_txt = ax.text(0.05, 0.95, '', transform=ax.transAxes, bbox=dict(fc='white', alpha=0.8))
    
    cores = ['blue', 'cyan', 'green', 'orange', 'red']
    linhas = []; massas = []; rastros = []
    h_x = [[] for _ in range(5)]; h_y = [[] for _ in range(5)]
    
    for j in range(5):
        l, = ax.plot([], [], '-', lw=2, color=cores[j], alpha=0.3)
        p, = ax.plot([], [], 'o', ms=6, color=cores[j], alpha=0.8)
        r, = ax.plot([], [], '-', lw=1, color=cores[j], alpha=0.6, label=f'Sim {j+1}')
        linhas.append(l); massas.append(p); rastros.append(r)
    
    ax.legend(loc='lower right', fontsize='small')
    ani = animation.FuncAnimation(fig, update_sensibilidade, frames=len(tempo), 
                                  fargs=(linhas, massas, rastros, time_txt, todos_dados, tempo, L, h_x, h_y), 
                                  interval=dt, blit=True)
    print(">> Exibindo divergência (diferença inicial: 0.005 rad)...")
    plt.show()
    
    # Opção de salvar GIF (requer pillow)
    if input("Salvar animação? (s/n): ").lower() == 's':
        try: 
            nome_arq = f"anim_{datetime.datetime.now().strftime('%H%M')}.gif"
            print(">> Gerando GIF (Aguarde,isso pode demorar)")
            ani.save(nome_arq, writer='pillow', fps=30)
            print(f">> Salvo: {nome_arq}")
        except Exception as e: print(f"Erro ao salvar: {e}")


# ===== 3.2 Update de Frames ======


def update_dual_simples(i, line, point, txt, phase_line, phase_pt, x, y, th, w, t, hist_theta, hist_w):
    """
    Atualiza frame a frame a animação do Pêndulo Simples/Forçado.
    
    Args:
        i (int): Índice do frame atual.
        line, point (matplotlib objects): Elementos gráficos do pêndulo físico.
        txt (matplotlib.text): Texto informativo do frame.
        phase_line, phase_pt (matplotlib objects): Elementos do gráfico de fase.
        x, y, th, w, t (numpy.array): Arrays com todos os dados pré-calculados.
        hist_theta, hist_w (list): Buffers para desenhar o rastro no espaço de fase.
        
    Returns:
        tuple: Objetos gráficos atualizados para o blitting.
    """
    # ...
    if i == 0: hist_theta.clear(); hist_w.clear()
    
    # Atualiza Física (Esquerda)
    line.set_data([0, x[i]], [0, y[i]])
    point.set_data([x[i]], [y[i]])
    
    # Atualiza Fase (Direita) - Com corte de linha para 'wrap around'
    th_norm = normalizar_angulo(th[i])
    if len(hist_theta) > 0 and abs(th_norm - hist_theta[-1]) > 3.0:
         hist_theta.append(np.nan); hist_w.append(np.nan)
    
    hist_theta.append(th_norm); hist_w.append(w[i])
    phase_line.set_data(hist_theta, hist_w)
    phase_pt.set_data([th_norm], [w[i]])
    
    txt.set_text(f't: {t[i]:.1f}s\nθ: {np.degrees(th_norm):.1f}°\nω: {w[i]:.2f}')
    return line, point, txt, phase_line, phase_pt


def update_dual_duplo(i, l1, p1, tr1, txt, ph_tr1, ph_pt1, ph_tr2, ph_pt2, x1, y1, x2, y2, theta1, omega1, theta2, omega2, t, h_x, h_y, h_theta1, h_omega1, h_theta2, h_omega2):
    """
    Atualiza frame a frame a animação do Pêndulo Duplo.
    
    Args:
        i (int): Índice do frame atual.
        l1, p1, tr1 (matplotlib objects): Linhas, massas e rastro da simulação física.
        ph_tr1, ph_pt1... (matplotlib objects): Rastros e pontos do espaço de fase.
        x1, y1, x2, y2... (numpy.array): Dados de posição pré-calculados.
        h_x, h_y... (list): Históricos (buffers) de posição para desenhar rastros.
        
    Returns:
        tuple: Objetos gráficos atualizados.
    """
    # ...

    if i == 0:
        h_x.clear(); h_y.clear(); h_theta1.clear(); h_omega1.clear(); h_theta2.clear(); h_omega2.clear()
    
    # Física: Hastes e Rastro
    l1.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
    p1.set_data([x1[i], x2[i]], [y1[i], y2[i]])
    
    h_x.append(x2[i]); h_y.append(y2[i])
    if len(h_x)>150: h_x.pop(0); h_y.pop(0) # Rastro limitado para não poluir
    tr1.set_data(h_x, h_y)
    
    # Fase: Normalização e Plotagem
    theta1_n = normalizar_angulo(theta1[i])
    theta2_n = normalizar_angulo(theta2[i])
    
    if len(h_theta1) > 0 and abs(theta1_n - h_theta1[-1]) > 3.0: h_theta1.append(np.nan); h_omega1.append(np.nan)
    if len(h_theta2) > 0 and abs(theta2_n - h_theta2[-1]) > 3.0: h_theta2.append(np.nan); h_omega2.append(np.nan)

    h_theta1.append(theta1_n); h_omega1.append(omega1[i])
    h_theta2.append(theta2_n); h_omega2.append(omega2[i])
    
    if len(h_theta1) > 400: h_theta1.pop(0); h_omega1.pop(0); h_theta2.pop(0); h_omega2.pop(0)
        
    ph_tr1.set_data(h_theta1, h_omega1); ph_pt1.set_data([theta1_n], [omega1[i]])
    ph_tr2.set_data(h_theta2, h_omega2); ph_pt2.set_data([theta2_n], [omega2[i]])
    
    txt.set_text(f't:{t[i]:.1f}s\nθ1:{np.degrees(theta1_n):.0f}°\nθ2:{np.degrees(theta2_n):.0f}°')
    return l1, p1, tr1, txt, ph_tr1, ph_pt1, ph_tr2, ph_pt2


def update_sensibilidade(i, linhas, massas, rastros, textos, todos_dados, tempo, L, h_x, h_y):
    """Atualiza 5 pêndulos simultaneamente."""
    if i == 0:
        for k in range(5): h_x[k].clear(); h_y[k].clear()

    for idx, dados in enumerate(todos_dados):
        theta1, theta2 = dados[1], dados[3]
        x1 = L * np.sin(theta1[i]); y1 = -L * np.cos(theta1[i])
        x2 = x1 + L * np.sin(theta2[i]); y2 = y1 - L * np.cos(theta2[i])
        
        linhas[idx].set_data([0, x1, x2], [0, y1, y2])
        massas[idx].set_data([x1, x2], [y1, y2])
        h_x[idx].append(x2); h_y[idx].append(y2)
        if len(h_x[idx]) > 200: h_x[idx].pop(0); h_y[idx].pop(0)
        rastros[idx].set_data(h_x[idx], h_y[idx])
        
    textos.set_text(f'Tempo: {tempo[i]:.1f}s')
    return linhas + massas + rastros + [textos]






# ================= 4. MENU E INTERFACE =================

def print_logo():
    logo = """
    ###########################
    #  PyPendulum Versão 0.91 #
    ###########################
           ___________
                |
                .
               /|
              / |
           L /  | g
            /   | 
           /    v
          O 

    Feito por: Thales Serafim
    UFRJ - Universidade Federal do Rio de Janeiro
    Escola Politecnica
    """
    print(logo)

def obter_float(msg, padrao):
    """
    Helper para input seguro de números decimais no console.
    
    Args:
        msg (str): Texto a ser exibido para o usuário.
        padrao (float): Valor default caso o usuário dê 'Enter' vazio.
        
    Returns:
        float: O valor inserido ou o padrão.
    """
    # ...
    val = input(f"{msg} [{padrao}]: ").strip()
    return float(val) if val else padrao

def menu_casos_famosos(cfg):
    """
    Carrega predefinições educativas (presets) para demonstração.
    
    Args:
        cfg (dict): Dicionário de configurações atual para ser atualizado.
        
    Returns:
        tuple: (tipo [int], parametros [list]). Retorna (None, None) se cancelado.
    """
    # ...
    print("\n--- GALERIA DE CASOS FAMOSOS ---")
    print("1. Equilíbrio Instável (Simples)")
    print("2. Espiral da Morte (Amortecido)")
    print("3. Sobreamortecido (Viscoso)")
    print("--- PÊNDULO DUPLO ---")
    print("4. Modo Normal Simétrico")
    print("5. Modo Normal Anti-Simétrico")
    print("6. Batimento (Energia)")
    print("7. Queda Horizontal")
    
    escolha = input(">>> ")
    cfg["L"] = 1.0; cfg["massa"] = 1.0; cfg["tf"] = 20.0
    
    if escolha == '1': return 1, [np.radians(179.9), 0, 1.0, 1.0]
    elif escolha == '2': 
        cfg["amortecimento"]=0.5; cfg["amp_forca"]=0.0; cfg["freq_forca"]=0.0
        return 2, [np.radians(160), 0, 1.0, 1.0, 0.5, 0.0, 0.0]
    elif escolha == '3': 
        cfg["amortecimento"]=5.0; cfg["amp_forca"]=0.0; cfg["freq_forca"]=0.0
        return 2, [np.radians(160), 10.0, 1.0, 1.0, 5.0, 0.0, 0.0]
    elif escolha == '4': return 3, [np.radians(10), 0, np.radians(10), 0, 1.0, 1.0]
    elif escolha == '5': return 3, [np.radians(15), 0, np.radians(-15), 0, 1.0, 1.0]
    elif escolha == '6': cfg["tf"]=40; return 3, [np.radians(20), 0, 0, 0, 1.0, 1.0]
    elif escolha == '7': return 3, [np.radians(90), 0, np.radians(90), 0, 1.0, 1.0]
    else: return None, None


def main():
    cfg = carregar_estado()
    print_logo()
    
    while True:
        print("Menu principal do PyPendulum: Digite qual caso você quer simular.")
        print("\n1. Simples | 2. Forçado | 3. Duplo | 4. CASOS FAMOSOS | 5. Sair")
        op = input(">>> ")
        if op == '5': salvar_estado(cfg); sys.exit()
        
        dados = []; tipo = 0
        
        # Opção 4: Carrega predefinições
        if op == '4':
            tipo, params = menu_casos_famosos(cfg)
            if tipo is None: continue
            print(">> Calculando predefinição...")
            dados = resolver_sistema(tipo, params, cfg["tf"])
            
        # Opções Manuais
        elif op in ['1', '2', '3']:
            cfg["L"] = obter_float("L", cfg.get("L", 1.0))
            cfg["massa"] = obter_float("m", cfg.get("massa", 1.0))
            cfg["tf"] = obter_float("Tempo", cfg.get("tf", 40.0)) # Default 40s para caos
            
            if op == '1':
                p = [np.radians(obter_float("Th0",45)), obter_float("W0",0), cfg["L"], cfg["massa"]]
                print(">> Calculando...")
                dados = resolver_sistema(1, p, cfg["tf"]); tipo = 1
                
            elif op == '2':
                p = [np.radians(obter_float("Th0",45)), obter_float("W0",0), cfg["L"], cfg["massa"], 
                     obter_float("b",0.5), obter_float("F",1.2), obter_float("w_f",0.66)]
                print(">> Calculando...")
                dados = resolver_sistema(2, p, cfg["tf"]); tipo = 2
                
            elif op == '3':
                # Inputs específicos para o pêndulo duplo
                theta1 = obter_float("Th1",120); omega1 = obter_float("W1",0)
                theta2 = obter_float("Th2",30); omega2 = obter_float("W2",0)
                
                # Verifica se o usuário quer análise de caos
                if input("Verificar sensibilidade (Borboleta)? (s/n): ").lower() == 's':
                    print(">> Gerando 5 simulações (pode demorar)...")
                    lst = []
                    pertubacao = 0.005 # Aprox 0.3 graus de diferença
                    for i in range(5):
                        # Cria pequena variação na condição inicial
                        p_mod = [np.radians(theta1+i*0.28), omega1, np.radians(theta2), omega2, cfg["L"], cfg["massa"]]
                        lst.append(resolver_sistema(3, p_mod, cfg["tf"]))
                    animar_sensibilidade(lst, cfg["L"])
                    continue # Pula a animação normal e volta pro menu
                
                print(">> Calculando único...")
                dados = resolver_sistema(3, [np.radians(theta1), omega1, np.radians(theta2), omega2, cfg["L"], cfg["massa"]], cfg["tf"])
                tipo = 3
        else: continue

        # Exibe a animação principal (Dual View)
        animar_dual(tipo, dados, cfg)
        
        # Sub-menu de pós-processamento
        while True:
            sub = input("[A] Plotar grafico de energia | [B] Exportar dados (.dat) | [C] Voltar: ").lower()
            if sub == 'a':
                t = dados[0]
                if len(dados)==6: K,U,E = dados[3], dados[4], dados[5]
                else: K,U,E = dados[5], dados[6], dados[7]
                plt.figure(figsize=(10,5))
                plt.plot(t,K,'b--',label='Cinética'); plt.plot(t,U,'g:',label='Potencial'); plt.plot(t,E,'r-',label='Total')
                plt.legend(); plt.grid(True); plt.title("Análise de Energia"); plt.show()
            elif sub == 'b': salvar_dados_dat("t, theta, omega, K   ,U   ,E   ", dados)
            elif sub == 'c': break

if __name__ == "__main__":      #Caso utilize import PyPendulum
    main()