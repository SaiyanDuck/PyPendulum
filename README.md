PyPendulum v002
Autor: Thales Serafim Santore Disciplina: Computação 1 - POLI UFRJ (2025.2) Data: 23/09/2025

-=[Objetivo do Projeto]=-
O intuito desse programa é mostrar de maneira interativa com o usuário as propriedades físicas de um pêndulo, indo dos casos mais simples (pêndulo simples) até os mais avançados (pêndulo duplo e pêndulo amortecido forçado).

O programa realiza uma simulação animada do pêndulo com os inputs do usuário e plota dados importantes que detalham o movimento (energias, retratos de fase, etc.).


-=[Requisitos]=-
Este projeto utiliza as seguintes bibliotecas:
* Python 3.x
* Matplotlib 
* NumPy 
* SciPy


-=[Roadmap do Projeto]=-
Grupo 1: Casos Básicos
CASO 1.1: PÊNDULO SIMPLES

Status: (70% COMPLETO)

O caso mais simples, utiliza o algoritmo de Runge-Kutta para plotar com precisão o movimento do pêndulo. É também plotado o retrato de fase e as energias.

A fazer: Realizar alguns ajustes que vão ser importantes para o programa base.

CASO 1.2: PÊNDULO CÔNICO (Caso extra!)

Status: (ESPECULATIVO)

Esse é um algoritmo que não sei se vou levar para a versão final, pois é um caso tridimensional, e eu queria focar nos problemas bidimensionais.

Grupo 2: Casos Intermediários
CASO 2.1: PÊNDULO AMORTECIDO

Status: (0% COMPLETO)

Um caso famoso da física básica. Quero usar isso para detalhar os regimes de amortecimento (subcrítico, crítico e supercrítico).

(Sugestão): Incluir a discussão sobre modelos de arrasto (linear vs. quadrático) e a analogia com circuitos RLC da Engenharia Elétrica.

CASO 2.2: PÊNDULO FORÇADO

Status: (0% COMPLETO)

Trabalhar com a simulação do pêndulo com atuação de forças externas.

CASO 2.3: PÊNDULO AMORTECIDO FORÇADO

Status: (0% COMPLETO)

A combinação dos dois casos anteriores. Algumas coisas interessantes acontecem quando lidamos com forças externas e atrito (ex: ressonância e caos).

CASO 2.4: PÊNDULO FÍSICO (Caso extra!)

Status: (0% COMPLETO)

Um caso que talvez eu trabalhe. Trata-se do movimento de oscilação de objetos rígidos como uma régua ou uma porta. Por questão de complexidade, eu iria me limitar a casos bem básicos (uma régua ou um celular).

Grupo 3: Casos Avançados
CASO 3.1: PÊNDULO DUPLO

Status: (0% COMPLETO)

Um caso famoso que trabalha com movimento caótico. Tem várias propriedades interessantes que podem ser vistas realizando simulações de pêndulo.

CASO 3.2 PÊNDULO ELÁSTICO

Status: (0% COMPLETO)

Um caso que utiliza as equações de Euler-Lagrange para descrever um movimento complexo do pêndulo sendo atuado por uma mola.