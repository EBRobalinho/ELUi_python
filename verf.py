from Classes import *
from funcoes import *



# Dimensões da seção retangular
b = 0.2  # em metros
h = 0.5  # em metros

# Posições extremas da seção transversal
yt = h / 2
yb = -h / 2

# Definição do tipo de concreto
concreto = Concreto(20)

# Definição do tipo de aço
steel = Aco(50)  # CA50, colocar , 50 mas o aço aguenta 500 Mpa

# Posição das camadas de aço
dt = 0.025  # distância da armadura à borda da seção
nc = 2  # número de camadas de aço
yc = [h / 2 - dt, -h / 2 + dt]  # posição das camadas de aço
n_barras = np.array([3, 3])  # número de barras por camada
Area = 0.25 * (16 ** 2) * (10 ** -6) * np.pi * n_barras  # área de aço por camada

Area_num = 0.4 * b * h * concreto.sigma_cd / steel.fyd
Area = [(Area_num) / 2, (Area_num) / 2]

# Diagrama no espaço de deformações da seção
dig = Diagrama(concreto, yc, yt, yb, h)

# Cálculo da Verificação via Newton Raphson
nr = NewtonRaphson(10 ** -10, 1000, 10 ** -10, [0, 0],2)

# Adimensionais
nu = 0.3
mi = 0.1

# Verificação por diferenças finitas
M_0 = concreto.sigma_cd * b * h * h * mi
m = 100
L = 5
N_0 = concreto.sigma_cd * b * h * nu
f_inic = 0
tol_f = 10 ** -10
sol = Dif_fin_ver(nr, concreto, steel, yc, Area, b, h, yt, yb, m, f_inic, L, N_0, M_0, tol_f)

# Plotar resultados
plt.plot(100 * sol[:, 2] / h, sol[:, 4] * L / m, '*')
plt.show()