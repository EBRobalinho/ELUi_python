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
yc = [h / 2 - dt, -h / 2 + dt]  # posição das camadas de aço
n_barras = np.array([3, 3])  # número de barras por camada

omega=0.4 #Adimensional da área de aço total na estrutura

Area_num = omega * b * h * concreto.sigma_cd / steel.fyd #Área de aço an estrutura em m²

Area = [(Area_num) / 2, (Area_num) / 2]

# Diagrama no espaço de deformações da seção
dig = Diagrama(concreto, yc, yt, yb, h)

# Cálculo da Verificação via Newton Raphson
nr = NewtonRaphson(10 ** -10, 1000, 10 ** -10, [0, 0],2)

# Adimensionais
nu = 0.3 #Adimensional do esforço solicitante normal
mi = 0.1 #Adimensional do momento fletor solicitante

# Verificação por diferenças finitas
M_0 = concreto.sigma_cd * b * h * h * mi
N_0 = concreto.sigma_cd * b * h * nu

m = 10 #Número de elementos na diferença finita
L = 5 #Comprimento do pilar

f_inic = 0
tol_f = 10 ** -10
sol = Dif_fin_ver(nr, concreto, steel, yc, Area, b, h, yt, yb, m, f_inic, L, N_0, M_0, tol_f)


d=dig.diagrama_def(12000)
N=np.zeros(len(d))
M=np.zeros(len(d))
for i in range(0,len(d)):
    [N[i],M[i]]=vetor_esfor(concreto, b, h, yt, yb, steel, yc, Area, d[i,0], d[i,1])

plt.plot(M[:]/(concreto.sigma_cd * b * h * h ), N[:]/(concreto.sigma_cd * b * h) , '--')
plt.plot(sol[:, 1]/(concreto.sigma_cd * b * h*h),sol[:, 0]/(concreto.sigma_cd * b * h), '*')
plt.show()

# Plotar resultados
plt.plot(100 * sol[:, 2] / h, sol[:, 4] * L / m, '*')
plt.show()