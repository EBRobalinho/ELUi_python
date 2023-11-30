import numpy as np
import matplotlib.pyplot as plt

def Pot_I(concreto, m, e):
    n = concreto.n
    pot_const = concreto.sigma_cd * ((e**(m + 1)) / (m + 1))
    prod = concreto.sigma_cd

    if e < 0:
        pot = 0
    elif 0 <= e <= concreto.e_c2:
        pot = pot_const
        for i in range(1, m + 2):
            pot += (-1)**(i - 1) * np.math.comb(m, i - 1) * (((concreto.e_c2 - e)**(n + i) * concreto.sigma_cd * concreto.e_c2**(m - n - i + 1)) / (n + i))
        for i in range(1, m + 2):
            prod *= (concreto.e_c2) / (n + i)
        pot = pot - np.math.factorial(m) * prod
    elif e > concreto.e_c2:
        pot = pot_const
        for i in range(1, m + 2):
            prod *= (concreto.e_c2) / (n + i)
        pot = pot - np.math.factorial(m) * prod

    return pot

def Pot_J(concreto, m, def_, e_0, k):
    return concreto.Tensao_deformacao(def_) * (def_ - e_0)**m

def esf_res_aco(aco, posicao_aco, area_aco, e_0, k):
    N_s = 0
    M_s = 0

    for i in range(len(posicao_aco)):
        N_s += area_aco[i] * aco.Tensao_deformacao(e_0 + k * posicao_aco[i])
        M_s += area_aco[i] * aco.Tensao_deformacao(e_0 + k * posicao_aco[i]) * posicao_aco[i]

    return N_s, M_s

def esf_res_concreto(concreto, b, h, yt, yb, e_0, k):
    S = 0  # momento estático usado no momento fletor resistente

    if abs(k) * h < 1e-3:
        N_c = concreto.Tensao_deformacao(e_0) * b * h  # b * h representa a área da secção retangular
        M_c = concreto.Tensao_deformacao(e_0) * S
    else:
        N_c = b * (Pot_I(concreto, 0, e_0 + yt * k) - Pot_I(concreto, 0, e_0 + yb * k)) / k
        M_c = b * (
            (Pot_I(concreto, 1, e_0 + yt * k) - Pot_I(concreto, 1, e_0 + yb * k)) -
            e_0 * (Pot_I(concreto, 0, e_0 + yt * k) - Pot_I(concreto, 0, e_0 + yb * k))
        ) / k**2

    return N_c, M_c

def esf_estrutura(concreto, b, h, yt, yb, steel, yc, Area, e_0, k):
    N_s, M_s = esf_res_aco(steel, yc, Area, e_0, k)
    N_c, M_c = esf_res_concreto(concreto, b, h, yt, yb, e_0, k)
    N_r = N_s + N_c
    M_r = M_s + M_c

    return N_r, M_r


def vetor_esfor(concreto, b, h, yt, yb, steel, yc, Area, vec_e0, vec_k):
    N, M = esf_estrutura(concreto, b, h, yt, yb, steel, yc, Area, vec_e0, vec_k)

    return N, M

def Jacobian_s(aco, posicao_aco, area_aco, e_0, k):
    J_s = np.zeros((2, 2))

    for i in range(len(area_aco)):
        J_s[0, 0] -= 1 * aco.D_s(e_0 + k * posicao_aco[i]) * area_aco[i]
        J_s[1, 0] -= 1 * aco.D_s(e_0 + k * posicao_aco[i]) * area_aco[i] * posicao_aco[i]
        J_s[0, 1] = J_s[1, 0]
        J_s[1, 1] -= 1 * aco.D_s(e_0 + k * posicao_aco[i]) * area_aco[i] * posicao_aco[i]**2

    return J_s

def Jacobian_c(concreto, b, h, yt, yb, e_0, k):
    J_c = np.zeros((2, 2))

    if abs(k) * h < 1e-5:
        J_c = -1 * concreto.D_c(e_0) * np.array([[b * h, 0], [0, b * h**3 / 12]])
    else:
        J_c[0, 0] = -1 * b * (Pot_J(concreto, 0, e_0 + k * yt, e_0, k) - Pot_J(concreto, 0, e_0 + k * yb, e_0, k)) / k
        J_c[1, 0] = -1 * b * (
            Pot_J(concreto, 1, e_0 + k * yt, e_0, k) - Pot_J(concreto, 1, e_0 + k * yb, e_0, k) -
            (Pot_I(concreto, 0, e_0 + yt * k) - Pot_I(concreto, 0, e_0 + yb * k))
        ) / k**2
        J_c[0, 1] = J_c[1, 0]
        J_c[1, 1] = -1 * b * (
            (Pot_J(concreto, 2, e_0 + k * yt, e_0, k) - Pot_J(concreto, 2, e_0 + k * yb, e_0, k)) -
            2 * ((Pot_I(concreto, 1, e_0 + yt * k) - Pot_I(concreto, 1, e_0 + yb * k)) -
                 e_0 * (Pot_I(concreto, 0, e_0 + yt * k) - Pot_I(concreto, 0, e_0 + yb * k)))) / k**3

    return J_c

def NR(nr, concreto, aco, posicao_aco, area_aco, b, h, yt, yb, Ns, Ms):
    e_0 = nr.chute_inicial[0]
    k = nr.chute_inicial[1]

    # Criação do vetor do vetor dos esforços iniciais
    N, M = vetor_esfor(concreto, b, h, yt, yb, aco, posicao_aco, area_aco, e_0, k)
    # Verificação da norma admensionalizada para os esforços iniciais
    f_ad = np.sqrt(((N - Ns) / (concreto.sigma_cd * b * h))**2 + ((M - Ms) / (concreto.sigma_cd * b * h**2))**2)

    if f_ad <= nr.tol_norm:
        print("Chute inicial já serve como solução para o par de deformações")
    else:
        # Criação da matriz jacobiana com as deformações iniciais
        Js = Jacobian_s(aco, posicao_aco, area_aco, e_0, k)
        Jc = Jacobian_c(concreto, b, h, yt, yb, e_0, k)
        J = Js + Jc
        det_J_ad = 100
        # Critérios do início do loop
        i = 1
        Sol = np.array([[i, f_ad, e_0, k]])
        # Loop
        while (f_ad >= nr.tol_norm) and (i <= nr.tol_int):
            f_x = np.array([Ns, Ms]) - np.array([N, M])
            # Para avalia a matriz invertida na situação por pontos
            invJ = np.array([[J[1, 1], -1 * J[0, 1]], [-1 * J[0, 1], J[0, 0]]]) / np.linalg.det(J)
            var_t = np.dot(invJ, f_x)

            # Interação sobre as variáveis
            var = np.linalg.solve(J, f_x)
            e_0 = e_0 - var[0]
            k = k - var[1]

            # Cálculo do novo jacobiano
            Js = Jacobian_s(aco, posicao_aco, area_aco, e_0, k)
            Jc = Jacobian_c(concreto, b, h, yt, yb, e_0, k)
            J = Js + Jc
            # Determinante da matriz adimensionalizada para interar no loop
            det_J_ad = (
                J[0, 0] * J[1, 1] / (((concreto.sigma_cd * b * h)**2) * (h**2)) -
                (J[0, 1] / (concreto.sigma_cd * b * h**2))**2
            )
            # Testar condicionamento, ver se o sistema tem solução
            if det_J_ad <= nr.tol_det:
                print('Não Existe Solução!')
                Sol = np.array([[0, 0, 0, 0]])
                break
            else:
                # Criação do vetor do vetor dos esforços
                N, M = vetor_esfor(concreto, b, h, yt, yb, aco, posicao_aco, area_aco, e_0, k)
                # Norma adimensionalizada do vetor para interar no loop
                f_ad = np.sqrt(((N - Ns) / (concreto.sigma_cd * b * h))**2 + ((M - Ms) / (concreto.sigma_cd * b * h**2))**2)

                # Guardar interações
                i = i + 1
                Sol = np.vstack([Sol, [i, f_ad, e_0, k]])

    return Sol


def deformacoes(nr, concreto, steel, yc, Area, b, h, yt, yb, N_0, M_0):
    v = NR(nr, concreto, steel, yc, Area, b, h, yt, yb, N_0, M_0)
    e0 = v[-1, 2]
    k = v[-1, 3]
    return e0, k


def deformacoes(nr, concreto, steel, yc, Area, b, h, yt, yb, N_0, M_0):
    v = NR(nr, concreto, steel, yc, Area, b, h, yt, yb, N_0, M_0)
    e0 = v[-1, 2]
    k = v[-1, 3]
    return e0, k

def Dif_fin_ver(nr, concreto, steel, yc, Area, b, h, yt, yb, m, f_inic, L, N_0, M_0, tol_f):
    y = np.zeros(m)
    curvatura = np.zeros(m)
    M = np.zeros(m)
    delta_L = L / m

    # Cálculo dos esforços na secção 0
    M[0] = M_0
    e0, k = deformacoes(nr, concreto, steel, yc, Area, b, h, yt, yb, N_0, M[0])
    iterador = 1

    while iterador == 1:
        if (e0 != 0) or (k != 0):
            for i in range(1, m+1):
                if (e0 != 0) or (k != 0):
                    if i == 1:
                        curvatura[0] = k / 1000  # curvatura 0
                        y[0] = curvatura[0] * (delta_L**2) * 0.5
                        M[0] = M_0 + N_0 * (f_inic - y[0])
                        sol = [N_0, M[0], y[0], curvatura[0], 1]
                    else:
                        curvatura[i-1] = k / 1000
                        if i == 2:
                            y[i-1] = curvatura[i-1] * (delta_L**2) + 2 * y[i-2]
                        else:
                            y[i-1] = curvatura[i-1] * (delta_L**2) + 2 * y[i-2] - y[i-3]
                        M[i-1] = M_0 + N_0 * (f_inic - y[i-1])
                        sol = np.vstack([sol, [N_0, M[i-1], y[i-1], curvatura[i-1], i]])
                else:
                    print(f"Secção {i} não aguentou as solicitações")
                    sol = np.array([])
                    iterador = 0
                    break

                e0, k = deformacoes(nr, concreto, steel, yc, Area, b, h, yt, yb, N_0, M[i-1])

        else:
            print("Secção inicial não aguentou as solicitações")
            sol = np.array([])
            iterador = 0

        if abs(y[-1] - f_inic) < tol_f:
            iterador = 0
        else:
            f_inic = y[-1]
            y = np.zeros(m)
            curvatura = np.zeros(m)
            M = np.zeros(m)
            M[0] = M_0 + N_0 * (f_inic)
            e0, k = deformacoes(nr, concreto, steel, yc, Area, b, h, yt, yb, N_0, M[0])

    return sol

def Plot_Trajetoria_eq(nr, concreto, steel, yc, Area, b, h, yt, yb, m, f_inic, L, tol_f,e,nu_min,nu_max):
    N=np.zeros(100)
    M=np.zeros(100)
    f=np.zeros(100)
    i=0
    for nu in np.linspace(nu_min,nu_max,100):    # Adimensionais
        mi = nu*e #Adimensional do momento fletor solicitante
        # Verificação por diferenças finitas
        M_0 = concreto.sigma_cd * b * h * h * mi
        N_0 = concreto.sigma_cd * b * h * nu
        sol=Dif_fin_ver(nr, concreto, steel, yc, Area, b, h, yt, yb, m, f_inic, L, N_0, M_0, tol_f)
        if sol!=[]:
            N[i]=sol[-1, 0]
            M[i]=sol[-1, 1]
            f[i]=sol[-1, 2]
            i=i+1
        else:
            break
    return [N,M,f]

def Plot_Trajetoria_elui(nr, concreto, steel, yc, Area, b, h, yt, yb, m, f_inic, L, tol_f,e_min,e_max):
    N=np.zeros(100)
    M=np.zeros(100)
    f=np.zeros(100)
    i=0
    for e in np.linspace(e_min,e_max,50):    # Adimensionais
        [Npt, Mpt,fpt] = Plot_Trajetoria_eq(nr, concreto, steel, yc, Area, b, h, yt, yb, m, f_inic, L, tol_f,e,0.001,6)
        N[i]=ultimo_valor_nao_nulo(Npt)
        M[i]=ultimo_valor_nao_nulo(Mpt)
        f[i]=e
        i=i+1
    return [N,M,f]


def ultimo_valor_nao_nulo(vetor):
    for i in range(len(vetor) - 1, -1, -1):
        if vetor[i] != 0:
            return vetor[i]