import numpy as np
import matplotlib.pyplot as plt



class Aco:
    def __init__(self, Res_aco):
        self.fyk = Res_aco * 10
        self.fyd = self.fyk / 1.15
        self.E_s = 210 * 10**3
        self.gamma_s = 1.15
        self.e_yd = 1000 * self.fyd / self.E_s
        self.Densidade_aco = 7.850  # g/cm^3

    def Tensao_deformacao(self, deform):
        if 0 <= abs(deform) <= self.e_yd:
            sigma = self.fyd * deform / self.e_yd
        elif abs(deform) >= self.e_yd:
            sigma = self.fyd * deform / abs(deform)
        else:
            sigma = 0
        return sigma

    def D_s(self, deform):
        if 0 <= abs(deform) <= self.e_yd:
            sigma = self.E_s / 1000
        elif abs(deform) >= self.e_yd:
            sigma = 0
        else:
            sigma = None
        return sigma


class Concreto:
    def __init__(self, res_concreto):
        self.fck = float(res_concreto)
        self.gamma_c = 1.4
        self.Densidade_concreto = 2.5  # g/cm^3
        self.sigma_cd = 0.85 * self.fck / self.gamma_c

        if self.fck <= 50:
            self.e_c2 = 2
            self.e_cu = 3.5
            self.n = 2
        else:
            self.e_c2 = 2 + 0.085 * (self.fck - 50)**0.53
            self.e_cu = 2.6 + 35 * ((90 - self.fck) / 100)**4
            self.n = 1.4 + 23.4 * ((90 - self.fck) / 100)**4

    def Tensao_deformacao(self, deform):
        if deform < 0:
            sigma = 0
        elif 0 <= deform <= self.e_c2:
            sigma = self.sigma_cd * (1 - (1 - (deform / self.e_c2))**(self.n))
        else:
            sigma = self.sigma_cd
        return sigma

    def D_c(self, deform):
        if deform < 0:
            sigma = 0
        elif 0 <= deform <= self.e_c2:
            sigma = self.sigma_cd * (self.n / self.e_c2) * ((1 - (deform / self.e_c2))**(self.n - 1))
        else:
            sigma = 0
        return sigma

    def Plot_Tensao_deformacao(self):
        size = 100000
        sigma_plot = np.zeros(size)
        e = np.zeros(size)

        for k in range(1, 201):
            e[k] = k / size
            if e[k] < 0:
                sigma_plot[k] = 0
            elif 0 <= e[k] <= self.e_c2:
                sigma_plot[k] = self.sigma_cd * (1 - (1 - (e[k] / self.e_c2))**(self.n))
            else:
                sigma_plot[k] = self.sigma_cd

        plt.plot(e, sigma_plot, '*')
        plt.show()


class Diagrama:
    def __init__(self, concreto, posicao_aco, yt, yb, h):
        u = concreto.e_cu
        c = concreto.e_c2
        self.A = np.array([0, c])
        self.B = np.array([u/h, u*yt/h])
        self.C = np.array([(u+10)/(yt-min(posicao_aco)), u-yt*(u+10)/(yt-min(posicao_aco))])
        self.D = np.array([0, -10])
        self.E = np.array([(u+10)/(yb-max(posicao_aco)), u-yb*(u+10)/(yb-max(posicao_aco))])
        self.F = np.array([-u/h, -u*yb/h])

    def Testar_pontos(self, concreto, posicao_aco, h, yt, yb, e_0, k):
        i = 0
        print(f'Teste dos pontos com e_0 = {e_0:.3f} e k = {k:.3f}')
        if k >= 0:
            if e_0 + 10 + min(posicao_aco) * k >= 0:
                print('Passou no limite de ductibilidade')
                i = 1
            else:
                print('Não passou no limite de ductibilidade')
                
            if concreto.e_cu - yt * k - e_0 >= 0:
                print('Passou no limite de deformação máxima do concreto e_cu')
                if i == 1:
                    i = 2
            else:
                print('Não passou no limite de deformação máxima do concreto e_cu')
                
            if concreto.e_c2 - (yb + (concreto.e_c2 / concreto.e_cu) * h) * k - e_0 >= 0:
                print('Passou no limite de deformação máxima do concreto e_c2')
                if i == 2:
                    i = 3
            else:
                print('Não passou no limite de deformação máxima do concreto e_c2')
                
            if i == 3:
                print("A estrutura aguenta no ELU as deformações solicitadas")
            else:
                print("A estrutura não aguenta no ELU as deformações solicitadas")

        if k < 0:
            if e_0 + 10 + max(posicao_aco) * k >= 0:
                print('Passou no limite de ductibilidade')
                i = 1
            else:
                print('Não passou no limite de ductibilidade')

            if concreto.e_cu - yb * k - e_0 >= 0:
                print('Passou no limite de deformação máxima do concreto e_cu')
                if i == 1:
                    i = 2
            else:
                print('Não passou no limite de deformação máxima do concreto e_cu')

            if concreto.e_c2 - (yt + (concreto.e_c2 / concreto.e_cu) * h) * k - e_0 >= 0:
                print('Passou no limite de deformação máxima do concreto e_c2')
                if i == 2:
                    i = 3
            else:
                print('Não passou no limite de deformação máxima do concreto e_c2')

            if i == 3:
                print("A estrutura aguenta no ELU as deformações solicitadas")
            else:
                print("A estrutura não aguenta no ELU as deformações solicitadas")

    def diagrama_def(self, n_pontos):
        mesh = np.linspace(0, 1, n_pontos)
        retas = np.zeros((len(mesh), 2))

        for i in range(len(mesh)):
            if 0 <= mesh[i] <= 1/6:
                retas[i, :] = 6 * self.A * (1/6 - mesh[i]) + 6 * self.B * mesh[i]

            if 1/6 < mesh[i] <= 2/6:
                retas[i, :] = 6 * self.B * (2/6 - mesh[i]) + 6 * self.C * (mesh[i] - 1/6)

            if 2/6 < mesh[i] <= 3/6:
                retas[i, :] = 6 * self.C * (3/6 - mesh[i]) + 6 * self.D * (mesh[i] - 2/6)

            if 3/6 < mesh[i] <= 4/6:
                retas[i, :] = 6 * self.D * (4/6 - mesh[i]) + 6 * self.E * (mesh[i] - 3/6)

            if 4/6 < mesh[i] <= 5/6:
                retas[i, :] = 6 * self.E * (5/6 - mesh[i]) + 6 * self.F * (mesh[i] - 4/6)

            if 5/6 < mesh[i] <= 6/6:
                retas[i, :] = 6 * self.F * (6/6 - mesh[i]) + 6 * self.A * (mesh[i] - 5/6)

        return retas

class NewtonRaphson:
    def __init__(self, tol_norm, tol_int, tol_det, chute_inicial, dimension):
        self.tol_norm = tol_norm
        self.tol_int = tol_int
        self.tol_det = tol_det
        self.dimension = dimension
        self.chute_inicial = chute_inicial
