%% Instituto Tecnológico de Aeronaútica
% Aluno: Eduardo Bezerra Robalinho Dantas da Gama;
% *Verificação na flexão normal composta*
% 
% 
% Reiniciando as variáveis de ambiente

clear all
%% 
% Dimensões da secção Retangular:

b=0.2; %em metros
h=0.5; %em metros
%% 
% Posições extremas da secção transversal

yt=h/2;
yb=-h/2; 
%% 
% Definição do tipo de concreto

concreto=Concreto(20);
%% 
% Definição do tipo de Aço:

steel=Aco(50);% CA50, colocar , 50 mas o aço aguenta 500 Mpa
diametro=16; % em milimetros;
qtd_aco=4; %quantidade total de barras de aço
%% 
% Posição das camadas de aço:

dt=0.025; %distância da armadura a borda da secção 
nc=2; %número de camadas de aço
yc=[h/2-dt,-h/2+dt]; %posição das camadas de aço
n_barras=[3,3]; %número de barras por camada;
Area=0.25*(diametro^2)*(10^-6)*pi*n_barras; %área de aço por camada

Area_num=0.4*b*h*concreto.sigma_cd/steel.fyd;
Area=[(Area_num)/2,(Area_num)/2];
%% 
% Diagrama no espaço de deformações da secção

dig=Diagrama(concreto,yc,yt,yb,h); %yt é positivo, yb é negativo
%% 
% Cálculo da Verificação via Newton Raphson

nr=Newton_Raphson(10^-10,1000,10^-10,[0,0],2);
%% 
% Adimensionais:

nu=0.3;
mi=0.1;
%% 
% Verificação por diferenças finitas:

M_0=concreto.sigma_cd*b*h*h*mi;
m=5;
L=5;
N_0=concreto.sigma_cd*b*h*nu;
f_inic=0;
tol_f=10^-10;
sol=Dif_fin_ver(nr,concreto,steel,yc,Area,b,h,yt,yb,m,f_inic,L,N_0,M_0,tol_f);
plot(100*sol(:,3)/h,sol(:,5)*L/m,'*');
%% 
% Função deformada do Pilar:

function sol=Dif_fin_ver(nr,concreto,steel,yc,Area,b,h,yt,yb,m,f_inic,L,N_0,M_0,tol_f)
    y=zeros(m,1);
    curvatura=zeros(m,1);
    M=zeros(m,1);
    delta_L=L/m;
    %Cálculo dos esforços na secção 0
    M(1)=M_0;
    [e0,k]=deformacoes(nr,concreto,steel,yc,Area,b,h,yt,yb,N_0,M(1));
    iterador=1;
while iterador==1

    if (e0~=0) || (k~=0) 
        for i=1:m
            if (e0~=0) || (k~=0) 
                if i==1
                    curvatura(1)=k/1000; %curvatura 0
                    y(1)=curvatura(1)*(delta_L^2)*0.5;
                    M(1)=M_0 + N_0*(f_inic-y(1));
                    sol=[N_0,M(1),y(1),curvatura(1),1];
                else
                    curvatura(i)=k/1000;
                    if i==2
                        y(i)=curvatura(i)*delta_L^2 + 2*y(i-1);
                    else
                        y(i)=curvatura(i)*delta_L^2 + 2*y(i-1) - y(i-2);
                    end
                    M(i)=M_0 + N_0*(f_inic-y(i));
                    sol=[sol;N_0,M(i),y(i),curvatura(i),i];
                end
            else
                fprintf("Secção %d não aguentou as solicitações",i);
                sol=[];
                iterador=0;
                break
            end
            [e0,k]=deformacoes(nr,concreto,steel,yc,Area,b,h,yt,yb,N_0,M(i));
        end
    else
        fprintf("Secção inicial não aguentou as solicitações");
        sol=[];
        iterador=0;
    end  

    if abs(y(m,1)-f_inic)<tol_f
        iterador=0;
    else
        f_inic=y(m,1);
        y=zeros(m,1);
        curvatura=zeros(m,1);
        M=zeros(m,1);
        M(1)=M_0 + N_0*(f_inic);
        [e0,k]=deformacoes(nr,concreto,steel,yc,Area,b,h,yt,yb,N_0,M(1));
    end

end
end
%% Funções chamadas no cálculo dos esforços

function [e0,k]=deformacoes(nr,concreto,steel,yc,Area,b,h,yt,yb,N_0,M_0)
    v=NR(nr,concreto,steel,yc,Area,b,h,yt,yb,N_0,M_0);
    e0=v(end,3);
    k=v(end,4);
end
%% 
% Função Potencial I

function pot=Pot_I(concreto,m,e)
    n=concreto.n;
    pot_const=concreto.sigma_cd*((e^(m+1))/(m+1));
    prod=concreto.sigma_cd;
    if e<0
        pot=0;
    end
    if (e>=0) && (e<=concreto.e_c2)
        pot=pot_const;
        for i=1:m+1
            pot=pot+(-1)^(i-1)*nchoosek(m,i-1)*(((concreto.e_c2-e)^(n+i)*concreto.sigma_cd*concreto.e_c2^(m-n-i+1))/(n+i));
        end
        for i=1:m+1
            prod=prod*(concreto.e_c2)/(n+i);
        end
        pot=pot-factorial(m)*prod;
    end
    if (e>concreto.e_c2)
        pot=pot_const;
        for i=1:m+1
            prod=prod*(concreto.e_c2)/(n+i);
        end
        pot=pot-factorial(m)*prod;
    end
end

%% 
% Função potencial J:

function potj=Pot_J(concreto,m,def,e_0,k)
    potj=concreto.Tensao_deformacao(def)*(def-e_0)^m;
end

%% 
% Esforços resistentes no aço:

function [N_s,M_s]=esf_res_aco(aco,posicao_aco,area_aco,e_0,k)
    N_s=0;
    M_s=0;
    for i=1:length(posicao_aco)
        N_s=N_s+area_aco(i)*aco.Tensao_deformacao(e_0+k*posicao_aco(i));
        M_s=M_s+area_aco(i)*aco.Tensao_deformacao(e_0+k*posicao_aco(i))*posicao_aco(i);
    end
end
%% 
% Esforços resistentes no concreto:

function [N_c,M_c]=esf_res_concreto(concreto,b,h,yt,yb,e_0,k)
    S=0; %momento estático usado no momento fletor resistente
    if abs(k)*h<10^-3
        N_c=concreto.Tensao_deformacao(e_0)*b*h;  %b*h representa a área da secção retangular
        M_c=concreto.Tensao_deformacao(e_0)*S;
    else
        N_c=b*(Pot_I(concreto,0,e_0+yt*k)-Pot_I(concreto,0,e_0+yb*k))/k;
        M_c=b*((Pot_I(concreto,1,e_0+yt*k)-Pot_I(concreto,1,e_0+yb*k))-e_0*(Pot_I(concreto,0,e_0+yt*k)-Pot_I(concreto,0,e_0+yb*k)))/k^2;
    end
end
%% 
% Esforços resistentes na estrutura:

function [N_r,M_r]=esf_estrutura(concreto,b,h,yt,yb,steel,yc,Area,e_0,k)
    [N_s,M_s]=esf_res_aco(steel,yc,Area,e_0,k);
    [N_c,M_c]=esf_res_concreto(concreto,b,h,yt,yb,e_0,k);
    N_r=(N_s+N_c);
    M_r=(M_s+M_c);
end
%% 
% Esforços resistentes total por deformação:

function [N,M]=vetor_esfor(concreto,b,h,yt,yb,steel,yc,Area,vec_e0,vec_k)
    for i=1:length(vec_e0)
        [N(i),M(i)]=esf_estrutura(concreto,b,h,yt,yb,steel,yc,Area,vec_e0(i),vec_k(i));
    end
end
%% 
% Contribuição de esforços resistentes no aço:

function J_s=Jacobian_s(aco,posicao_aco,area_aco,e_0,k)
    J_s=zeros(2,2);
    for i=1:length(area_aco)
            J_s(1,1)=J_s(1,1)-1*aco.D_s(e_0 + k*posicao_aco(i))*area_aco(i);
            J_s(2,1)=J_s(2,1)-1*aco.D_s(e_0 + k*posicao_aco(i))*area_aco(i)*posicao_aco(i);
            J_s(1,2)=J_s(2,1);
            J_s(2,2)=J_s(2,2)-1*aco.D_s(e_0 + k*posicao_aco(i))*area_aco(i)*posicao_aco(i)^2;     
    end
end
%% 
% Contribuição de esforços resistentes no concreto:

function J_c=Jacobian_c(concreto,b,h,yt,yb,e_0,k)
    J_c=zeros(2,2);
    if abs(k)*h<10^-5
        J_c=-1*concreto.D_c(e_0)*[b*h,0;0,b*h^3/12];
    else
        J_c(1,1)=-1*b*(Pot_J(concreto,0,e_0+k*yt,e_0,k)-Pot_J(concreto,0,e_0+k*yb,e_0,k))/k;
        J_c(2,1)=-1*b*(Pot_J(concreto,1,e_0+k*yt,e_0,k)-Pot_J(concreto,1,e_0+k*yb,e_0,k)-(Pot_I(concreto,0,e_0+yt*k)-Pot_I(concreto,0,e_0+yb*k)))/k^2;
        J_c(1,2)=J_c(2,1);
        J_c(2,2)=-1*b*((Pot_J(concreto,2,e_0+k*yt,e_0,k)-Pot_J(concreto,2,e_0+k*yb,e_0,k))-2*((Pot_I(concreto,1,e_0+yt*k)-Pot_I(concreto,1,e_0+yb*k))-e_0*(Pot_I(concreto,0,e_0+yt*k)-Pot_I(concreto,0,e_0+yb*k))))/k^3;
    end
end
%% 
% Newton-Raphson

function Sol = NR(nr,concreto,aco,posicao_aco,area_aco,b,h,yt,yb,Ns,Ms)

    e_0 = nr.chute_inicial(1);
    k = nr.chute_inicial(2);

    %Criação do vetor do vetor dos esforços iniciais
    [N,M]=vetor_esfor(concreto,b,h,yt,yb,aco,posicao_aco,area_aco,e_0,k);
    %Verificação da norma admensionalizada para os esforços iniciais
    f_ad=sqrt( ((N-Ns)/(concreto.sigma_cd*b*h))^2  +  ((M-Ms)/(concreto.sigma_cd*b*h^2))^2);


    if f_ad<=nr.Tol_norm
        fprintf("Chute inicial já serve como solução para o par de deformações");
    else
        %Criação da matriz jacobiana com as deformações iniciais
        Js=Jacobian_s(aco,posicao_aco,area_aco,e_0,k);
        Jc=Jacobian_c(concreto,b,h,yt,yb,e_0,k);
        J=(Js+Jc);
        det_J_ad=100;
        %Critérios do início do loop
        i=1;
        Sol=[i,f_ad,e_0,k];
        %Loop
        while (f_ad>=nr.Tol_norm) && (i<=nr.Tol_int)     
            f_x=(([Ns,Ms]-[N,M])');
            %Para avalia a matriz invertida na situação por pontos 
            invJ=[J(2,2),-1*J(1,2);-1*J(1,2),J(1,1)]/det(J);
            var_t=invJ*f_x;
            
            %Interação sobre as variáveis
            var=(J)\(([Ns,Ms]-[N,M])');         %0.73*10^-2* %fator de escala
            e_0=e_0-var(1,1);
            k=k-var(2,1);
    
            %Cálculo do novo jacobiano
            Js=Jacobian_s(aco,posicao_aco,area_aco,e_0,k);
            Jc=Jacobian_c(concreto,b,h,yt,yb,e_0,k);
            J=(Js+Jc);
            %Determinante da matriz adimensionalizada para interar no loop
            det_J_ad=( J(1,1)*J(2,2)/(((concreto.sigma_cd*b*h)^2)*(h^2)) - (J(1,2)/(concreto.sigma_cd*b*h^2))^2 );
            %Testar condicionamento, ver se o sistema tem solução
            if det_J_ad<=nr.Tol_det
                fprintf('Não Existe Solução!');
                Sol=[0,0,0,0];
                break
            else
                %Criação do vetor do vetor dos esforços
                [N,M]=vetor_esfor(concreto,b,h,yt,yb,aco,posicao_aco,area_aco,e_0,k);
                %Norma adimensionalizada do vetor para interar no loop
                f_ad=sqrt( ((N-Ns)/concreto.sigma_cd*b*h)^2  +  ((M-Ms)/concreto.sigma_cd*b*h^2)^2);
    
                %Guardar interações
                i=i+1;
                Sol=[Sol;i,f_ad,e_0,k];
            end
        end
    
    
    end

    end
%% 
%
