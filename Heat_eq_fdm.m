%EM570 - Projeto 2
%Leonardo Santiago, Gabriel Beltrami, Eduardo Smarieri Xavier

clear variables;
clc;

%% Iremos resolver a equação de difusão do calor da forma:
%  
% dT/dt = alfa*d2T/dx2
%
% Pelo método de Euler explícito
%
% Problema 1: - Placa plana, temperaturas fixas nas extremidades;
%             - Condição inicial: T(x,t=0) = Ti;
%             - Condições de contorno: T(x=0,t) = T1 , T(x=L,t) = T2
%
% Problema 2: - Placa plana, convecção nas superfícies;
%             - Condição inicial: T(x,t=0) = Ti;
%             - Condições de contorno: h(T_inf-T(x)) = -KdT/dx, em x=0,x=L
%
% Problema 3: - Placa plana, semi-infinita;
%             - Condição inicial: T(x,t=0) = Ti;
%             - Condições de contorno: T(x=0,t) = T1 , T(x=L,t) = T2
%
% Problema 4: - Placa plana, Distribuição de temp. inicial senoidal;
%             - Condição inicial: T(x,t=0) = sen(mx) + T1 + (T2-T1)*(x/L);
%             - Condições de contorno: Escolhida
%%
%Propriedades:
h = 1100; %[W/m^2.K]
k = 30; %[W/m.K]
cp = 600; %[J/kg.K]
ro = 1e4; %[kg/m^3]
alfa = k/(ro*cp); %[m^2/s]
Ti = 25; %[ºC]
T_inf = 250; %[ºC]
q = 0; %[W]  Geração de energia

%Dimensões do problema:
L = 0.01; %[m]
n = 21; % Numero de pontos da malha
dx = L/(n-1); % Discretização da dimensão x
x = 0:dx:L; % Malha discretizada

%Parâmetros adimensionais: (o critério mais restritivo é Fo*(1+Bi)<= 1/2)
Bi = h*dx/k; %               Número de Biot
Fo_max = 0.5/(1+Bi);%        Valor limite para o número de Fourier
if Fo_max >= 0.5
   Fo_max = 0.49; 
end
dt_max = Fo_max*dx^2/alfa; % Valor limite para delta t em função de Fo
dt = 0.01; %[s]            Discretização da malha temporal
Fo = alfa*dt/dx^2; %         Número de Fourier
tf = 300; %[s]               Intervalo de tempo analisado
p = tf/dt; %                 Numero de pontos na malha temporal

%%
%Distribuição de temperatura em t = 0 (Condição inicial):
type_T0 = 'uniforme'; % Tipo de distribuição inicial a ser escolhida

switch type_T0
    case 'uniforme'
        T0 = 300*ones([1,n]); %[ºC] Distribuição inicial de temperatura
    case 'senoidal'
        T1 = Ti; %[ºC]
        T2 = T1; %[ºC]
        m = 1000000*L; % Fator da senoide
        T0 = T1*(sin(m*x/L) +1.5); %[ºC]
    case 'parabolica' %Nesse caso, há geração de energia na parede
        qi = 1e7; %[W] Fonte de geração de energia inicial
        Ts1 = 340.91; %[ºC] Temperatura em x = 0
        Ts2 = Ts1; %[ºC] Temperatura na extremidade x=L
        T0 = qi*(L/2)^2/(2*k)*(1-(x-L/2).^2/(L/2)^2) + (Ts2-Ts1)*(x-(L/4))/(2*L) + (Ts1+Ts2)/2;
end
%%
%Determinando T(x,t): 
type_cc = 'T_fixa'; % Condição de contorno escolhida
T = T0;

switch type_cc
    case 'T_fixa'     
        T(1) = 0;
        T(n) = T(1); 
        for i=1:p
             Ta = T;
             for m = 2:n-1
                  T(m) = Fo*(Ta(m-1)+Ta(m+1)+q*dx^2/k)+Ta(m)*(1-2*Fo); 
             end
        plot(x,T,'-oc','LineWidth',3)
        grid on;
        xlabel('x [m]','fontweight','bold','fontsize',12)
        ylabel('T [ºC]','fontweight','bold','fontsize',12)
        ylim([0,1000])
        pause( 0.1 );       
        end 
    case 'convecção'
        for i=1:p
             Ta = T;
             T(1) = Ta(1)*(1-2*Fo-2*Fo*Bi)+2*Fo*(Bi*T_inf+Ta(2)+q*dx^2/(2*k));
             for m = 2:n-1
                  T(m) = Fo*(Ta(m-1)+Ta(m+1)+q*dx^2/k)+Ta(m)*(1-2*Fo); 
             end
             T(n) = Ta(n)*(1-2*Fo-2*Fo*Bi)+2*Fo*(Bi*T_inf+Ta(n-1)+q*dx^2/(2*k));
        plot(x,T,'-om','LineWidth',3)
        grid on;
        xlabel('x [m]','fontweight','bold','fontsize',12)
        ylabel('T [ºC]','fontweight','bold','fontsize',12)
        ylim([300,500])
        pause( 1 );     
        end

end
        