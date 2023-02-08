%EM570 - Projeto 2
%Leonardo Santiago, Gabriel Beltrami, Eduardo Smarieri Xavier

clear variables;
clc;

%% Iremos resolver a equa��o de difus�o do calor da forma:
%  
% dT/dt = alfa*d2T/dx2
%
% Pelo m�todo de Euler expl�cito
%
% Problema 1: - Placa plana, temperaturas fixas nas extremidades;
%             - Condi��o inicial: T(x,t=0) = Ti;
%             - Condi��es de contorno: T(x=0,t) = T1 , T(x=L,t) = T2
%
% Problema 2: - Placa plana, convec��o nas superf�cies;
%             - Condi��o inicial: T(x,t=0) = Ti;
%             - Condi��es de contorno: h(T_inf-T(x)) = -KdT/dx, em x=0,x=L
%
% Problema 3: - Placa plana, semi-infinita;
%             - Condi��o inicial: T(x,t=0) = Ti;
%             - Condi��es de contorno: T(x=0,t) = T1 , T(x=L,t) = T2
%
% Problema 4: - Placa plana, Distribui��o de temp. inicial senoidal;
%             - Condi��o inicial: T(x,t=0) = sen(mx) + T1 + (T2-T1)*(x/L);
%             - Condi��es de contorno: Escolhida
%%
%Propriedades:
h = 1100; %[W/m^2.K]
k = 30; %[W/m.K]
cp = 600; %[J/kg.K]
ro = 1e4; %[kg/m^3]
alfa = k/(ro*cp); %[m^2/s]
Ti = 25; %[�C]
T_inf = 250; %[�C]
q = 0; %[W]  Gera��o de energia

%Dimens�es do problema:
L = 0.01; %[m]
n = 21; % Numero de pontos da malha
dx = L/(n-1); % Discretiza��o da dimens�o x
x = 0:dx:L; % Malha discretizada

%Par�metros adimensionais: (o crit�rio mais restritivo � Fo*(1+Bi)<= 1/2)
Bi = h*dx/k; %               N�mero de Biot
Fo_max = 0.5/(1+Bi);%        Valor limite para o n�mero de Fourier
if Fo_max >= 0.5
   Fo_max = 0.49; 
end
dt_max = Fo_max*dx^2/alfa; % Valor limite para delta t em fun��o de Fo
dt = 0.01; %[s]            Discretiza��o da malha temporal
Fo = alfa*dt/dx^2; %         N�mero de Fourier
tf = 300; %[s]               Intervalo de tempo analisado
p = tf/dt; %                 Numero de pontos na malha temporal

%%
%Distribui��o de temperatura em t = 0 (Condi��o inicial):
type_T0 = 'uniforme'; % Tipo de distribui��o inicial a ser escolhida

switch type_T0
    case 'uniforme'
        T0 = 300*ones([1,n]); %[�C] Distribui��o inicial de temperatura
    case 'senoidal'
        T1 = Ti; %[�C]
        T2 = T1; %[�C]
        m = 1000000*L; % Fator da senoide
        T0 = T1*(sin(m*x/L) +1.5); %[�C]
    case 'parabolica' %Nesse caso, h� gera��o de energia na parede
        qi = 1e7; %[W] Fonte de gera��o de energia inicial
        Ts1 = 340.91; %[�C] Temperatura em x = 0
        Ts2 = Ts1; %[�C] Temperatura na extremidade x=L
        T0 = qi*(L/2)^2/(2*k)*(1-(x-L/2).^2/(L/2)^2) + (Ts2-Ts1)*(x-(L/4))/(2*L) + (Ts1+Ts2)/2;
end
%%
%Determinando T(x,t): 
type_cc = 'T_fixa'; % Condi��o de contorno escolhida
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
        ylabel('T [�C]','fontweight','bold','fontsize',12)
        ylim([0,1000])
        pause( 0.1 );       
        end 
    case 'convec��o'
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
        ylabel('T [�C]','fontweight','bold','fontsize',12)
        ylim([300,500])
        pause( 1 );     
        end

end
        