% Progetto C1 di Controlli Automatici T 
% Enrico Tornabene, Cristiano Pelliccia, Filippo Valmaggi
% 
% Ultima Modifica: 16/01/24 18:00
%
%
% SPECIFICHE
% 
% -) Errore a regime in risposta a un gradino w(t) = 1(t) e d(t) = 1(t)
%    pari a 0.02
%
% -) Attenuazione di almeno 30dB per d(t)
%       con [omega_d_min, omega_d_MAX] = [0,0.05]
%
% -) Attenuazione di almeno 75dB per n(t)
%       con [omega_n_min, omega_n_MAX] = [8*10^4,8*10^6]
%
% -) Epsilon% = 5%
% -) S% <= 5%
% -) Ta,1 <= 0.01 s
%
% M_f >= 30 gradi
%--------------------------------------------------------------------------

clear all; close all; clc

%dati per la visualizzione, pulsazione minima e massima
omega_plot_min = 1e-4;
omega_plot_max = 1e6;


% Parametri Fisici:
m = 2500;   % Massa
beta = 0.1; % Coefficiente di non linearità della molla
n = 3;      % Coefficiente di non linearità della molla
k = 5 * 1e5;% Costante elastica lineare della molla
b = 350;    % Coefficiente di attrito dinamico
g = 9.81;   % Acc. Gravitazionale 

% Calcolo coppia di equilibrio:
x1_e = 0.35; 
x2_e = 0; 
u_e = m*9.81 + k*x1_e;

% Linearizzazione e calcolo della FdT
A = [0 1; -k/m -b/m]
B = [0 ; 1/m]
C = [1 0]
D = 0

s = tf('s');
G = C * inv(s*eye(2) - A) * B + D


%% Definizione Specifiche

% errore a regime
WW = 1;
DD = 1;
e_star = 0.02;

% attenuazione disturbo sull'uscita
A_d = 30; %decibel max
omega_d_min = 1e-4;
omega_d_MAX = 0.05;

% attenuazione disturbo di misura
A_n = 75; %decibel max 
omega_n_min = 8*1e4; 
omega_n_MAX = 8*1e6;

% Sovraelongazione massima e tempo d'assestamento al 5%
S_star = 5;
T_star = 1e-2; %Ta,epsilon

% Margine di fase
Mf_esp = 30; %in gradi

%% Regolatore statico - senza poli nell'origine, agendo solo sul guadagno statico 

% valore minimo prescritto per L(0)
mu_s_error = (DD+WW)/e_star; %errore
mu_s_dist  = 10^(A_d/20);    %vincolo disturbo uscita

% guadagno minimo del regolatore ottenuto come L(0)/G(0)
G_0 = abs(evalfr(G,0));
G_omega_d_MAX = abs(evalfr(G,j*omega_d_MAX));
R_s = max(mu_s_error/G_0, mu_s_dist/G_omega_d_MAX); 

% Sistema esteso
G_e = R_s*G;

%% Diagrammi di Bode di Ge con specifiche
figure(2);
hold on;

% Calcolo specifiche S% => Margine di fase
xi_star = abs( log(S_star/100))/sqrt(pi^2 + log(S_star/100)^2);
Mf = max(xi_star*100, Mf_esp); %Mf_esp da specifica

% Specifiche su d (disturbo in uscita) 
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n (rumore di misura) 
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione critica)
omega_Ta_min = 1e-4;           % per plottare
omega_Ta_max = 460/(Mf*T_star) % omega_c >= 460/(Mf*T^*) ~ 4.6
Bnd_Ta_x = [omega_Ta_min; omega_Ta_max; omega_Ta_max; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150]; 
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(G_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% Specifiche sovraelongazione (margine di fase), patch verde
omega_c_min = omega_Ta_max;
omega_c_max = omega_n_min;

%range
phi_up = Mf - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y = [phi_up; phi_up; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);


% siamo nello scenario B 

%% Regolatore dinamico

Mf_star = Mf+5;   % per stare larghi
omega_c_star = 800; % loop-shaping

%per calcolo dei coefficienti del regolatore dinamico
mag_omega_c_star_dB = abs(evalfr(G_e,j*omega_c_star));              
arg_omega_c_star    = rad2deg(angle(evalfr(G_e,j*omega_c_star)));   

M_star = 1/mag_omega_c_star_dB;
phi_star = Mf_star - 180 - arg_omega_c_star

%formule d'inversione
tau = (M_star - cosd(phi_star))/(omega_c_star*sind(phi_star));          %cosd(x) = cos(x*pi/2)
alpha_tau = (cosd(phi_star) - 1/M_star)/(omega_c_star*sind(phi_star));
alpha = alpha_tau / tau;

if min(tau,alpha) < 0
    fprintf('Errore: parametri rete anticipatrice negativi');
    return;
end

%Per rispettare vincolo su disturbo uscita e margine di fase
% polo ad altra frequenza
R_high_frequency = 1/(1 + s/(3e4));

% regolatore dinamico = rete anticipatrice * polo ad alta frequenza
R_d = (1 + tau*s)/(1 + alpha * tau*s)*R_high_frequency;

% regolatore 
RR = R_s*R_d;

LL = RR*G; % funzione ad anello aperto

% grafico per vedere che le specifiche vengano rispettate
figure(3);
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_arg);

%% Check prestazioni in anello chiuso 

% Funzione di sensitività complementare
FF = LL/(1+LL);

% Risposta al gradino
figure(4);

T_simulation = 2*T_star; % T_simulation = 0.02
[y_step,t_step] = step(WW*FF, T_simulation); % test segnale a gradito su anello chiuso
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

LV = evalfr(WW*FF,0)

% vincolo sovraelongazione 
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'5% 
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% Check disturbo in uscita

% Funzione di sensitività
SS = 1/(1+LL);
figure(5);

% Simulazione disturbo in uscita a pulsazione 0.01
omega_d = 0.01;
tt = 0:1e-2:2e3;
% dd = DD*sin(omega_d*1*tt) + DD*sin(omega_d*2*tt) + DD*sin(omega_d*3*tt) + DD*sin(omega_d*4*tt);
dd = DD * sum( sin(omega_d * (1:4)' * tt) , 1); % sommatoria 
y_d = lsim(SS,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('d(t)','y_d(t)')

%% Check disturbo di misura
figure(6);

% Simulazione disturbo di misura 
omega_n = omega_n_min; % da range specifica disturbo misura 8*10^4
NN      = 0.2;         % ampiezza 
tt = 0:1e-5:2*1e-3; 
% nn = NN*sin(omega_n*tt) + NN*sin(omega_n*2*tt) + NN*sin(omega_n*3*tt) + NN*sin(omega_n*4*tt);
nn = NN * sum( sin(omega_n * (1:4)' * tt), 1 ); %sommatoria

y_n = lsim(FF,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('n(t)','y_n(t')
%
%
%
%% Dati per Simulink: 
kk = 1:4;
%
%
%
%% PUNTO OPZIONALE 1
% figure(100);
% T_simulation = 5;
% [Y, ~] = step(WW * (G_0 / (1 + G_0)), T_simulation * 5);
% clf;
% y_veicolo = [1 1 0.9 0.9];
% x_veicolo = [-2 2 2 -2];
% S_100=0.05;
% for i=1:2:50
%     clf;
%     axis([-10 10  0 6]);
%     patch(x_veicolo, Y(i) + y_veicolo, 'r');
%     patch([-10, 10, 10, -10], [WW * (1 + S_100) + 1, WW * (1 + S_100) + 1, WW * (1 - S_100) + 1, WW * (1 - S_100) + 1], 'y', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.3);
%     grid on, zoom on, hold on;
%      plot([-10, -9, -9:((+9 +Y(i)-4)/9):Y(i)-4, Y(i)-4, Y(i)-2],...
%         [-10, -9, -9:((+9 +Y(i)-4)/9):Y(i)-4, Y(i)-4, Y(i)-2],'r','LineWidth',2) % spring
%     title('Risposta al gradino W 1(t)');
% 
%     pause(0.03);
% end
%
%
%
%% PUNTO OPZIONALE 2 
% Parametri della simulazione
simTime = 0.2;  % Tempo di simulazione per ciascuna condizione iniziale
numSteps = 1;   % Passi
% range plot [plotMin, plotMax]
plotMin = -4;   
plotMax = 4;

% TEST SU x1
% Inizializza la figura
figure(7);
hold on;

% Aggiungiamo il valore con x1 = 0.35 per avere un riferimento
x1_value = 0.35;
x2_value = 0;
WW = 0; 
simOut = sim('schema_di_controlloC1_Pelliccia_Tornabene_Valmaggi_PtiOpz.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');
outputData = simOut.y_nlin; %prelevo il valore di ritorno (matrice con coppia tempo:valore) 
plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('x1 = %g', x1_value)); %plotto sul valore di ritorno

for i = plotMin:numSteps:plotMax
    % Aumenta gradualmente il valore di x1
    x1_value = i;
    x2_value=x2_e;
    WW=0;

    % Simulazione con il valore di x1 corrente
    simOut = sim('schema_di_controlloC1_Pelliccia_Tornabene_Valmaggi_PtiOpz.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');

    % Valore di ritorno
    outputData = simOut.y_nlin;

    % Plotto il grafico per il valore di x1 attuale
    plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('x1 = %d', x1_value));
end

% per grafico
title('Variazione di x1 nel tempo');
xlabel('Time');
ylabel('Amplitude');
legend('show');
grid on;
hold off;
%
%
%
% TEST su x2
numSteps = 5;  
plotMin = -20;
plotMax = 20;
simTime = 0.5;
% Inizializza la figura
figure(8);
hold on;

for i = plotMin:numSteps:plotMax
    % Aumenta gradualmente il valore di x2 mantenendo x1 = x1e
    x1_value = x1_e;
    x2_value=i;
    WW=0;
    
    %Simulo e plotto
    simOut = sim('schema_di_controlloC1_Pelliccia_Tornabene_Valmaggi_PtiOpz.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');
    outputData = simOut.y_nlin;
    plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('x2 = %d', x2_value));
end

% Per grafico
title('Variazione di x2 nel tempo');
xlabel('Time');
ylabel('Amplitude');
legend('show');
grid on;
hold off;
%
%
%
%% Punto opzionale 3
% Test su WW (ampiezza gradino input) 
% Inizializzo la figura
figure(9);
hold on;
simTime = 0.05
x1_value = x1_e;
x2_value=x2_e;
plotMin = 0;
plotMax = 10;
numSteps = 1;
for i = plotMin:numSteps:plotMax
    % Aumento WW
    WW=i;

    simOut = sim('schema_di_controlloC1_Pelliccia_Tornabene_Valmaggi_PtiOpz.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');
    outputData = simOut.y_nlin;
    plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('WW = %d', WW));
end

% Per grafico
title('Variazione di WW nel tempo');
xlabel('Time');
ylabel('Amplitude');
legend('show');
grid on;
hold off;


WW = 1;