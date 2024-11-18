%% Progetto Controlli Automatici T B1
% Cremonini Davide, Nanni Gabriele, Tarlassi Federica

clear all; close all; clc;
% Paramentri Tecnologici
r_s=1.5; %tasso riproduzione cellule suscettibili
r_r=1.3; %tasso riproduzione cellule resistenti
K=200; %num massimo cellule dell'ambiente
gamma=0.3; %costante di passaggio da r a s
beta=0.3; %costante di passaggio da s a r
alfa=0.2; %costante di mutazione a seguito del farmaco
m_s=0.7; %mortalità suscettibili
m_r=0.15; %mortalità resistenti
n_s_e=100; %equilibrio suscettibili
n_r_e=100; %equilibrio resistenti

% Equazioni Differenziali
%n_s'=-r_s*ln((n_s+n_r)/K)*n_s-m_s*c_f*n_s-beta*n_s+gamma*n_r-alfa*c_f*n_s
%n_r'=-r_r*ln((n_s+n_r)/K)*n_r-m_r*c_f*n_r+beta*n_s-gamma*n_r+alfa*c_f*n_s

%% Punto1

%x1'(t)=f1(x(t),u(t))=-r_s*ln((x1(t)+x2(t))/K)*x1(t)-m_s*u(t)*x1(t)-beta*x1(t)+gamma*x2(t)-alfa*u(t)*x1(t)
%x2'(t)=f2(x(t),u(t))=-r_r*ln((x1(t)+x2(t))/K)*x2(t)-m_r*u(t)*x2(t)+beta*x1(t)-gamma*x2(t)+alfa*u(t)*x1(t)
%y(t)=h(x(t),u(t))=x2(t)

%u(t)=(-r_s*ln((x1(t)+x2(t))/K)*x1(t)-beta*x1(t)+gamma*x2(t))/((m_s+alfa)*x1(t))
%u(t)=(-r_r*ln((x1(t)+x2(t))/K)*x2(t)+beta*x1(t)-gamma*x2(t))/((m_r*x2(t)-alfa*x1(t))
%con xe1=n_s_e=100   xe2=n_r_e=100
%u_e_1=(-r_s*log((n_s_e+n_r_e)/K)*n_s_e-beta*n_s_e+gamma*n_r_e)/((m_s+alfa)*n_s_e)
%u_e_2=(-r_r*log((n_s_e+n_r_e)/K)*n_r_e+beta*n_s_e-gamma*n_r_e)/(m_r*n_r_e-alfa*n_s_e)
%abbiamo verificato che esiste u_e=u_e_1=u_e_2=0
u_e=0;
xe1=n_s_e;
xe2=n_r_e;

%linearizzazione nell'intorno dell'equilibrio

A=[(-r_s/2)-beta, (-r_s/2)+gamma; -(r_r/2)+beta, -(r_r/2)-gamma];
B=[-m_s*xe1-alfa*xe1; -m_r*xe2+alfa*xe1];
C=[0, 1];
D=0;

%% Punto2
%funzione di trasferimento
s=tf('s');
GG=C*(inv(s*eye(2)-A))*B-D;

%% Punto3

% SPECIFICHE
% 
% -) Errore a regime in risposta a un gradino w(t) = -2*1(t), d(t) = sommatoria k=1->4 0.3sen(0.0125kt) e n(t) = sommatoria k=1->4 0.2sen(10^4kt) 2*1(t) pari a 0
%
% -) Attenuazione di almeno 60dB per d(t)
%       con [omega_d_min, omega_d_MAX] = [0,0.05]
%
% -) Attenuazione di almeno 90dB per n(t)
%       con [omega_n_min, omega_n_MAX] = [10^4,10^6]
%
% -) S% <= 7%
% -) Ta,5 <= 1 s
%

%%Scrivere specifiche

% ampiezze gradini
WW = -2;

% errore a regime
e_star =0;

% attenuazione disturbo sull'uscita
A_d = 60;
omega_d_min =1e-4;
omega_d_MAX = 0.05;

% attenuazione disturbo di misura
A_n = 90;
omega_n_min = 1e4;
omega_n_MAX = 1e6;

% Sovraelongazione massima e tempo d'assestamento all'5%
S_100_spec = 0.07;
T_a5_spec = 1;

%% Diagramma di Bode

figure(1);
bode(GG);
grid on, zoom on;
if 0        %%%%%%%%%%%%%%%%%%%%---SEMAFORO Figure 1---%%%%%%%%%%%%%%%%%%%%
    return;
end
   
%% Patch specifiche

figure(2);
grid on; zoom on; hold on;

%scrivere calcolo specifiche S% => Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100;
omega_Ta_MAX = 300/(Mf_spec*T_a5_spec); 

% Specifiche su d
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n (massima pulsazione di taglio)
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = 1e-4; % lower bound per il plot
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -300; -300];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}", "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
omega_plot_min = 1e-4;
omega_plot_max = 1e6;
margin(GG, {omega_plot_min, omega_plot_max});
grid on; zoom on;



% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

if 0        %%%%%%%%%%%%%%%%%%%%---SEMAFORO Figure 2---%%%%%%%%%%%%%%%%%%%%
    return;
end

%% Regolatore statico

%scrivere valore minimo prescritto per L(0)
mu_s = 200;

%scrivere guadagno minimo del regolatore ottenuto come L(0)/G(0)
G_0 = abs(evalfr(GG,0));
RR_s = (mu_s/G_0)/s; 

% Sistema esteso
GG_e=RR_s*GG;


%% Diagrammi di Bode di Ge con specifiche

figure(3);
hold on;

%scrivere calcolo specifiche S% => Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100;
omega_Ta_MAX = 300/(Mf_spec*T_a5_spec); 

% Specifiche su d
omega_d_min = 1e-4; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n (massima pulsazione di taglio)
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = 1e-4; % lower bound per il plot
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -300; -300];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}", "Ge(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
omega_plot_min = 1e-4;
omega_plot_max = 1e6;
margin(GG_e, {omega_plot_min, omega_plot_max});
grid on; zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["Ge(j\omega)"; "M_f"];
legend(Legend_arg);

% FINE regolatore statico

if 0 %%%%%%%%%%%%%%%%%%%%---SEMAFORO Figure 3---%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return;
end

% Inserimento di uno zero a bassa frequenza

RR_lf = (1 + 5*s);
GG_zero = GG_e * RR_lf;

figure(4);
hold on;

%scrivere calcolo specifiche S% => Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100;
omega_Ta_MAX = 300/(Mf_spec*T_a5_spec); 

% Specifiche su d
omega_d_min = 1e-4; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n (massima pulsazione di taglio)
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = 1e-4; % lower bound per il plot
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -300; -300];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}", "Ge(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
omega_plot_min = 1e-4;
omega_plot_max = 1e6;
margin(GG_zero, {omega_plot_min, omega_plot_max});
grid on; zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["Ge(j\omega)"; "M_f"];
legend(Legend_arg);

if 0        %%%%%%%%%%%%%%%%%%%%---SEMAFORO Figure 4---%%%%%%%%%%%%%%%%%%%%
    return;
end

% Inserimento di poli in alta frequenza

RR_hf = 1/(1 + (1/2000)*s)^4;
GG_zero_polo = GG_e * RR_lf * RR_hf;

figure(5);
hold on;

%scrivere calcolo specifiche S% => Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100;
omega_Ta_MAX = 300/(Mf_spec*T_a5_spec); 

% Specifiche su d
omega_d_min = 1e-4; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n (massima pulsazione di taglio)
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = 1e-4; % lower bound per il plot
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -300; -300];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}", "Ge(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
omega_plot_min = 1e-4;
omega_plot_max = 1e6;
margin(GG_zero_polo, {omega_plot_min, omega_plot_max});
grid on; zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["Ge(j\omega)"; "M_f"];
legend(Legend_arg);

RR = RR_s*RR_lf*RR_hf;
LL = GG_zero_polo;

if 0  %%%%%%%%%%%%%%%%%%%%---SEMAFORO Figure 5---%%%%%%%%%%%%%%%%%%%%
    return;
end

%% Check prestazioni in anello chiuso

% Funzione di sensitività complementare
FF = LL/(1+LL);

% Risposta al gradino

figure(6)
T_simulation = 2;
[y_step,t_step] = step(WW*FF, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[WW*(1+S_100_spec),WW*(1+S_100_spec),WW-1,WW-1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
ylim([WW-1,0]);

% vincolo tempo di assestamento all'5%

LV = -2; 

patch([T_a5_spec,T_simulation,T_simulation,T_a5_spec],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5_spec,T_simulation,T_simulation,T_a5_spec],[LV*(1+0.05),LV*(1+0.05),LV-1,LV-1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

% Diagramma di bode sistema in anello chiuso
figure(7);
bode(FF)
grid on, zoom on

if 0  %%%%%%% SEMAFORO controllo gradino %%%%%%%
    return;
end

%% Check disturbo in uscita


% Funzione di sensitività
SS = 1/(1+LL);
figure(8);

% Simulazione disturbo d(t)
omega_d_test = 0.0125;
tt = (0:1e-2:1e3)';
dd = 0.3*(sin(omega_d_test*tt)+sin(omega_d_test*2*tt)+sin(omega_d_test*3*tt)+sin(omega_d_test*4*tt));

y_d = lsim(SS,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('dd','y_d')

% Simulazione disturbo n(t)
omega_n_test = 1e4;
tt = (0:1e-5:1e-2)';
nn = 0.2*(sin(omega_n_test*tt)+sin(omega_n_test*2*tt)+sin(omega_n_test*3*tt)+sin(omega_n_test*4*tt));

figure(9)

y_n = lsim(-FF,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('nn','y_n')

% simulazione y(t) completa

tt = (0:1e-4:1)';
% DD = 0.3*(omega_d_test/(s^2+omega_d_test^2)+omega_d_test*2/(s^2+(omega_d_test*2)^2)+omega_d_test*3/(s^2+(omega_d_test*3)^2)+omega_d_test*4/(s^2+(omega_d_test*4)^2));
% NN = 0.2*(omega_n_test/(s^2+omega_n_test^2)+omega_n_test*2/(s^2+(omega_n_test*2)^2)+omega_n_test*3/(s^2+(omega_n_test*3)^2)+omega_n_test*4/(s^2+(omega_n_test*4)^2));

dd = 0.3*(sin(omega_d_test*tt)+sin(omega_d_test*2*tt)+sin(omega_d_test*3*tt)+sin(omega_d_test*4*tt));
nn = 0.2*(sin(omega_n_test*tt)+sin(omega_n_test*2*tt)+sin(omega_n_test*3*tt)+sin(omega_n_test*4*tt));
y_d_test = lsim(SS,dd,tt);
y_n_test = lsim(-FF,nn,tt);
y_w_test = -2*step(FF,tt);
 
yy = y_w_test + y_d_test + y_n_test;
 
 figure(10)
 hold on, grid on, zoom on
 rif = -2*heaviside(tt); % Heaviside comando per la rappresentazione del gradino nel tempo
 plot(tt,rif,'b');
 plot(tt,yy,'r');
 hold on, grid on, zoom on
 legend('W','y')

%% Punto 5

% Andamento sistema non linearizzato non regolato
% input applicato: costante e pari a 0
inp = @(t) 0;

% intervallo di tempo
interv = [0 10]; % da 0 a 10 secondi

dyn = @(t, x) [-r_s*log((x(1)+x(2))/K)*x(1)-m_s*inp(t)*x(1)-beta*x(1)+gamma*x(2)-alfa*inp(t)*x(1);
    -r_r*log((x(1)+x(2))/K)*x(2)-m_r*inp(t)*x(2)+beta*x(1)-gamma*x(2)+alfa*inp(t)*x(1)];

% vettore delle condizioni iniziali
x0 = [95; 102];

% risolviamo l'equazione differenziale
[time, traj] = ode45(dyn, interv, x0);

figure(11)
plot(time,traj)
title('Traiettoria delle cellule cancerose')
xlim(interv)
xlabel('tempo [s]')
ylabel('stato')
legend('cellule sensibili', 'cellule resistenti')
grid on; zoom on; box on;

%RR antitrasformato restituisce delta ingresso e delta risposta in
%relazione all'ingresso e all'uscita di equilibrio, mentre il sistema non
%linearizzato ragiona in modo assoluto, quindi è necessario ricordare la
%discrepanza di riferimento.

%% continua su Simulink

sample_time = 0.0001;