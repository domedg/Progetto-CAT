%% Progetto Controlli Automatici T B1
% Barone Leonardo, Del Giudice Domenico, Galli Francesco, Guzzonato Leonardo

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Paramentri Tecnologici %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_s=1.7; %tasso riproduzione cellule suscettibili
r_r=1.4; %tasso riproduzione cellule resistenti
K=500; %num massimo cellule dell'ambiente
gamma=0.2; %costante di passaggio da r a s
beta=0.8; %costante di passaggio da s a r
alfa=0.5; %costante di mutazione a seguito del farmaco
m_s=0.95; %mortalità suscettibili
m_r=0.05; %mortalità resistenti
n_s_e=100; %equilibrio suscettibili
n_r_e=400; %equilibrio resistenti

u_e=0;
xe1=n_s_e;
xe2=n_r_e;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Linearizzazione nell'intorno dell'equilibrio %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[r_s*(1-(2*xe1+xe2)/K)-m_s*u_e-beta-alfa*u_e, (r_s*xe1)/K+gamma;
    (r_r*xe2)/K+beta+alfa*u_e, r_r*(1-(xe1+2*xe2)/K)-m_r*u_e-gamma];
A
B=[-m_s*xe1-alfa*xe1; -m_r*xe2+alfa*xe1];
B
C=[0, 1];
C
D=0;
D




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PUNTO 2 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%funzione di trasferimento
modello = ss(A,B,C,D);
GG = tf(modello);

GG
% Calcolo dei poli e degli zeri della funzione di trasferimento
p = pole(GG);
z = zero(GG);

% Stampa dei poli
fprintf('Poli della funzione di trasferimento:\n');
disp(p);

% Stampa degli zeri
fprintf('Zeri della funzione di trasferimento:\n');
disp(z);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PUNTO 3 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Patch specifiche
% 3.1 Errore a regime nullo 
WW = -110;
DD = -20;
e_star = 0;

% 3.2 Margine di fase
Mf_esp = 45;

% 3.3 Sovraelongazione percentuale massima 10%  
S_100_spec = 0.1;

% 3.4 Tempo d'assestamento al epsilon% = 5%
epsilon = 0.05;
T_star = 0.2;

% 3.5attenuazione disturbo sull'uscita
A_d = 45; %valore in Db
omega_d_min = 0;
omega_d_MAX = 0.1;

% attenuazione disturbo di misura
A_n = 80; %% valore in Db
omega_n_min = 5*10^3;
omega_n_MAX = 5*10^6;


%% Diagramma di Bode
figure(1);
bode(GG);
grid on, zoom on;

%3.2 scrivere calcolo specifiche S% => Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100;
% Uso 300 perchè ln(1/0.05) = circa 3, deriva dalla formula del punto 4.
omega_Ta_MAX = 300/(Mf_spec*T_star); 
omega_Ta_min = 1e-4; % lower bound per il plot

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

%% Patch specifiche
figure(2);
hold on;

% Specifiche su d
omega_plot_min = 1e-4;
omega_plot_max = 1e6;
Bnd_d_x = [omega_plot_min; omega_d_MAX; omega_d_MAX; omega_plot_min];
Bnd_d_y = [A_d;A_d;-100;-100];


% Specifiche su n (massima pulsazione di taglio)
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -300; -300];


patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG, {omega_plot_min, omega_plot_max});
grid on; zoom on;

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);


% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

% Imposta limiti sugli assi
xlim([1e-4, 1e6]); % Limiti della frequenza (asse X)
ylim([-200, 200]); % Limiti della magnitudine (asse Y)


%print('-depsc', 'bode_patch.eps');
%Fine figura 4: Diagramma di Bode con patch delle zone proibite 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Regolatore statico %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scrivere valore minimo prescritto per L(0)
mu_s = 1.5e4;

%scrivere guadagno minimo del regolatore ottenuto come L(0)/G(0)
G_0 = abs(evalfr(GG,0));
s=tf('s');
RR_s = (mu_s/G_0)/s; 
%RR_s = mu_s/s; 

% Sistema esteso
GG_e=RR_s*GG;

%% Diagrammi di Bode di Ge con specifiche
figure(3);
hold on;

%3.2 scrivere calcolo specifiche S% => Margine di fase
logsq = (log(S_100_spec))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = xi*100;

% Uso 300 perchè ln(1/0.05) = circa 3, deriva dalla formula del punto 4.
omega_Ta_MAX = 300/(Mf_spec*T_star); 
omega_Ta_min = 1e-4; % lower bound per il plot


% Specifiche su d
omega_plot_min = 1e-4;
omega_plot_max = 1e6;
Bnd_d_x = [omega_plot_min; omega_d_MAX; omega_d_MAX; omega_plot_min];
Bnd_d_y = [A_d;A_d;-100;-100];


% Specifiche su n (massima pulsazione di taglio)
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -300; -300];


patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e, {omega_plot_min, omega_plot_max});
grid on; zoom on;

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);


% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

% Imposta limiti sugli assi
xlim([1e-4, 1e6]); % Limiti della frequenza (asse X)
ylim([-200, 200]); % Limiti della magnitudine (asse Y)

print('-depsc', 'ge_bode_patch.eps');

%%%%%%%% FINE regolatore statico %%%%%%%%



