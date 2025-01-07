%% Progetto Controlli Automatici T A3
%% GRUPPO 21 
%% Partecipanti: Barone Leonardo, Del Giudice Domenico, Galli Francesco, Guzzonato Leonardo

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

%Valori d'equilibrio
u_e=0;
xe1=n_s_e;
xe2=n_r_e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Linearizzazione nell'intorno dell'equilibrio %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[r_s*(1-(2*xe1+xe2)/K)-m_s*u_e-beta-alfa*u_e, -(r_s*xe1)/K+gamma;
    -(r_r*xe2)/K+beta+alfa*u_e, r_r*(1-(xe1+2*xe2)/K)-m_r*u_e-gamma];
A;
B=[-m_s*xe1-alfa*xe1; -m_r*xe2+alfa*xe1];
B;
C=[0, 1];
C;
D=0;
D;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PUNTO 2 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%funzione di trasferimento
modello = ss(A,B,C,D);
GG = tf(modello);
fprintf('Funzione di trasferimento G:\n');
GG;

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
Mf_star_esp = 45;

% 3.3 Sovraelongazione percentuale massima 10%  
S_star = 0.1;

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


%% Diagramma di Bode con poli e zeri
figure(1);

% Diagramma di Bode
[mag, phase, w] = bode(GG);
mag = squeeze(mag);
phase = squeeze(phase);

% Magnitudine
subplot(2, 1, 1);
semilogx(w, 20*log10(mag), 'b', 'LineWidth', 1.5); hold on;

% Aggiunta di poli e zeri (corretti)
poli_handle = [];  % Per tracciare i poli in legenda
zeri_handle = [];  % Per tracciare gli zeri in legenda

for k = 1:length(p)
    h = semilogx(abs(p(k)), interp1(w, 20*log10(mag), abs(p(k)), 'linear', 'extrap'), 'xr', 'MarkerSize', 10, 'LineWidth', 2); % Poli
    if isempty(poli_handle)
        poli_handle = h; % Salva il primo oggetto per legenda
    end
end

for k = 1:length(z)
    h = semilogx(abs(z(k)), interp1(w, 20*log10(mag), abs(z(k)), 'linear', 'extrap'), 'og', 'MarkerSize', 10, 'LineWidth', 2); % Zeri
    if isempty(zeri_handle)
        zeri_handle = h; % Salva il primo oggetto per legenda
    end
end

xlabel('Frequency (rad/s)');
ylabel('Ampiezza [dB]');
title('Diagramma di Bode - Magnitudine');
grid on;
ylim([-40 60]); % Scala simile alla figura

% Fase
subplot(2, 1, 2);
semilogx(w, phase, 'r', 'LineWidth', 1.5); hold on;

for k = 1:length(p)
    semilogx(abs(p(k)), interp1(w, phase, abs(p(k)), 'linear', 'extrap'), 'xr', 'MarkerSize', 10, 'LineWidth', 2); % Poli
end

for k = 1:length(z)
    semilogx(abs(z(k)), interp1(w, phase, abs(z(k)), 'linear', 'extrap'), 'og', 'MarkerSize', 10, 'LineWidth', 2); % Zeri
end

xlabel('Frequency (rad/s)');
ylabel('Fase [°]');
title('Diagramma di Bode - Fase');
grid on;
subplot(2, 1, 1);

% Inizializzazione della legenda
legend_items = {'Magnitudine'};
legend_handles = [semilogx(NaN, NaN, 'b', 'LineWidth', 1.5)]; % Handle per la Magnitudine

legend_handles(end+1) = semilogx(NaN, NaN, 'r', 'LineWidth', 1.5); % Handle per la Fase
legend_items{end+1} = 'Fase';

% Aggiunta dei poli
if ~isempty(poli_handle)
    legend_items{end+1} = 'Poli';
    legend_handles(end+1) = poli_handle;
end

% Aggiunta degli zeri 
if ~isempty(zeri_handle)
    legend_items{end+1} = 'Zeri';
    legend_handles(end+1) = zeri_handle;
end

legend(legend_handles, legend_items, 'Location', 'best');


%3.2 scrivere calcolo specifiche S% => Margine di fase
xi_star = abs(log(S_star))/sqrt(pi^2 + log(S_star)^2);
Mf_star = max(xi_star*100,Mf_star_esp);

% Uso 300 perchè ln(1/0.05) = circa 3, deriva dalla formula del punto 4.
omega_Ta_MAX = 300/(Mf_star*T_star); 
omega_Ta_min = 1e-4; % lower bound per il plot

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_star - 180;
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Regolatore statico %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=tf('s');
mu_s_d=(10^(A_d/20))/abs(evalfr(GG/s, 1i*omega_d_MAX));
RR_s = mu_s_d/s;
RR_s;

% Sistema esteso
GG_e=RR_s*GG;

%% Diagrammi di Bode di Ge con specifiche
figure(3);
hold on;

%3.2 scrivere calcolo specifiche S% => Margine di fase
xi_star = abs(log(S_star))/sqrt(pi^2 + log(S_star)^2);
Mf_star = max(xi_star*100,Mf_star_esp);

% Uso 300 perchè ln(1/0.05) = circa 3, deriva dalla formula del punto 4.
omega_Ta_MAX = 300/(Mf_star*T_star); 
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

phi_spec = Mf_star - 180;
phi_low = -270; % lower bound per il plot

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

% Imposta limiti sugli assi
xlim([1e-4, 1e6]); % Limiti della frequenza (asse X)
ylim([-200, 200]); % Limiti della magnitudine (asse Y)


%%%%%%%%%%%%%%%% FINE regolatore statico %%%%%%%%%%%%%%%%

% if 0 %%%%%%%%%%%%%%%%%%%%---SEMAFORO Figura 3---%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     return;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Regolatore Dinamico %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_c_star = 30;

mag_omega_c_star = abs(evalfr(GG_e,1i*omega_c_star));
arg_omega_c_star = rad2deg(angle(evalfr(GG_e,1i*omega_c_star)));

M_star = 1/mag_omega_c_star;
% Aggiunta di 5 gradi introdotta per margine di sicurezza
phi_star = Mf_star + 5 - 180 - arg_omega_c_star;

%%%%%%%%% Verifica Valori %%%%%%%%%
%Verifica che: Mf_star = 180 + arg(G_e(j*omega_c_star)) + phi_star
arg_Ge_omega_c_star = rad2deg(angle(evalfr(GG_e, 1j*omega_c_star)));
%Rimozione dei 5 gradi aggiunti, per verificare la condizione
tot = arg_Ge_omega_c_star + phi_star + 180 - 5;
condition1 = ['Verifica della condizione Mf_star = 180 + arg(G_e(j*omega_c_star)) + phi_star: ',num2str(Mf_star),' = ',num2str(tot)];
disp(condition1);

% Verifica che: |G_e(j*omega_c_star)|db + 20logM* = 0
Ge_dB = mag2db(abs(evalfr(GG_e, 1j*omega_c_star)));
M_star_db = 20*log10(M_star);
condition = Ge_dB + M_star_db; % Vera se essere uguale a 0 con tolleranza per numeri molto piccoli
condition2 = ['Verifica della condizione |G_e(j*omega_c_star)|_dB + 20logM* = 0: ',num2str(condition),'= 0'];
disp(condition2);

%FORMULE DI INVERSIONE
tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
alpha = alpha_tau / tau;

if min(tau,alpha) < 0
    fprintf('Errore: parametri rete anticipatrice negativi');
    return;
end

%Regolatore dinamico
RR_d = (1 + tau*s)/(1 + alpha * tau*s);
RR_d;

RR = RR_s*RR_d;
RR;

%Funzione d'anello
LL = RR*GG;
LL;
pole_LL = pole(LL);
%Stampa di debug(usata anche successivamente per verificare che le funzioni
%abbiano solo poli a parte reale negativa denotandone la stabilità)
disp('I poli della funzione L sono:');
disp(pole_LL);

figure(4);
hold on;

%[DEBUG] Indicatore che evidenzia la posizione di omega_c_star(deve essere
% sull' asse 0db
%xline(omega_c_star, '--r', '\omega_{c*}', 'LabelVerticalAlignment', 'bottom', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);

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


%%%%%%%% FINE Regolatore Dinamico %%%%%%%%

% SEMAFORO
if 0
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TEST SUL SISTEMA LINEARIZZATO %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prestazioni in anello chiuso

% Funzione di sensitività complementare
FF = LL/(1+LL);
FF;
pole_FF = pole(FF);
disp('I poli della funzione F sono:');
disp(pole_FF);

% Risposta al gradino
figure(5)
T_simulation = 2*T_star;
[y_step,t_step] = step(WW*FF, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

%Ampiezza del gradino w(t)
LV = -110;

% Imposta i limiti degli assi specificati
ylim([-250, 50]); % Limiti asse Y

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star),LV*(1+S_star),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

% Diagramma di bode sistema in anello chiuso
figure(6);
bode(FF)
grid on, zoom on

if 0  %%%%%%% SEMAFORO controllo gradino %%%%%%%
    return;
end

%% Check disturbo in uscita

% Funzione di sensitività
SS = 1/(1+LL);
figure(7);

% Simulazione disturbo d(t)
omega_d_test = 0.025;
%L’apice trasforma questo vettore riga in un vettore colonna (column vector).
tt = (0:1e-2:1e3)';
%Sommatoria da k=1 a 4
dd = -2*(sin(omega_d_test*tt)+sin(omega_d_test*2*tt)+sin(omega_d_test*3*tt)+sin(omega_d_test*4*tt));

y_d = lsim(SS,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('dd','y_d')

% Simulazione disturbo n(t)
omega_n_test = 5e3;
tt = (0:1e-5:1e-2)';

nn = 0.1*(sin(omega_n_test*tt)+sin(omega_n_test*2*tt)+sin(omega_n_test*3*tt)+sin(omega_n_test*4*tt));

figure(8)

y_n = lsim(-FF,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('nn','y_n')


% simulazione y(t) completa
tt = (0:1e-4:1)';
dd = -2*(sin(omega_d_test*tt)+sin(omega_d_test*2*tt)+sin(omega_d_test*3*tt)+sin(omega_d_test*4*tt));
nn = 0.1*(sin(omega_n_test*tt)+sin(omega_n_test*2*tt)+sin(omega_n_test*3*tt)+sin(omega_n_test*4*tt));
y_d_test = lsim(SS,dd,tt);
y_n_test = lsim(-FF,nn,tt);
y_w_test = WW*step(FF,tt);

yy = y_w_test + y_d_test + y_n_test;

figure(9)
hold on, grid on, zoom on

rif = WW*heaviside(tt); %Heaviside comando per la rappresentazione del gradino (di ampiezza WW) nel tempo
plot(tt,rif,'b');
plot(tt,yy,'r');
hold on, grid on, zoom on
legend('W','y')

%%%%%%%%% FINE TEST SISTEMA LINEARIZZATO %%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PUNTO OPZIONALE 1 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Andamento sistema non linearizzato non regolato
% input applicato (concentrazione del farmaco): 
%   -costante e pari a 0, numero sensibili e resistenti invariati
%   -costante e pari a 20, le cellule resistenti diventano 100(quanto le
%       sensibili iniziali e le sensibili si azzerano   
%   -costante e pari a 100, cellule sensibili e resistenti si annullano.

%% Grafico a barre "Evoluzione densità delle cellule"
inp = @(t, dose) dose;

% Intervallo di tempo e condizioni iniziali
interv = [0 10];

%Condizione iniziale
x0 = [100; 400];

% Dosaggi di farmaco da testare
dosi = [0, 2, 5, 7, 10, 13, 15, 17, 20];

% Array per memorizzare i risultati finali
densita_finale = zeros(length(dosi), 2); % [Sensibili, Resistenti]

for i = 1:length(dosi)
    dose = dosi(i);
    dyn = @(t, x) [
        r_s*(1-(x(1)+x(2))/K)*x(1) - m_s*inp(t, dose)*x(1) - beta*x(1) + gamma*x(2) - alfa*inp(t, dose)*x(1);
        r_r*(1-(x(1)+x(2))/K)*x(2) - m_r*inp(t, dose)*x(2) + beta*x(1) - gamma*x(2) + alfa*inp(t, dose)*x(1)
    ];

% Risoluzione dell'equazione differenziale
[time, traj] = ode45(dyn, interv, x0);

% Memorizzazione delle popolazioni finali
densita_finale(i, :) = traj(end, :);
end

figure(10)

bar(dosi, densita_finale, 'grouped');
title('Evoluzone della densità delle cellule')
xlabel('Dose di farmaco')
ylabel('Popolazione finale')
legend({'Cellule sensibili', 'Cellule resistenti'}, 'Location', 'northeast')
grid on; zoom on; box on;


%% Grafico con slider "Evoluzione densità delle cellule"

% Creazione della figura principale
f = figure('Name', 'Figura 11', ...
           'NumberTitle', 'off', ...
           'Position', [200, 200, 900, 600]);

% Creazione degli assi per il grafico
ax = axes('Parent', f, 'Position', [0.1 0.3 0.8 0.6]);
hold on; grid on;
sensibili_plot = plot(ax, nan, nan, 'b', 'LineWidth', 2); % Cellule sensibili
resistenti_plot = plot(ax, nan, nan, 'r', 'LineWidth', 2); % Cellule resistenti
legend(ax, 'Cellule sensibili', 'Cellule resistenti', 'Location', 'best');
xlabel(ax, 'Tempo [s]');
ylabel(ax, 'Densità delle cellule');
title(ax, 'Evoluzione della densità delle cellule');
xlim(ax, [0 5]);
ylim(ax, [0 500]);

% Creazione del testo dinamico per il valore dello slider
label = uicontrol('Parent', f, ...
                  'Style', 'text', ...
                  'Position', [370, 80, 160, 20], ...
                  'String', 'Concentrazione farmaco = 20', ...
                  'FontSize', 8, ...
                  'HorizontalAlignment', 'center');

% Creazione dello slider
slider = uicontrol('Parent', f, ...
                   'Style', 'slider', ...
                   'Min', 0, ...
                   'Max', 45, ...
                   'Value', 20, ...
                   'Position', [150, 50, 600, 20], ...
                   'Callback', @(src, event) aggiornaGrafico(src, sensibili_plot, resistenti_plot, label, interv, x0));

% Aggiornamento iniziale del grafico
aggiornaGrafico(slider, sensibili_plot, resistenti_plot, label, interv, x0);

% Funzione di aggiornamento del grafico
function aggiornaGrafico(slider, sensibili_plot, resistenti_plot, label, interv, x0)
    r_s=1.7; %tasso riproduzione cellule suscettibili
    r_r=1.4; %tasso riproduzione cellule resistenti
    K=500; %num massimo cellule dell'ambiente
    gamma=0.2; %costante di passaggio da r a s
    beta=0.8; %costante di passaggio da s a r
    alfa=0.5; %costante di mutazione a seguito del farmaco
    m_s=0.95; %mortalità suscettibili
    m_r=0.05; %mortalità resistenti
    % Valore corrente dello slider
    inp_val = get(slider, 'Value');
        
    % Aggiornamento dell'etichetta dinamica
    set(label, 'String', sprintf('Concentrazione farmaco = %.1f', inp_val));
    
    % Definizione della dinamica
    dyn = @(t, x) [r_s*(1-(x(1)+x(2))/K)*x(1) - m_s*inp_val*x(1) - beta*x(1) + gamma*x(2) - alfa*inp_val*x(1);
                   r_r*(1-(x(1)+x(2))/K)*x(2) - m_r*inp_val*x(2) + beta*x(1) - gamma*x(2) + alfa*inp_val*x(1)];
    
    % Risoluzione dell'equazione differenziale
    [time, traj] = ode45(dyn, interv, x0);
    
    % Aggiornamento dei dati del grafico
    set(sensibili_plot, 'XData', time, 'YData', traj(:, 1));
    set(resistenti_plot, 'XData', time, 'YData', traj(:, 2));
end

%%%% Fine grafico animazione "Evoluzione densità delle cellule" %%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PUNTO OPZIONALE 2 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Supponendo un riferimento w(t) ≡ 0, esplorare il range di condizioni iniziali dello stato del sistema
%non lineare (nell'intorno del punto di equilibrio) tali per cui l'uscita del sistema in anello chiuso
%converga a h(xe, ue).

%% Varizione della posizione inizale x1

% Parametri della simulazione
simTime = 0.5;  % Tempo di simulazione per ciascuna condizione iniziale
numSteps = 100;  % Numero di passi

figure(12);
hold on;

x1_value = xe1;
x2_value = xe2;
WW_opz = 0;
simOut = sim('simulink_puntiOpzionali.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');
outputData = simOut.y_nlin;
% Plot separato per migliorare la visualizzazione delle leggenda
plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('x1 = %d', 1));

for i = 100:numSteps:1000
    % Variazione di x1 con x2 costante
    x1_value = i;
    x2_value=xe2;
    WW_opz=0;

    % Simulazione con il valore di x1 corrente
    simOut = sim('simulink_puntiOpzionali.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');
    outputData = simOut.y_nlin;

    % Plot del grafico per il valore di x1 corrente
    plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('x1 = %d', x1_value));
end

title('Variazione di x1 nel tempo');
xlabel('Time');
ylabel('Amplitude');
xlim([0, 0.2]); % Limiti asse X
legend('show');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% Varizione della velocità inizale x2

% Parametri della simulazione
simTime = 0.25;  % Tempo di simulazione per ciascuna condizione iniziale
numSteps = 10;  % Numero di passi

figure(13);
hold on;
% Plot per il primo valore di x2 (386)
x1_value = xe1;  
x2_value = 386;  
WW_opz = 0;      
simOut = sim('simulink_puntiOpzionali.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');
outputData = simOut.y_nlin;
plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('x2 = %d', x2_value));

for i = 390:numSteps:500
    % Variazione di x2 con x1 costante
    x1_value=xe1;
    x2_value = i;
    WW_opz=0;

    % Simulazione con il valore di x2 corrente
    simOut = sim('simulink_puntiOpzionali.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');
    outputData = simOut.y_nlin;

    % Plot del grafico per il valore di x2 corrente
    plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('x2 = %d', x2_value));
end


title('Variazione di x2 nel tempo');
xlabel('Time');
ylabel('Amplitude');
legend('show');
grid on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PUNTO OPZIONALE 3 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Esplorare il range di ampiezza di riferimenti a gradino tali per cui il controllore rimane efficace sul
%% sistema non lineare.

% Parametri della simulazione
simTime = 0.5;  % Tempo di simulazione per ciascuna condizione iniziale
numSteps = 150;  % Numero di passi

figure(14);
hold on;

for i = -1110:numSteps:1000
    % Aumenta gradualmente il valore di WW
    x1_value =xe1;
    x2_value=xe2;
    WW_opz=i;

    % Simulazione con il valore di WW corrente
    simOut = sim('simulink_puntiOpzionali.slx', 'StopTime', num2str(simTime), 'LoadInitialState', 'off', 'SaveOutput', 'on', 'OutputSaveName', 'yout');
    outputData = simOut.y_nlin;

    % Plot del grafico per il valore di WW corrente
    plot(outputData.time, outputData.Data(:, 1), 'DisplayName', sprintf('WW = %d', WW_opz));
end

title('Risposte al variare del riferimento');
xlabel('Time');
ylabel('WW');
legend('show');
grid on;
hold off;

open("simulink_progetto.slx") %apre il progetto simulink
open("simulink_puntiOpzionali.slx") %apre il progetto simulink per i punti opzionali

