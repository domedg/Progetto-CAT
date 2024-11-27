%% Progetto Controlli Automatici T B1
% Barone Leonardo, Del Giudice Domenico, Galli Francesco, Guzzonato Leonardo

clear all; close all; clc;
% Paramentri Tecnologici
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

A=[r_s*(1-(2*xe1+xe2)/K)-m_s*u_e-beta-alfa*u_e, (r_s*xe1)/K+gamma;
    (r_r*xe2)/K+beta+alfa*u_e, r_r*(1-(xe1+2*xe2)/K)-m_r*u_e-gamma];
A
B=[-m_s*xe1-alfa*xe1; -m_r*xe2+alfa*xe1];
B
C=[0, 1];
C
D=0;
D
%% Punto2
%funzione di trasferimento
% s=tf('s');
% GG1=C*(inv(s*eye(2)-A))*B-D;
modello = ss(A,B,C,D);
GG = tf(modello);

GG
p=pole(GG);
z=zero(GG);
p
z

bode(GG);

%% Visualizzazione del diagramma di Bode
% Creazione della figura
figure('Name', 'Diagramma di Bode', 'NumberTitle', 'off', 'Color', 'w');

% Disegna il diagramma di Bode
[mag, phase, wout] = bode(GG);

% Personalizzazione
subplot(2, 1, 1); % Magnitudine
semilogx(wout, 20*log10(squeeze(mag)), 'LineWidth', 1.5, 'Color', '#0072BD'); % Linea blu
grid on;
title('Diagramma di Bode - Magnitudine', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequenza (\omega) [rad/s]', 'FontSize', 12);
ylabel('Ampiezza [dB]', 'FontSize', 12);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);

subplot(2, 1, 2); % Fase
semilogx(wout, squeeze(phase), 'LineWidth', 1.5, 'Color', '#D95319'); % Linea arancione
grid on;
title('Diagramma di Bode - Fase', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequenza (\omega) [rad/s]', 'FontSize', 12);
ylabel('Fase [°]', 'FontSize', 12);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);

% Personalizzazione generale
sgtitle('Diagramma di Bode della Funzione di Trasferimento', 'FontSize', 16, 'FontWeight', 'bold');









% % Confronta le risposte in frequenza
% % Confronta le risposte in frequenza con linee piene e tratteggiate
% figure;
% bodeplot(GG, 'b', GG1, 'r--'); % 'b' per blu (linea continua) e 'r--' per rosso (linea tratteggiata)
% 
% % Aggiungi la legenda
% legend('GG', 'GG1');
% 
% % Aggiungi il titolo
% title('Confronto tra le risposte in frequenza');
% grid on;

