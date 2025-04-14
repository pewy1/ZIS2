

clc;
close all;
clear all;
%inicializace parametrů systému
R1=10;
R2=20;
C=0.001;
Ts=0.01;


%g1
T=((R1*R2)/(R1+R2))*C;
K=R2/(R1+R2);

%g2
a1k=-exp(-Ts/T);
b1k=K*(1-exp(-(Ts/T)));

theta_k=[a1k;b1k];
p=[R1;R2;C];

numerator = K;
denominator = [T 1];
G = tf(numerator, denominator); %přenos systému

Gd = c2d(G, Ts, 'zoh'); %diskrétní přenos


figure;
step(G);
grid on;
title('Spojitý přenos systému - odezva na jednotkový skok');
xlabel('Čas [s]');
ylabel('Výstup');

figure;
step(Gd);
grid on;
title('Diskretizovaný přenos c2d - odezva na jednotkový skok');
xlabel('Čas [s]');
ylabel('Výstup');



%% simulace - zadani (1) - inicializace parametrů


N = 15;             % počet kroků simulace
Ts = 0.01;           % vzorkovací perioda
u = ones(1, N);      % vstup (např. skoková funkce)
y = zeros(1, N);     % výstup
%y(1)=1;



%% simulace -zadani (1)- ARX modelu bez šumu



for k = 2:N

    %u(k)=10*sign(sin(2*pi*0.05*Ts*(k-1)));

    y(k) = -a1k * y(k-1) + b1k * u(k-1);

end

% Zobrazení výsledků
figure;
time = (0:N-1) * Ts;
stairs(time, y);
xlabel('čas [s]');
ylabel('y[k]');
title('Simulace diskrétního ARX modelu - jednotkový skok');
grid on;


%% simulace ARX modelu s šumem

whiteNoise = sqrt(0.02) * randn(N, 1);


for k = 2:N

    y(k) = -a1k * y(k-1) + b1k * u(k-1)+whiteNoise(k);

end

% Zobrazení výsledků
figure;
time = (0:N-1) * Ts;
stairs(time, y);
xlabel('čas [s]');
ylabel('y[k]');
title('Simulace diskrétního ARX modelu se šumem - jednotkový skok');
grid on;

figure;
plot(whiteNoise, 'ro', 'LineStyle','none','MarkerFaceColor','r');
hold on;
yline(0, 'k-', 'LineWidth', 1.5);  % Vodorovná čára na y=0
title('Bílý šum');




%% nejmenší čtverce



R1_promenne=0;
F = 7000;      
a_true=a1k;
b_true=b1k;
sigma = 0.01;           % standardní odchylka šumu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 0.99;          % zapomínací faktor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ts=0.01;

% Příprava paměti pro u, y, atd.
u = zeros(1,F);
y = zeros(1,F);
y_bezSumu=zeros(1,F);
y_est = zeros(1,F);
e = sqrt(0.02) * randn(F, 1); 
whiteNoise2= -10 + (20) * rand(F, 1);

R_matice=zeros(1,F);
thetaTrue = [a_true; b_true]; % "skutečné" parametry

% Odhad: 2 parametry => theta = [a; b]
thetaEstNoForget = zeros(2,F);   % RLS bez zapomínání
thetaEstForget   = zeros(2,F);   % RLS s expon. zapomínáním

T_real_simul=zeros(1,F); %reálné T při simulaci s proměnným R1;
K_real_simul=zeros(1,F); %reálné K při simulaci s proměnným R1;

T_est=zeros(1,F);   %hodnoty získané z odhadnutých parametrů
K_est=zeros(1,F);

err_ef=zeros(1,F);


% Kovarianční matice (pro klasické RLS a RLS s EF)
P_no = 100 * eye(2);    % počáteční velká nejistota
P_ef = 100 * eye(2);

% Pro uložení chyb odhadu
thetaErrNoForget = zeros(2,F);
thetaErrForget   = zeros(2,F);

% Počáteční podmínky pro y(k)
% (Je-li model 1. řádu, stačí y(1)=0, u(1)=0 atp.)
y(1) = 0;
u(1) = 1;  % třeba začneme s jednotkovým vstupem




% Simulace a odhad
for k = 2:F



    %změna parametru R1:
     if k < 3000
         R1_promenne = 10;

     elseif k < 5000
         R1_promenne = 10 + 0.1 * ((k * Ts) - 30);

     else
        R1_promenne= 12;

     end

    R_matice(k)=R1_promenne;


    T_real_simul(k)=((R1_promenne*R2)/(R1_promenne+R2))*C;
    K_real_simul(k)=R2/(R1_promenne+R2);

    T=((R1_promenne*R2)/(R1_promenne+R2))*C;
    K=R2/(R1_promenne+R2);
    a1k=-exp(-Ts/T);
    b1k=K*(1-exp(-(Ts/T)));
    a_true=a1k;
    b_true=b1k;





     u(k)=whiteNoise2(k); %šum jako vstupní signál
    
    %u(k)=10*sign(sin(2*pi*0.05*Ts*(k-1))); %obdelníkový vstupní signál


    % 2) Výpočet výstupu systému (skutečný model)
    y(k) = -a1k * y(k-1) + b1k * u(k-1)+e(k);
    y_bezSumu(k) = -a1k * y_bezSumu(k-1) + b1k * u(k-1);

    
    phi_k = [-y(k-1); u(k-1)];  % Regresor (2×1)


    % 4) RLS s exponenciálním zapomínáním
    denom_ef = lambda*sigma^2 + phi_k' * P_ef * phi_k;
    K_ef = (P_ef * phi_k) / denom_ef;
    % aktualizace odhadu
    err_ef(k) = y(k) - phi_k' * thetaEstForget(:,k-1); %rezidua
    thetaEstForget(:,k) = thetaEstForget(:,k-1) + K_ef * err_ef(k);
    % aktualizace kovarianční matice
    P_ef = (1/lambda) * (P_ef - K_ef * phi_k' * P_ef);
    % uložení chyby odhadu
    thetaErrForget(:,k) = thetaTrue - thetaEstForget(:,k); %chyba odhadu

     y_est(k) = -thetaEstForget(1,k)*y(k-1) + thetaEstForget(2,k)*u(k-1);

     a_est=thetaEstForget(1,k);
     b_est=thetaEstForget(2,k);


     T_est(k)=-((Ts)/log(-a_est));
     K_est(k)=b_est/(1-exp(-Ts/T_est(k)));


     

end

% Zobrazení výsledků
t = 1:F;
grid on;
figure;
plot(t, u, t, y, 'LineWidth', 1.2);
legend('Skutečný výstup se šumem','Vstup u(k)');
xlabel('k');
title('Skutečný výstup se šumem a vstup');
grid on;

figure;
plot(t, err_ef ,'LineWidth', 1.2);
legend('Reziduum r(k)');
xlabel('k');
title('Reziduum r(k)');
grid on;



figure;
plot(t, thetaEstForget(1,:), t, thetaEstForget(2,:), '-','LineWidth',1.2);
legend('Odhad a1','Odhad b1');
xlabel('k');
title('Odhadnuté parametry - LS s exponenciálním zapomínáním');
grid on;



figure;
hold on;
plot(t, thetaEstForget(1,:), 'LineWidth', 1.5);                 % a (odhadnutý)
plot(t, -exp(-Ts ./ T_real_simul), 'LineWidth', 1.5); hold on; % a (reálný)

legend('a1 reálný','a1 odhadnutý');
xlabel('k');
ylabel('Hodnota a');
title('Reálný a odhadnutý parametr a1');
grid on;


figure;
hold on;
plot(t, thetaEstForget(2,:), '-', 'LineWidth', 1.5);                                      % b (odhadnutý)
plot(t, K_real_simul .* (1 - exp(-Ts ./ T_real_simul)), '-', 'LineWidth', 1.5); hold on; % b (reálný)

legend('b1 reálný','b1 odhadnutý');
xlabel('k');
ylabel('Hodnota b1');
title('Reálný a odhadnutý parametr b1');
grid on;




figure;
plot(t, T_est(1,:), t, T_real_simul(1,:), '-','LineWidth',1.2);
legend('T odhadnutý','T reálný');
xlabel('k');
title('Reálný a odhadnutý parametr T');
grid on;

figure;
plot(t, K_est(1,:), t, K_real_simul(1,:), '-','LineWidth',1.2);
legend('K odhadnutý','K reálný');
xlabel('k');
title('Reálný a odhadnutý parametr K');
grid on;


figure;
plot(t, thetaErrForget(1,:), t, thetaErrForget(2,:),'LineWidth',1.2);
legend('Chyba a1','Chyba b1');
xlabel('k');
title('Chyba odhadu parametrů - LS s exponenciálním zapomínáním');
grid on;

figure;
plot(t,y,t, y_est, 'LineWidth', 1.2);
legend('Skutečný výstup se šumem','Odhadovaný výstup');
xlabel('k');
title('Skutečný výstup se šumem a odhadovaný výstup');
grid on;

figure;
plot(t, y_bezSumu, 'LineWidth', 1.2);
legend('Výstup bez šumu');
xlabel('k');
title('Výstup bez šumu');
grid on;


figure;
plot(t, R_matice, 'LineWidth', 1.2);
legend('R1');
xlabel('k');
title('Hodnota parametru R1 v čase');
grid on;




%% Popis proměnných ve skriptu

% --- Parametry systému ---
% R1, R2 .......... hodnoty odporů v ohmech
% C ............... kapacita v faradech
% Ts .............. vzorkovací doba (sampling period)

% --- Přenos systému ---
% T ............... časová konstanta systému (odvozená z R1, R2 a C)
% K ............... statické zesílení systému
% numerator ....... čitatel přenosové funkce (pro spojitý model)
% denominator ..... jmenovatel přenosové funkce (pro spojitý model)
% G ............... spojitý přenos systému
% Gd .............. diskretizovaný přenos systému pomocí ZOH

% --- Parametry ARX modelu ---
% a1k ............. diskrétní parametr vyjadřující vliv minulého výstupu
% b1k ............. diskrétní parametr vyjadřující vliv minulého vstupu
% theta_k ......... vektor nominálních (počátečních) parametrů [a1k; b1k]
% p ............... vektor fyzikálních parametrů [R1; R2; C]

% --- Simulace ARX bez a se šumem ---
% u ............... vstupní signál (např. skok nebo šum)
% y ............... výstup systému
% y_bezSumu ....... výstup systému bez přidaného šumu
% whiteNoise ...... generovaný bílý šum (přidávaný do výstupu)
% whiteNoise2 ..... náhodný vstupní signál (šum pro u)

% --- Adaptivní odhad parametrů (RLS) ---
% F ............... počet kroků simulace (časových vzorků)
% e ............... šum měření (rušení na výstupu)
% thetaTrue ....... skutečné (aktuální) hodnoty parametrů v čase [a; b]
% thetaEstForget .. odhadnuté parametry pomocí RLS s exponenciálním zapomínáním
% thetaErrForget .. chyba odhadu parametrů (thetaTrue - thetaEstForget)
% y_est ........... odhadnutý výstup systému na základě thetaEstForget
% phi_k ........... regrese [−y(k−1); u(k−1)]
% P_ef ............ kovarianční matice odhadu pro RLS-EF
% K_ef ............ zisk (gain) RLS algoritmu
% err_ef .......... reziduum – rozdíl mezi skutečným a odhadnutým výstupem
% lambda .......... zapomínací faktor RLS algoritmu (blízko 1 = pomalé zapomínání)
% sigma ........... odhad směrodatné odchylky šumu (pro RLS)
% R1_promenne ..... časově proměnná hodnota R1 (simulace změn systému)
% R_matice ........ historie hodnot R1 v čase (pro vykreslení)

% --- Grafy ---
% time ............ časová osa pro diskrétní simulace
% t ............... index pro grafické zobrazení vektorů (1:F)

% --- Volitelné ---
% thetaEstNoForget .... připraveno pro RLS bez zapomínání (nepoužívá se)
% thetaErrNoForget .... chyba odhadu bez zapomínání (nepoužívá se)
