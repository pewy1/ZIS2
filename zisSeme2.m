
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
lambda = 0.9;          % zapomínací faktor
Ts=0.01;

% Příprava paměti pro u, y, atd.
u = zeros(1,F);
y = zeros(1,F);
y_bezSumu=zeros(1,F);
y_est = zeros(1,F);
e = sqrt(0.02) * randn(F, 1); 
whiteNoise2= -10 + (20) * rand(N, 1);

R_matice=zeros(1,F);
thetaTrue = [a_true; b_true]; % "skutečné" parametry

% Odhad: 2 parametry => theta = [a; b]
thetaEstNoForget = zeros(2,F);   % RLS bez zapomínání
thetaEstForget   = zeros(2,F);   % RLS s expon. zapomínáním


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

    T=((R1_promenne*R2)/(R1_promenne+R2))*C;
    K=R2/(R1_promenne+R2);
    a1k=-exp(-Ts/T);
    b1k=K*(1-exp(-(Ts/T)));
    a_true=a1k;
    b_true=b1k;





    % u(k)=whiteNoise2(k); %šum jako vstupní signál
    
    u(k)=10*sign(sin(2*pi*0.05*Ts*(k-1))); %obdelníkový vstupní signál


    % 2) Výpočet výstupu systému (skutečný model)
    y(k) = -a1k * y(k-1) + b1k * u(k-1)+e(k);
    y_bezSumu(k) = -a1k * y_bezSumu(k-1) + b1k * u(k-1);

    
    phi_k = [-y(k-1); u(k-1)];  % Regresor (2×1)


    % 4) RLS s exponenciálním zapomínáním
    denom_ef = lambda*sigma^2 + phi_k' * P_ef * phi_k;
    K_ef = (P_ef * phi_k) / denom_ef;
    % aktualizace odhadu
    err_ef = y(k) - phi_k' * thetaEstForget(:,k-1);
    thetaEstForget(:,k) = thetaEstForget(:,k-1) + K_ef * err_ef;
    % aktualizace kovarianční matice
    P_ef = (1/lambda) * (P_ef - K_ef * phi_k' * P_ef);
    % uložení chyby odhadu
    thetaErrForget(:,k) = thetaTrue - thetaEstForget(:,k);

     y_est(k) = -thetaEstForget(1,k)*y(k-1) + thetaEstForget(2,k)*u(k-1);
end

% Zobrazení výsledků
t = 1:F;

figure;
plot(t, y, t, u, 'LineWidth', 1.2);
legend('skutečný výstup se šumem','Vstup u(k)');
xlabel('k');
title('skutečný výstup se šumem a vstup');


figure;
plot(t, thetaEstForget(1,:), t, thetaEstForget(2,:), '--','LineWidth',1.2);
legend('a - RLS-EF','b - RLS-EF');
xlabel('k');
title('Parametry - RLS s exponenciálním zapomínáním');


figure;
plot(t, thetaErrForget(1,:), t, thetaErrForget(2,:),'LineWidth',1.2);
legend('Chyba a - EF','Chyba b - EF');
xlabel('k');
title('Chyba odhadu - s exponenciálním zapomínáním');

figure;
plot(t,y,t, y_est, 'LineWidth', 1.2);
legend('skutečný výstup se šumem','odhadovaný výstup');
xlabel('k');
title('skutečný (+šum) a odhadovaný výstup');
grid on;

figure;
plot(t, y_bezSumu, 'LineWidth', 1.2);
legend('výstup bez šumu');
xlabel('k');
title('výstup bez sumu');
grid on;


figure;
plot(t, R_matice, 'LineWidth', 1.2);
legend('R1');
xlabel('k');
title('R1 v čase');
grid on;
