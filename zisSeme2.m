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
G = tf(numerator, denominator); %přenos 

Gd = c2d(G, Ts, 'zoh'); %diskrétní přenos


figure;
step(G);
grid on;
title('Odezva na jednotkový skok - přenos');
xlabel('Čas [s]');
ylabel('Výstup');

figure;
step(Gd);
grid on;
title('Odezva na jednotkový skok - přenos c2d');
xlabel('Čas [s]');
ylabel('Výstup');

%% simulace - inicializace parametrů


N = 20;             % počet kroků simulace
Ts = 0.01;           % vzorkovací perioda
u = ones(1, N);      % vstup (např. skoková funkce)
y = zeros(1, N);     % výstup
%y(1)=1;



%% simulace ARX modelu bez šumu


for k = 2:N
    y(k) = -a1k * y(k-1) + b1k * u(k-1);

end

% Zobrazení výsledků
figure;
time = (0:N-1) * Ts;
stairs(time, y);
xlabel('čas [s]');
ylabel('y[k]');
title('Simulace diskrétního ARX modelu');
grid on;


%% simulace ARX modelu s šumem



%noise = zeros(1, N);
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
title('Simulace diskrétního ARX modelu se šumem');
grid on;

figure;
plot(whiteNoise, 'ro', 'LineStyle','none','MarkerFaceColor','r');
hold on;
yline(0, 'k-', 'LineWidth', 1.5);  % Vodorovná čára na y=0
title('Bílý šum');




