% Ejercicios de Monte Carlo - Métodos Numéricos Estocásticos


%% Ejercicio 1

% Generar una uniforme U entre a y b con el comando rand de Matlab para un
% número N de realizaciones.
clear all; clc

a = -10; b = 10; N = 10;
X = rand(N,1);
Y = ( b - a )*X - a;
hist(Y); pause


% Generar, con una U(0,1), una distribución de variable discreta que cumpla
% unas determinadas condiciones y probabilidades.
hist(DiscreteRand(rand(1,10000))); pause


% Generar una función PDF exponencial a partir de la U(0,1), usando para
% ello el teorema de la función inversa.
clear all; clc

a = 0; b = 1; N = 200; k = 2;
X = a + (b - a)*rand(N,1);
Y = (-1/k) * log(X);

hold on; histogram( Y,'Normalization','pdf')
y = [0:0.001:10]; plot(y,k*exp(-k*y),'-r'); hold off; pause


% Box-Muller (Elapsed time is 0.084815 seconds.)
clear all; clc

a = 0; b = 1; N = 10000;
X = a + (b - a)*rand(N,2);
Y = [sqrt(-2*log(X(:,1))).*cos(2*pi*X(:,2)),...
     sqrt(-2*log(X(:,1))).*sin(2*pi*X(:,2))];

figure(1); hold on; histogram( Y(:,1),'Normalization','pdf')
y = [-10:0.001:10]; plot(y,exp(-0.5.*y.*y)./sqrt(2*pi),'-r'); hold off
figure(2); hold on; histogram( Y(:,2),'Normalization','pdf')
y = [-10:0.001:10]; plot(y,exp(-0.5.*y.*y)./sqrt(2*pi),'-r'); hold off
pause


% Marsaglia (Elapsed time is 0.007269 seconds.)
clear all; clc

a = -1; b = 1; N = 1.27*10000; %(almost 1000)
X = a + (b - a)*rand(N,2); X = X(find(X(:,1).^2 + X(:,2).^2 < 1),:);
W = X(:,1).^2 + X(:,2).^2;
Y = [X(:,1).*sqrt(-2.*log(W)./W), X(:,2).*sqrt(-2.*log(W)./W)];

figure(1); hold on; histogram( Y(:,1),'Normalization','pdf')
y = [-10:0.001:10]; plot(y,exp(-0.5.*y.*y)./sqrt(2*pi),'-r'); hold off
figure(2); hold on; histogram( Y(:,2),'Normalization','pdf')
y = [-10:0.001:10]; plot(y,exp(-0.5.*y.*y)./sqrt(2*pi),'-r'); hold off
pause



%% Ejercicio 2 

% Generar una uniforme N entre a y b con el comando randn de Matlab para un
% número N de realizaciones.
clear all; clc

M = [1 2]; S = [1 0.5; 0.5 2]; N = 1e4; D = length(M);
M = repmat(M,N,1); S = chol(S);

X = randn(N,D);
Y = M + X*S;


%% Ejercicio 3






%% Ejercicio 4

% Calcular la integral de exp(x) desde 0 hasta 1 usando el rand de Matlab y
% em Método Monte Carlo.
clear all; clc

a = 0; b= 1; N = 1e4; cl = 0.95;                         % confidence level
f = @(x) exp(x);
X = a + (b-a)*rand(N,1); Y = exp(X);
fprintf('I = %.8f +/- %.8f\n',mean(Y),tinv(cl,N-1)*std(Y)/sqrt(N))



% Calcular la integral anterior usando variables antitéticas.

N = N/2;                                   % because 2 samples will be done
X = a + (b-a)*rand(N,1); X = [X,1-X];
Y = f(X); C = cov(Y(:,1),Y(:,2));
SY = 0.5*(var(Y(:,1)) +  C(1,2));                         % Variance sample
sY = sqrt(SY)/sqrt(N);                                 % Std sample/sqrt(N)
fprintf('I = %.8f +/- %.8f\n',0.5*mean(Y(:,1)+Y(:,2)),tinv(cl,N-1)*sY)
pause


%% Ejercicio 5

% Calcular F(y) a partir de f(y) y también la función h(x) que define una
% muestra de puntos tales que siguen la función F(u).
clear all; clc

f = @(x) (x<-15)*(0)  + ((x>=-15)&(x<0))*(1/60)           + ((x>=0)&(x<15))*(1/30)            + ((x>=15)&(x<30))*(1/60)                         + ((x>=30)&(x<0))*(0)
F = @(x) ((x<-15)*(0) + ((x>=-15)&(x<0)).*(0 + (x+15)/60) + ((x>=0)&(x<15)).*(1/4             + (x)/30)   + ((x>=15)&(x<30)).*(3/4 + (x-15)/60) + ((x>=30))*(1))
h = @(x)                ((x>=0)&(x<1/4)).*(60*x -15)      + ((x>=1/4)&(x<3/4)).*((60*x-15)/2) + ((x>=3/4)&(x<1)).*(60*x -30) 

% Obtener una muestra de la función del apartado anterior.
a=0; b=1; N =1e4;
X = a + (b-a)*rand(N,1);
histogram( h(X),'Normalization','pdf'); hold on
plot([-30:0.01:45],f([-30:0.01:45]))

% Obtener media para una muestra como la anterior de 10, 100 y 1000
% realizaciones.
sample10   = h(rand(10,1));    mean10   = sum(sample10)/10;
sample100  = h(rand(100,1));   mean100  = sum(sample100)/100;
sample1000 = h(rand(1000,1));  mean1000 = sum(sample1000)/1000;

sigma210   = sum( (sample10 - mean10).^2 )/(10-1);
sigma2100  = sum( (sample100 - mean100).^2 )/(100-1);
sigma21000 = sum( (sample1000 - mean1000).^2 )/(1000-1);

sigma10mean = sqrt(sigma210/10);
sigma100mean = sqrt(sigma2100/100);
sigma1000mean = sqrt(sigma21000/1000);

cl = 0.95
fprintf('  10 items - %.8f +/- %.8f\n',mean10,tinv(cl,N-1)*sigma10mean)
fprintf(' 100 items - %.8f +/- %.8f\n',mean100,tinv(cl,N-1)*sigma100mean)
fprintf('1000 items - %.8f +/- %.8f\n',mean1000,tinv(cl,N-1)*sigma1000mean)











%% Functions

function Y = DiscreteRand(X)
    for l = 1:1:length(X)
        if     0.0 < X(l) < 0.1; Y(l) = 1;
        elseif 0.1 < X(l) < 0.3; Y(l) = 3;
        elseif 0.3 < X(l) < 0.7; Y(l) = 5;
        elseif 0.7 < X(l) < 1.0; Y(l) = 9;
        end
    end
end



