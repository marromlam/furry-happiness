clear all; clc

%% Modelo de población biológico

% Resolver numericamente con Euler-Maruyama el modelo.
  a = @(t,Y) 2*t*Y + exp(t+t*t);
  b = @(t,Y) 3*exp(t*t);
  T0 = 0; Y0 = 1; N = 100; T = 1;
  X = Brownian(T0,N);
  [t, Yh] = EulerMaruyama( a, b, T, N, Y0, T0, X );

% Solución analítica
  Y = Y0*ones(size(X)); int = zeros(size(X));
  for l = 2:1:N
    int(l) = int(l-1) + 3*(X(l)-X(l-1));
  end
  Y = exp(t.*t).*(Y0 + (exp(t)-1) + int);
  figure(1); hold on; plot(t,Yh,'r*-'); plot(t,Y); hold off

% Tabla del ejercicio
  matAux1 = zeros(6,3); 
  N = [10,100]; n = [10,100,1000]; 
  for j = 1:1:length(N)
    for k = 1:1:length(n)
      Naux = N(j); naux = n(k);
      for l = 1:1:naux
        vecAux1 = zeros(naux,1); vecAux2 = zeros(naux,1);
        X = Brownian(T0,Naux);
        [t, Yh] = EulerMaruyama( a, b, T, Naux, Y0, T0, X );
        Y = Y0*ones(size(X)); int = zeros(size(X));
        for m = 2:1:Naux
          int(m) = int(m-1) + 3*(X(m)-X(m-1));
        end
        Y = exp(t.*t).*(Y0 + (exp(t)-1) + int);
        vecAux2(l) = Yh(Naux); vecAux1(l) = Y(Naux);
      end
      fprintf('$%i$ & $%i$ & ',Naux,naux);
      fprintf('$%.4f$ & $%.4f$ & ',mean(vecAux1),mean(vecAux2));
      fprintf('$%.4f$ \\\\ \\hline \n',mean(abs(vecAux1-vecAux2)))
    end
  end
  
%% El modelo de tipos de interes a corto propuesto por CIR
clear all; clc; close all
% Parámetros
  a = @(t,Y) 0.3*(0.04-Y);
  b = @(t,Y) 0.2*sqrt(Y);
  T0 = 0; Y0 = 0; N = 1000; T = 1;

% Tabla del ejercicio
  matAux1 = zeros(6,1); 
  N = [10,100]; n = [10,100,1000]; 
  for j = 1:1:length(N)
    for k = 1:1:length(n)
        figure(k); hold on
      Naux = N(j); naux = n(k);
      for l = 1:1:naux
        vecAux1 = zeros(naux,1); vecAux2 = zeros(naux,1);
        X = Brownian(T0,Naux);
        [t, Yh1] = EulerMaruyama( a, b, T, Naux, Y0, T0, X );
        vecAux1(l) = Yh1(Naux);
        [t, Yh2] = MilsteinI( a, b, T, Naux, Y0, T0, X ); 
        vecAux2(l) = Yh2(Naux);
      end
       plot(t,real(Yh1));  plot(t,real(Yh2)); hold off
      fprintf('$%i$ & $%i$ & $%.6f$ &',Naux,naux,mean(vecAux1));
      fprintf(' $%.6f$ \\\\ \\hline\n',mean(vecAux2));
    end
  end
