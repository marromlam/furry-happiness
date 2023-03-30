function [ t, Y ] = EulerMaruyama( a, b, T, N, Y0, T0, X )

  t = linspace(T0,T,N)';
  h = (T-T0)/N;
  Y = Y0*ones(size(X));

  for l = 1:1:N-1
    Y(l+1) = Y(l) + a(t(l),Y(l))*h + b(t(l),Y(l))*(X(l+1)-X(l));
  end

end

