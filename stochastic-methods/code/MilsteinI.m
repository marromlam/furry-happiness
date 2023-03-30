function [ t, Y ] = MilsteinI( a, b, T, N, Y0, T0, X )

  t = linspace(T0,T,N)';
  h = (T-T0)/N;
  Y = Y0*ones(size(X));

  for l = 1:1:N-1
    Y(l+1) = Y(l) + a(t(l),Y(l))*h + b(t(l),Y(l))*(X(l+1)-X(l)) + ...
    0.5*b(t(l),Y(l))*(NumJac(@(w) b(t(l),w),Y(l)))*((X(l+1)-X(l))^2-h);
  end

end

% Cálculo de derivada numérica de una función (matriz jacobiana)

function [ J ] = NumJac( f, x0 )

  [r,c] = size(f(zeros(length(x0))));
  c = max(r,c); r = max(r,length(x0));
  J = zeros(r,c); 
  for l = 1:1:r
    if x0(l) ~= 0
      h = sqrt(eps)*x0(l);
    else 
      h = 1e-14;
    end
    xh1 = x0; xh1(l) = x0(l) + h; 
    xh2 = x0; xh2(l) = x0(l) - h;
    J(l,:)  = (feval(f,xh1) - feval(f,xh2))./(2*h);
  end
  J = transpose(J);
  
end

