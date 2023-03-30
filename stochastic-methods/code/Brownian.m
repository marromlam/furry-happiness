


function [ S, t ] = Brownian(t0,N)

    h = 1/N;
    X = randn(N,1);
    S = zeros(N,1); t = zeros(N,1);
    S(1) = 0; t(1) = t0;
    for i = 1:1:N-1; 
        S(i+1) = S(i)+sqrt(h)*X(i);
        t(i+1) = t(i)+h;
    end

end