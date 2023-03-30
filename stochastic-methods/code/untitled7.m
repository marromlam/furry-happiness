


for l = 1:1:1e2
    hold on
    [X,Y] = Brownian(0,250);
    plot(X, Y)
    csvwrite('Brow' + string(l) + '.csv',[X,Y])
end