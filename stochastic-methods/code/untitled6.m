
X = zeros(1e6,1);
f = @(x) (x<-15)*(0)  + ((x>=-15)&(x<0))*(1/60)           + ((x>=0)&(x<15))*(1/30)            + ((x>=15)&(x<30))*(1/60)                         + ((x>=30)&(x<0))*(0)


tic
for l = 1:1:1e6
    X(l) = AnyRand(f,-30,45);
end
toc
%histogram(X)