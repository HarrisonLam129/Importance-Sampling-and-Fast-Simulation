n = 10000;
lambda = 1/5;
%Generate n Unif[0,1] samples
u = rand(1,n);
%Transform into n Exp(lambda) samples
explambda = -(1/lambda)*log(u);

gammahat = zeros(1,n);
for i = 1:n
    if explambda(i) > 30
        gammahat(i) = (1/(3*lambda))*exp((lambda-1/3)*explambda(i));
    end
end
estimate = mean(gammahat);


N = 10000;
lambda = 1/50;
result = zeros(1,N);

for n = 1:N
    u = rand(1,n);
    explambda = -(1/lambda)*log(u);

    gammahat = zeros(1,n);
    for i = 1:n
        if explambda(i) > 30
            gammahat(i) = (1/(3*lambda))*exp((lambda-1/3)*explambda(i));
        end
    end
    result(i) = mean(gammahat);
end
plot(result)
title(append('$\lambda$ = 1/', num2str(1/lambda)), interpreter = 'latex');
xlabel('n');
ylabel('Estimate')

%Fix N
N = 5200;
rep = 10000;
resultmatrix = zeros(1,rep);
for j = 1:rep
    u = rand(1,N);
    explambda = -(1/lambda)*log(u);

    gammahat = zeros(1,N);
    for i = 1:N
        if explambda(i) > 30
            gammahat(i) = (1/(3*lambda))*exp((lambda-1/3)*explambda(i));
        end
    end
    resultmatrix(j) = mean(gammahat);
end
%Estimate the probability of the result being 'good' for fixed N 
frac = sum(resultmatrix<=5e-5 & resultmatrix>=4.09e-5)/rep


%Question 8
p1 = 1/4;
p2 = 0.75;
C = 30;
N = 10000;
rep = 10000;
resultmatrix = zeros(1,rep);

for j = 1:rep
    total = 0;
    for i = 1:N
        x = 1;
        Lratio = 1;
        while x ~= 0 && x ~= C
            if rand < p2
                x = x+1;
                Lratio = Lratio*p1/p2;
            else
                x = x-1;
                Lratio = Lratio*(1-p1)/(1-p2);
            end
        end
    
        if x == C
            total = total + Lratio;
        end
    end  
    resultmatrix(j) = total/N;
end

frac = sum(resultmatrix<=1.02e-14 & resultmatrix>=9.228e-15)/rep;

truevalue = 9.714e-15;
varestimate = mean((resultmatrix-truevalue).^2)
