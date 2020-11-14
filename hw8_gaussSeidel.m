function hw8_gaussSeidel(n)

% for plotting the unit circle in pole zero plots: 
r = 1;
th = 0:pi/50:2*pi;
xunit = r * cos(th);
yunit = r * sin(th);

% get diagonally dominant matrix:
randMat = diagonallyDominantMaker(n);

% chose some b matrix (doesn't matter what it starts from)
b = randn(n,1);

% calculated the Atildat and Btilda for the Jacobi algorithm:
A = randMat;

% first step for Gauss Seidal = decompose the matrix into 
% lower and upper matrices: 
L = tril(A);
U = A - L;

% calculate the Jacobi tran
Atil = -inv(L)*U;
Btil = inv(L);

numIter = 100;

x = randn(n,1);

% define the number of iterations it goes through:
for cnt1 = 1:(numIter - 1)
    x(:,cnt1+1) = Atil*x(:,cnt1) + Btil*b;
    err(cnt1) = norm(x(:,cnt1) - inv(Atil)*b);
end

str = sprintf('Gauss Seidel method \n for diagonal dominant n = %g',n);

figure,
hold on
for i = 1:n
    plot(real(x(i,:)))
end
title(str)
xlabel('iteration')
ylabel('error')
hold off

str2 = sprintf('Gauss Seidel method \n error for diagonal dominant n = %g',n);

figure,
plot(err)
title(str2)
xlabel('iteration')
ylabel('error')

% (iv): relaxation of constraints:
alpha = [0.1,0.5,0.9];

for i = 1:3
    for cnt1 = 1:(numIter - 1)
        x(:,cnt1+1) = (1-alpha(i))*(x(:,cnt1)) + alpha(i)*(Atil*x(:,cnt1) + Btil*b);
        err(cnt1) = norm(x(:,cnt1) - inv(A)*b);
    end
    
    % some string manipulation:
    str3 = sprintf('Gauss Seidel method \n method for diagonally dominant relaxed w/ alpha = %g, n = %g', alpha(i),n);
    
    figure,
    hold on
    for j = 1:n
        plot(real(x(j,:)))
    end
    title(str3)
    xlabel('iteration')
    ylabel('value')
    hold off
    
    str4 = sprintf('Gauss Seidel\n method error for diagonally \n dominant relaxed state system with alpha = %g, n = %g', alpha(i),n);

    figure,
    plot(err)
    title(str4)
    xlabel('iteration')
    ylabel('error')

end

%% Plotting eigen values:
A = randMat;

eigenVals = eig(A);

limitAxis = max(eigenVals);

figure,
hold on
plot(eigenVals,'o')
h = plot(xunit, yunit);
xlim([-limitAxis,limitAxis])
ylim([-limitAxis, limitAxis])
xlabel('Real')
ylabel('Imaginary')
title('pole zero plot for Jacobi with\n diagonally dominant matrix')
hold off


%% generate symmetric positive definite matrices:

% easiest example: 0 matrix with diagonal as all positive values:
zeromat = zeros(n,n);

% add real values numbers: 
ranMat2 = rand(n);

% dot multiply to get only values on diag:
symmPosDef = ranMat2.*(eye(n));

b = randn(n,1);

% first step for Gauss Seidal = decompose the matrix into 
% lower and upper matrices: 
L = tril(A);
U = A - L;

% calculate the Jacobi tran
Atil = -inv(L)*U;
Btil = inv(L);

numIter = 100;

x = randn(n,1);
% define the number of iterations it goes through:
for cnt1 = 1:(numIter - 1)
    x(:,cnt1+1) = Atil*x(:,cnt1) + inv(L)*b;
    err(cnt1) = norm(x(:,cnt1) - inv(A)*b);
end

str = sprintf('Gauss Seidel method \n for diagonal dominant n = %g',n);

figure,
hold on
for i = 1:n
    plot(real(x(i,:)))
end
xlabel('iteration')
ylabel('value')
title(str)
hold off

str5 = sprintf('Gauss Seidel \n residual error for size n = %g',n);

figure,
plot(err)
title(str5)
xlabel('iteration')
ylabel('error')

% (iv): relaxation of constraints:

alpha = [0.1,0.5,0.9];

for i = 1:length(alpha)
    clear str1

    for cnt1 = 1:(numIter - 1)
        x(:,cnt1+1) = (1-alpha(i))*(x(:,cnt1)) + alpha(i)*(Atil*x(:,cnt1) + inv(L)*b);
        err(cnt1) = norm(x(:,cnt1) - inv(A)*b);
    end
    
    % some string manipulation:
    str1 = sprintf('Gauss Seidel \n relaxed state system with alpha = %g', alpha(i));
    
    figure,
    hold on
    for l = 1:n
        plot(real(x(l,:)))
    end
    title(str1)
    xlabel('iteration')
    ylabel('value')
    hold off

    str5 = sprintf('Gauss Seidel \n residual error relaxed state system with alpha = %g', alpha(i));

    figure,
    plot(err)
    title(str5)
    xlabel('iteration')
    ylabel('residual error')
    
end

%%
A = symmPosDef;

eigenVals = eig(A);
limitAxis = max(max(eigenVals(:,:)));
figure,
hold on
plot(eigenVals,'o')
h = plot(xunit, yunit);
%xlim([-limitAxis,limitAxis])
%ylim([-limitAxis, limitAxis])
xlabel('Real')
ylabel('Imaginary')
title('pole zero plot for Jacobi with\n positive semi-denite matrix')
hold off

nothing = 0;
end
