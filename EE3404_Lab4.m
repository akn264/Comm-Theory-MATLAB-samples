%% EE3404 Lab 4: Random Variables | Aimee Nogoy | akn264

%% Q1: Generate a random Gaussian var
clear; clc; close all

xmean = 2.0; %E(X)
xvar = 0.5; %var(X)
nx = 1e4; %number of random samples
u = randn(nx,1); 
%lin transformation parameters
a = sqrt(xvar); b = 2;
x = a*u + b;  

%% Q2: Compute expectations
close all;

% i. E(X) 
Ex = (1/nx)*sum(x); %empirical
true_Ex = mean(x); %theoretical

% ii. var(X)
x2 = power(x,2);
Varx =  (1/nx)*(sum(x2 - Ex^2)); %empirical
true_var = var(x); %theoretical

% iii. E|X-E(X)|
huh = (1/nx)*sum(abs(x - Ex)); %empirical
true_huh = mean(abs(x-(true_Ex))); %theoretical

% iv. P(X > 3) = 1 - P(X < 3)
% ^^ for left part use Q(z) for right part use Phi(z)

z = (3-xmean)/sqrt(xvar);
Q = 0.5*(erfc(z/sqrt(2)));
%empirical:
num = 0;
for i=1:length(x)
    if x(i) > 3
        num = num+1;
    end
end
P3 = num/nx;
%theoretical:
true_P3 = Q;

%% Q3: Plotting a CDF
close all;

Xn = sort(x);
n = 1:nx;
y = n/nx;
zvec = (Xn-xmean)./sqrt(xvar); 
yt = 1-0.5*(erfc(zvec./sqrt(2)));
figure
plot(Xn,y,...
    'b',...
    'LineWidth',2)
hold on
plot(Xn,yt,...
    'g',...
    'LineWidth',2)
grid on
title('Empirical and Theoretical CDFs')
xlabel('Samples z'); ylabel('CDF')
legend('Empirical','Theoretical','Location','northwest')
set(gca,'Fontsize',11)
axis([0 4 0 1])

%they are almost the same!!!!

%% Q4: Generating uniform discrete random variables
close all
clear; clc;
%%
N = 1e3;
M = 4;
i = 1:M;
X = randi(M,[1,N]);

% theres def a shorter way to do this but i aint no natural coder
g = find(X==i(1)); h=find(X==i(2)); j = find(X==i(3));...
    k = find(X==i(M));
phatg = (1/N)*numel(g);
phath = (1/N)*numel(h);
phatj = (1/N)*numel(j);
phatk = (1/N)*numel(k);
true_p = 1/M;

%display for comparison
disp('The empirical probabilities are: ')
G = ['phatg = ' num2str(phatg)]; H = ['  phath = ' num2str(phath)];
J = ['  phatj = ' num2str(phatj)]; K = ['  phatk = ' num2str(phatk)];
S = [G H J K '  true_p = ' num2str(1/M)];
disp(S)
disp('All the estimates are approximately 0.25')

%% Q5: Generating general discrete random variables 
clear; clc; close all

n = 1000;
k = 0:8;
C = 1/(sum(exp(-k/2)));
p = C*exp(-k/2); %desired PMF
q = [0 cumsum(p)]; %CDF. note the zero added at the beginning
u = rand(n,1); %generate uniform random vars
[~,x] = histcounts(u,q,'Normalization',...
    'probability'); 
%^^count number of times u(n) is in each bin. normalized

xax = 0:length(q)-1; %x axis
figure
stem(xax,q,'g','LineWidth',2)
title('Theoretical and Empirical CDF')
grid on
hold on
stem(xax,x,'k','diamond')
legend('Theoretical','Empirical','Location','northwest')

%% Q6: Generating an exponential random variable
clear; clc; close all

n = 1e4;
u = rand(n,1);
u = sort(u);
x = -log(1-u);
x = sort(x);
Fx = 1-exp(-x);
figure
plot(x,Fx,'--','LineWidth',3) %theoretical CDF
grid on ; hold on
plot(x,(1:n)/n,'Linewidth',1.5) %empirical CDF
title('Theoretical and Empirical CDF')
legend('Theoretical','Empirical')





