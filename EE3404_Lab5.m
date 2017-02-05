%% EE3404 Lab 5: Power Detection | Aimee Nogoy | akn264

%% Q1: Generate random data
clear; clc; close all

p = 0.1; %probability of TX
la0 = 1; %lambda 0
n = 1e5;
gamma = 10; % signal-to-noise ratio
la1 = (1+gamma)*la0; %lambda 1

u = rand(n,1); 
x = (u < p);
laval = [la0; la1];
la = laval(x+1);
y = exprnd(la); %generate random samples of y
% ^ the expectation of an exponential rv is 1/lambda but the matlab
% function just requires the lambda


%% Q2: ML detection
format; clc

t = (((1/la0)-(1/la1))^-1)*log(la1/la0);
xhat = (y >= t);

% false alarm
b = sum((xhat == 1) & (x == 0)); 
a = sum(x == 0);
Pfa = b / a;

% missed detection
d = sum((x == 1) & (xhat == 0));
c = sum(x == 1);
Pmd = d / c;

% error rate
Perr = sum(xhat ~= x) / n;

disp('For ML detection:')
txt1=['Probability of false alarm: ' num2str(Pfa)];
txt2=['Probability of missed detection: ' num2str(Pmd)];
txt3=['Probability of error: ' num2str(Perr)];
disp(txt1)
disp(txt2)
disp(txt3)


%% Q3: MAP detection
clc

t = (((1/la0)-(1/la1))^-1)*((log(1/la0))-(log(1/la1))...
    +(log(0.9))-(log(0.1)));
xhat = (y>=t);

% false alarm
b = sum((xhat==1) & (x==0)); 
a = sum(x==0);
Pfa = b / a;

% missed detection
d = sum((xhat==0) & (x==1));
c = sum(x==1);
Pmd = d / c;

% error rate
Perr = sum(xhat ~= x) / n;

disp('For MAP detection:')
txt1=['Probability of false alarm: ' num2str(Pfa)];
txt2=['Probability of missed detection: ' num2str(Pmd)];
txt3=['Probability of error: ' num2str(Perr)];
disp(txt1)
disp(txt2)
disp(txt3)

%% Q4: ROC curve
clc; close all
tvec = linspace(min(y),max(y),500);
len = length(tvec);
X = x(1:len)';
Y = y(1:len)';
Xhat = zeros(1,len);
pfa = zeros(1,len); pmd = zeros(1,len);
% false alarm and missed detection
for k = 1:len
   Xhat = (Y >= tvec(k));
   num1 = sum(Xhat==1 & X==0);
   den1 = sum(X==0);
   num2 = sum(Xhat==0 & X==1);
   den2 = sum(X==1);
   pfa(k) = num1/den1;
   pmd(k) = num2/den2;
end
figure
plot(pmd,pfa,'LineWidth',1.5)
title('P_M_D versus P_F_A ROC Curve')
ylabel('P_M_D'); xlabel('P_F_A')
grid on

shiz = round(pfa,2,'significant');
[~,ind] = find((shiz-0.01)==0.1);

txt1 = ['For a false alarm rate of pfa=0.1,the pmd is '...
    num2str(pmd(ind)) ' and the threshold is t=' num2str(tvec(ind))];
disp(txt1)

%%

% countMD = 0; countFA = 0; countR = 0;
% for i = 1:n
%     if Z(i,1)==1
%         if Z(i,2)==0 %missed detection
%             countMD = countMD+1;
%         else %1&1
%             countR = countR+1;
%         end
%     else
%         if Z(i,2)==1 % false alarm
%             countFA = countFA+1;
%         else %0&0
%             countR = countR+1;
%         end
%     end
% end
% countErr = n-countR;

% error rate
% Perr = countErr/n;
% Perr
% 
% false alarm
% Pfa = countFA/sum(Z(:,1)==0);
% Pfa
% 
% missed detection
% Pmd = countMD/sum(Z(:,1)==1);
% Pmd

