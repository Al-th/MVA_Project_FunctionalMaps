clc;
clear;

f = 1;
t = linspace(0,1,10000);

x1 = cos(2*pi*f*t);
x2 = sin(2*pi*f*t);

x1 = x1';
x2 = x2';



PHI1 = [x1,x2];
PHI2 = [x2,x1];




a = PHI1'*PHI1;
b = PHI2'*PHI2;

C = b/a;

%Newfunc in original space 
x3 = (2*cos(2*pi*f*t))';

x3ImageCoeffs = C*(PHI1'*x3);
x3ImageValues= zeros(10000,1);
for i = 1:2
    x3ImageValues = x3ImageValues + x3ImageCoeffs(i)/5000*PHI2(:,i);
end

x4 = cos(2*pi*10*f*t)';

x4ImageCoeffs = C*(PHI1'*x4);
x4ImageValues= zeros(10000,1);
for i = 1:2
    x4ImageValues = x4ImageValues + x4ImageCoeffs(i)/5000*PHI2(:,i);
end

alpha = 0.5;
x5 = alpha*cos(2*pi*f*t)'+(1-alpha)*sin(2*pi*f*t)';
x5ImageCoeffs = C*(PHI1'*x5);
x5ImageValues= zeros(10000,1);
for i = 1:2
    x5ImageValues = x5ImageValues + x5ImageCoeffs(i)/5000*PHI2(:,i);
end




clf;
subplot(2,4,1);
plot(x1);
ylim([-2,2]);
subplot(2,4,2);
ylim([-2,2]);
plot(x3,'r')
subplot(2,4,3);
ylim([-2,2]);
plot(x4,'c');
subplot(2,4,4);
ylim([-2,2]);
plot(x5,'b');
ylim([-2,2]);


subplot(2,4,5);
plot(x2);
ylim([-2,2]);
subplot(2,4,6);
plot(x3ImageValues,'r');
ylim([-2,2]);
subplot(2,4,7);
plot(x4ImageValues,'c');
ylim([-2,2]);
subplot(2,4,8);
plot(x5ImageValues,'b');
ylim([-2,2]);