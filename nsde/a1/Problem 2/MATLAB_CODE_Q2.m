
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SOLVING GIVEN EQUATION ANALYTICALLY AND PLOTTING IT
t = linspace(0, 2, 100);
ySol = exp((t.*(10*t.^2 - 33))/30);
% ySol(t) = exp((t*(10*t^2 - 33))/30);

figure
% fplot(ySol,[0,2]);
plot(t, ySol)
grid on;
grid minor;
xlabel('t','FontWeight','bold')
ylabel('y','FontWeight','bold')
title("Analytical Solution Plot")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOTTING EULER METHOD, ADAM BASHFORTH METHOD AND FORTH ORDER RUNGE KUTTA
%%METHOD FOR h =0.03125

load Euler2b1_0.031250.csv

mEE= zeros(5,1); %%Generating array for Maximum Euler Error
aEE= zeros(5,1); %%Generating array for Average Euler Error
plot (Euler2b1_0_031250(:,1),Euler2b1_0_031250(:,2));
Nt=size(Euler2b1_0_031250,1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mEE(1)= max(ynew-Euler2b1_0_031250(:,2));
aEE(1) =mean(ynew-Euler2b1_0_031250(:,2));


hold on
load AB2b2_0.031250.csv

mAB= zeros(5,1); %%Generating array for Maximum Adam Bashforth
aAB= zeros(5,1); %%Generating array for Average Adam Bashforth
plot (AB2b2_0_031250(:,1),AB2b2_0_031250(:,2));
Nt=size(AB2b2_0_031250(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mAB(1)= max(ynew-AB2b2_0_031250(:,2));
aAB(1)= mean(ynew-AB2b2_0_031250(:,2));


hold on
load RK2b3_0.031250.csv

mRK= zeros(5,1); %%Generating array for Maximum Runge Kutta
aRK= zeros(5,1); %%Generating array for Average Runge Kutta
plot (RK2b3_0_031250(:,1),RK2b3_0_031250(:,2));
Nt= size(RK2b3_0_031250(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mRK(1)= max(ynew-RK2b3_0_031250(:,2));
aRK(1)= mean(ynew-RK2b3_0_031250(:,2));


legend('Euler','Adam Bashforth','Runge Kutta')
xlabel('t','FontWeight','bold')
ylabel('y','FontWeight','bold')
title("Function plot for step size= 0.03125")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOTTING EULER METHOD, ADAM BASHFORTH METHOD AND FORTH ORDER RUNGE KUTTA
%%METHOD FOR h =0.0625

figure
load Euler2b1_0.062500.csv
plot (Euler2b1_0_062500(:,1),Euler2b1_0_062500(:,2));
Nt= size(Euler2b1_0_062500(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mEE(2)= max(ynew-Euler2b1_0_062500(:,2));
aEE(2) =mean(ynew-Euler2b1_0_062500(:,2));



hold on
load AB2b2_0.062500.csv
plot (AB2b2_0_062500(:,1),AB2b2_0_062500(:,2));
Nt= size(AB2b2_0_062500(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mAB(2)= max(ynew-AB2b2_0_062500(:,2));
aAB(2)= mean(ynew-AB2b2_0_062500(:,2));




hold on
load RK2b3_0.062500.csv
plot (RK2b3_0_062500(:,1),RK2b3_0_062500(:,2));
Nt= size(RK2b3_0_062500(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mRK(2)= max(ynew-RK2b3_0_062500(:,2));
aRK(2)= mean(ynew-RK2b3_0_062500(:,2));




legend('Euler','Adam Bashforth','Runge Kutta')
xlabel('t','FontWeight','bold')
ylabel('y','FontWeight','bold')
title("Function plot for step size= 0.0625") 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOTTING EULER METHOD, ADAM BASHFORTH METHOD AND FORTH ORDER RUNGE KUTTA
%%METHOD FOR h =0.125


figure
load Euler2b1_0.125000.csv
plot (Euler2b1_0_125000(:,1),Euler2b1_0_125000(:,2));
Nt= size(Euler2b1_0_125000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mEE(3)= max(ynew-Euler2b1_0_125000(:,2));
aEE(3) =mean(ynew-Euler2b1_0_125000(:,2));



hold on
load AB2b2_0.125000.csv
plot (AB2b2_0_125000(:,1),AB2b2_0_125000(:,2));
Nt= size(AB2b2_0_125000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mAB(3)= max(ynew-AB2b2_0_125000(:,2));
aAB(3)= mean(ynew-AB2b2_0_125000(:,2));


hold on
load RK2b3_0.125000.csv
plot (RK2b3_0_125000(:,1),RK2b3_0_125000(:,2));
Nt =size(RK2b3_0_125000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mRK(3)= max(ynew-RK2b3_0_125000(:,2));
aRK(3)= mean(ynew-RK2b3_0_125000(:,2));



legend('Euler','Adam Bashforth','Runge Kutta')
xlabel('t','FontWeight','bold')
ylabel('y','FontWeight','bold')
title("Function plot for step size= 0.125") 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOTTING EULER METHOD, ADAM BASHFORTH METHOD AND FORTH ORDER RUNGE KUTTA
%%METHOD FOR h =0.25


figure
load Euler2b1_0.250000.csv
plot (Euler2b1_0_250000(:,1),Euler2b1_0_250000(:,2));
Nt= size(Euler2b1_0_250000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mEE(4)= max(ynew-Euler2b1_0_250000(:,2));
aEE(4) =mean(ynew-Euler2b1_0_250000(:,2));


hold on
load AB2b2_0.250000.csv
plot (AB2b2_0_250000(:,1),AB2b2_0_250000(:,2));
Nt= size(AB2b2_0_250000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mAB(4)= max(ynew-AB2b2_0_250000(:,2));
aAB(4)= mean(ynew-AB2b2_0_250000(:,2));



hold on
load RK2b3_0.250000.csv
plot (RK2b3_0_250000(:,1),RK2b3_0_250000(:,2));
Nt= size(RK2b3_0_250000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mRK(4)= max(ynew-RK2b3_0_250000(:,2));
aRK(4)= mean(ynew-RK2b3_0_250000(:,2));


legend('Euler','Adam Bashforth','Runge Kutta')
xlabel('t','FontWeight','bold')
ylabel('y','FontWeight','bold')
title("Function plot for step size= 0.25") 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOTTING EULER METHOD, ADAM BASHFORTH METHOD AND FORTH ORDER RUNGE KUTTA
%%METHOD FOR h =0.5


figure
load Euler2b1_0.500000.csv
plot (Euler2b1_0_500000(:,1),Euler2b1_0_500000(:,2));
Nt= size(Euler2b1_0_500000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mEE(5)= max(ynew-Euler2b1_0_500000(:,2));
aEE(5) =mean(ynew-Euler2b1_0_500000(:,2));




hold on
load AB2b2_0.500000.csv
plot (AB2b2_0_500000(:,1),AB2b2_0_500000(:,2));
Nt= size(AB2b2_0_500000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mAB(5)= max(ynew-AB2b2_0_500000(:,2));
aAB(5)= mean(ynew-AB2b2_0_500000(:,2));



hold on
load RK2b3_0.500000.csv
plot (RK2b3_0_500000(:,1),RK2b3_0_500000(:,2));
Nt= size(RK2b3_0_500000(:,1),1);
tnew=linspace(0,2,Nt);
ynew = exp((tnew.*(10*tnew.^2 - 33))/30)';
mRK(5)= max(ynew-RK2b3_0_500000(:,2));
aRK(5)= mean(ynew-RK2b3_0_500000(:,2));



legend('Euler','Adam Bashforth','Runge Kutta')
xlabel('t','FontWeight','bold')
ylabel('y','FontWeight','bold')
title("Function plot for step size= 0.5") 


figure
h= [0.031250 0.062500 0.125000 0.250000 0.500000]; 
loglog(h,mEE,h,mAB,h,mRK)
% xlim([])
% ylim([])
xlabel("Logarithmic of stepsizes")
ylabel("Logarithmic of Maximum Error")
grid on
grid minor
xsl = [0.5 0.125];
ysl = [0.001536 6e-6];
hold on
loglog(xsl, ysl, Color='k')
legend('Max Euler Method Error', 'Max Adam Bashforth Method Error', 'Max Runge Kutta Method Error','Slope = 4 line', 'Location', 'southeast')
title("Max Error Comparison")


figure
loglog(h,aEE,h,aAB,h,aRK)
% xlim([])
% ylim([])
legend('Mean Euler Method Error', 'Mean Adam Bashforth Method Error', 'Mean Runge Kutta Method Error', 'Location', 'southeast')
xlabel("Logarithmic of stepsizes")
ylabel("Logarithmic of Mean Error")
grid on
grid minor
title("Mean Error Comparison")


