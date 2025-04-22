syms y(t)
ode = diff(y,t) == y*t*t -1.1*y;
cond = y(0)==1;
ySol(t)=dsolve (ode,cond)

fplot(ySol,[0,2]);
grid on;
grid minor;
xlabel('t','FontWeight','bold')
ylabel('y','FontWeight','bold')
Title("Analytical Solution Plot")