close all
clear all
clc

%% Parameters
global lambda gam kap 
lambda=.9;
gam=.4;
kap=.01;
S0=1;
global delta
global a
global b
delta=1e-6; a=4; b=4;
Sbar=exp(-(1-lambda)/(kap*lambda));





%% Loop
%%% Loop Index
index=0;
segment=1;
%%% velocity guesses
vd=.01; %% the lower guess
vu=.02; %% the upper guess
v=vd;%initial guess


tstart =cputime;%% Measures CPU time

while (v>1e-4)||(v<1)||(index<20)
index=index+1;


%% The computation part
    
 IC1=[1e-3,S0];   %% the initial condition is away from (0,1) since it is an equilibrium point
 
Model = @(y)[v*(y(1)-ell(y(2))/gam); -gam*y(1)*y(2)*Dd(y(1))/(v*(kap+y(2)));];%+ kap*log(max(y(2)/S0,eps))/gam
%[t,y] = ode45(Model,[0,1e3],IC);

% Runge Kutta 4th order
clear Mpos Spos wpos t1 y;
eps=1e-10;
y(:,1)=IC1;
t(1)=0;
n=1e9;
h=1e-5;

%%% 4-th order Runge-Kutta
for i = 1:n
    
   h=dtPos(y(:,i));
    
    k1 = Model(y(:,i));
    k2 = Model(y(:,i)+.5*k1*h);
    k3 = Model(y(:,i)+.5*k2*h);
    k4 = Model(y(:,i)+k3*h);
    y(:,i+1) = y(:,i)+((k1+2*k2+2*k3+k4)/6)*h;
    
    
    t(i+1)=t(i)+h;
    

 
    
 %%% Stopping Criteria   
    if (y(1,i+1)<=0)
        segment=1;
       break; 
    elseif y(2,i+1)<=Sbar
        segment=-1;
        break;
    elseif (y(2,i+1)<eps)&&(y(1,i+1)>ell(y(2,i+1))/gam)
        segment=-2;
        break;
    end
    
end

if (index==1)&&(segment<0)%%% If the lower bound is incorrect
       fprintf('ERROR:  incorrect lower v prescribed\n');
       plot(y(1,:),y(2,:));
       break;
elseif (index==2)&&(segment>0)%%% If the upper bound is incorrect
       fprintf('ERROR:  incorrect upper v prescribed\n');
       plot(y(1,:),y(2,:));
       break;
elseif y(2,end)>.95%%% If the initial condition is too close to the equilibrium point
    fprintf('ERROR:  initial condition too close to (0,1)\n');
    break;
end   



%% Update v



%%% The Bisection Algorithm
if index==1
v=vu;
elseif segment<0
vu=v;
yold=y;
else
vd=v;
%yold=y;
end
v=.5*(vu+vd);
dv=(vu-vd)/v;

fprintf('Iteration= %d, wave-speed v=%f, v relative error= %f, in time %f\n',index,v,dv,cputime-tstart);

    
if(abs(vu-vd)/v)<1e-6%%% stops if tolerance is less than 0.1%
    y=yold;
    fprintf('THE FINAL WAVESPEED v=%f\n',v);

%% Finishing    



%%% Assinging variables        
M=y(1,1:(end-1));
S=y(2,1:(end-1));
w=v*(S0/lambda + (S-S0)+kap*log(S/S0))/gam;


 %%The Zeta computations
 clear zeta1;
zeta(1)=0;
i=length(M);
for j=2:i
%zeta(j)=zeta(j-1)-.5*(Dd(M(j))+Dd(M(j-1)))*(t(j)-t(j-1));
zeta(j)=zeta(j-1)+v*(1+kap/S(j))*(S(j)-S(j-1))/(M(j)*gam);
end

%%% Breaking the Loop   
    break;
end    
    


end

%% Plotting
%%% against TW coordinate
figure(1);
plot(zeta,M,'-g','LineWidth',2);
hold on;
plot(zeta,S,'-.b','LineWidth',2);
plot(zeta,w,':k','LineWidth',3);
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.FontWeight = 'bold';

lgd =legend("$M$","$S$","$\omega$",'interpreter','latex','Location','northwest');
lgd.FontSize = 20;
     xh= xlabel('$\xi$','interpreter','latex','fontsize',30);

 name=strcat('$\delta=', num2str(delta), '$, $a=b=$', num2str(a),', $\lambda=$',num2str(lambda)); 
    title(name,'interpreter','latex','fontsize',24);
%%%PLotting the tail part
zeta2=min(zeta):-.001:-1;
Ml=M(end)*exp(lambda*(zeta2-max(zeta2))/v);
wl=v*Ml/lambda;
plot(zeta2,Ml,'-g','LineWidth',2,'HandleVisibility','off');
plot(zeta2,wl,':k','LineWidth',3,'HandleVisibility','off');


%%% In the phase plane

figure(2);
plot(M,S,'b','LineWidth',2);
hold on;
plot([0 1],[S0 S0],'--k','HandleVisibility','off');
ss=0:.005:1;
plot(ell(ss)/gam,ss,'-.r','LineWidth',.5,'HandleVisibility','off');

%%tail part
Sl=exp(-(1-lambda-gam*Ml)/(kap*lambda));
plot(Ml,Sl,'b','LineWidth',2,'HandleVisibility','off');
grid on;
box on;
% view(135,20)

 axis([0 1 0 1]);
   
   
     xh= xlabel('$M$','interpreter','latex','fontsize',24);
    yh= ylabel('$S$','interpreter','latex','fontsize',24);
%     zh= zlabel('$S$','interpreter','latex','fontsize',24);
    
    % Get handle to current axes.
ax = gca;
% The below would set everything: title, x axis, y axis, and tick mark label font sizes.
 ax.FontSize = 16;
ax.FontWeight = 'bold';




%% Functions
%%Diffusion function
function d=Dd(m)
global a
global b
global delta

m=max(min(m,1),0);
d=(m.^a)./((1-m).^b);
d=delta*d;
end

function m=ell(s)
global lambda
global kap

m=(1-lambda)*(1-s)+ kap*lambda*log(s);
end

function dt=dtPos(y)
%%% Time stepping of the positive half
fac=10;
if y(1)<.2
dt=fac*1e-2;
else
dt=fac*1e-3;    
end

end
