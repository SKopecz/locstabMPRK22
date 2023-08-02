clear, clc, close all

%% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear test problem
a = 20;
b = a;
f = @(t,y) [b*y(2) - a*y(1); a*y(1) - b*y(2)];

tspan = [0 100];
dt0 = 1;
alpha = -0.5;

delta1 = 0.23;
delta2 = 0.24;

figure(1)
hold on
y0 = 0.5 + delta1*[-1;1];
[t, y] = MPRK22(y0, dt0, tspan,alpha,a,b);
[tref,yref] = ode23s(f,tspan,y0);
plot(tref,yref,'--')
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,y,'o')
ylim([0 1])
lgd = legend('$y_1$ exact','$y_2$ exact','$y_1$','$y_2$');
set(lgd,'Interpreter','Latex','Fontsize',14)
set(gca,'FontSize',14)
hold off

% print -depsc2 ../stable.eps

figure(2)
hold on
y0 = 0.5 + delta2*[-1;1];
[t, y] = MPRK22(y0, dt0, tspan,alpha,a,b);
[tref,yref] = ode23s(f,tspan,y0);
plot(tref,yref,'--')
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,y,'o')
ylim([0 1])
lgd = legend('$y_1$ exact','$y_2$ exact','$y_1$','$y_2$');
set(lgd,'Interpreter','Latex','Fontsize',14,'Location','east')
set(gca,'FontSize',14)
hold off

% print -depsc2 ../unstable.eps
  
%% Figure 2 a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

% Linear test problem
a = 20;
b = a;

tspan = [0 1e4];
dt0 = 1;
alpha = -0.5;

ss_diff = [];
delta = [];
N = 200;
for i = 1:N-1
  delta(i) = 0.5*(i-1)/N;
  y0 = 0.5 + delta(i)*[1; -1];

  [t, y] = MPRK22(y0, dt0, tspan,alpha,a,b);
  ss_diff(i) = norm(y(:,end) - [0.5;0.5],'inf');
end

figure(3)
semilogy(delta,ss_diff,'.','MarkerSize',10)
set(gca,'FontSize',14)
xlabel('$\delta$','Interpreter','Latex','FontSize',20)
ylabel('$d(-0.5,\delta)$','Interpreter','Latex','FontSize',20)
xlim([0 0.5])
ylim([1e-16 1e2])
grid on
grid minor
shg

% print -depsc2 ../d_delta.eps

%% Figure 2 b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear test problem
a = 200;
b = a;

tspan = [0 1e4];
dt0 = 1;
alpha = -0.5;

ss_diff = [];
delta = [];
N = 200;
for i = 1:N-1
  delta(i) = 0.5*(i-1)/N;
  y0 = 0.5 + delta(i)*[1; -1];

  [t, y] = MPRK22(y0, dt0, tspan,alpha,a,b);
  ss_diff(i) = norm(y(:,end) - [0.5;0.5],'inf');
end

figure(4)
semilogy(delta,ss_diff,'.','MarkerSize',10)
set(gca,'FontSize',14)
xlabel('$\delta$','Interpreter','Latex','FontSize',20)
ylabel('$d(-0.5,\delta)$','Interpreter','Latex','FontSize',20)
xlim([0 0.5])
ylim([1e-16 1e2])
grid on
grid minor
shg

% print -depsc2 ../d_delta2.eps

%% Figure 3 a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear test problem
a = 200;
b = a;

tspan = [0 1e4];
dt0 = 1;

% Increase Na and Nd for a higher resolution
Na = 25; %Na = 241;
Nd = 16; %Nd = 160;

delta = linspace(0,0.5,Nd+2);
delta = delta(2:end-1);

alpha = zeros(1,Na);
for j = 1:Na
 alpha(j) = -1.5 + 3*(j-1)/(Na-1);
end

ss_diff2 = zeros(Na,Nd);
parfor j = 1:Na
  for i = 1:Nd
    y0 = 0.5 + 0.5*i/(Nd+1)*[1; -1]; %y0 = 0.5 + delta(i)*[1; -1];
    try
      [t, y] = MPRK22(y0, dt0, tspan,alpha(j),a,b);
      ss_diff2(j,i) = norm(y(:,end) - [0.5;0.5],'inf');
    catch
      ss_diff2(j,i) = Inf;
    end
  end
end

% save ss_diff2.mat ss_diff2

figure(5)
imagesc(delta,alpha,ss_diff2,[1e-16 1e-1])
%%% Uncomment to see grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on
% dx = delta(2)-delta(1);
% dy = alpha(2)-alpha(1);
% xg = linspace(delta(1) - dx/2,delta(end)+dx/2,Nd+1);
% yg = linspace(alpha(1) - dy/2,alpha(end)+dy/2,Na+1);
% hm = mesh(xg,yg,zeros(Na+1,Nd+1));
% hm.FaceColor = 'none';
% hm.EdgeColor = 'k';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontSize',14)
set(gca,'ColorScale','log')
set(gca,'YDir','normal')
colorbar
xlabel('$\delta$','Interpreter','Latex','FontSize',20)
ylabel('$\alpha$','Interpreter','Latex','FontSize',20)
hold off

% print -depsc2 ../d_alpha_delta.eps
%% Figure 3 b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear test problem
a = 200;
b = a;

tspan = [0 1e4];
dt0 = 1;

% Increase Na and Nd for a higher resolution
Na = 25; %Na = 241;
Nd = 16; %Nd = 160;

delta = linspace(0,0.5,Nd+2);
delta = delta(2:end-1);

alpha = zeros(1,Na);
for j = 1:Na
 alpha(j) = -0.6 + 0.1*(j-1)/(Na-1);
end

ss_diff2b = zeros(Na,Nd);
parfor j = 1:Na
  for i = 1:Nd
    y0 = 0.5 + 0.5*i/(Nd+1)*[1; -1]; %y0 = 0.5 + delta(i)*[1; -1];
    try
      [t, y] = MPRK22(y0, dt0, tspan,alpha(j),a,b);
      ss_diff2b(j,i) = norm(y(:,end) - [0.5;0.5],'inf');
    catch
      ss_diff2b(j,i) = Inf;
    end
  end
end

% save ss_diff2b.mat ss_diff2b

figure(6)
imagesc(delta,alpha,ss_diff2b,[1e-16 1e-1])
%%% Uncomment to see grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on
% dx = delta(2)-delta(1);
% dy = alpha(2)-alpha(1);
% xg = linspace(delta(1) - dx/2,delta(end)+dx/2,Nd+1);
% yg = linspace(alpha(1) - dy/2,alpha(end)+dy/2,Na+1);
% hm = mesh(xg,yg,zeros(Na+1,Nd+1));
% hm.FaceColor = 'none';
% hm.EdgeColor = 'k';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontSize',14)
set(gca,'ColorScale','log')
set(gca,'YDir','normal')
colorbar
xlabel('$\delta$','Interpreter','Latex','FontSize',20)
ylabel('$\alpha$','Interpreter','Latex','FontSize',20)
hold off

% print -depsc2 ../d_alpha_delta2.eps
