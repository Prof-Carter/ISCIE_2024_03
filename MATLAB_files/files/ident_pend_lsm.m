clear
format compact
close all

myred   = [232   71   70]/255;
myblue	= [  0  112  192]/255;
mygreen	= [ 12  162  116]/255;
mypink	= [255  102  153]/255;
myskyblue    = [  0  176  240]/255;
mydarkpink   = [227   45  145]/255;
mydarkyellow = [228  131   18]/255;

% -----------------------
g  = 9.81e+00;
Tf2 = 0.05;
% -----------------------
load ident_pend_data

tmin = 5;
tmax = 30;

% -----------------------
figure(1)
set(gcf,'Position',[80 80 800 300]) % [x0 y0 width height]
subplot('Position',[0.11 0.2 0.85 0.75])

area([0 tmin], 110*[1 1],'FaceColor',mydarkyellow,'FaceAlpha',0.1,'EdgeAlpha',0)
hold on
area([0 tmin],-110*[1 1],'FaceColor',mydarkyellow,'FaceAlpha',0.1,'EdgeAlpha',0)
plot([tmin tmin],[-110 110],'Color',mydarkyellow)

q1 = stairs(t,phi2*180/pi,'Color',myred,'LineWidth',1.5);
hold off

xlim([0 15])
ylim([-110 110])

set(gca,'XTick',0:1:15)
set(gca,'YTick',-90:30:90)

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$${\phi}_{2}(t)$$ [deg]', 'interpreter', 'latex','FontSize',18)

grid on

% -----
figure(2)
set(gcf,'Position',[120 120 800 300]) % [x0 y0 width height]
subplot('Position',[0.11 0.2 0.85 0.75])

area([tmax t(end)], 110*[1 1],'FaceColor',mydarkyellow,'FaceAlpha',0.1,'EdgeAlpha',0.1)
hold on
area([tmax t(end)],-110*[1 1],'FaceColor',mydarkyellow,'FaceAlpha',0.1,'EdgeAlpha',0.1)
plot([tmax tmax],[-110 110],'Color',mydarkyellow)

r1 = stairs(t,phi2*180/pi,'Color',myred,'LineWidth',1.5);

xlim([15 30])
ylim([-110 110])

set(gca,'XTick',15:1:30)
set(gca,'YTick',-90:30:90)

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$${\phi}_{2}(t)$$ [deg]', 'interpreter', 'latex','FontSize',18)

grid on

% =========================================
k = length(t);
% -----------------------
dphi2(1) = (- 3*phi2(1) + 4*phi2(2) - phi2(3))/(2*h);
for i = 2:k-1
    dphi2(i) = (phi2(i+1) - phi2(i-1))/(2*h);
end
dphi2(k) = (phi2(k-2) - 4*phi2(k-1) + 3*phi2(k))/(2*h);
% -----------------------
ddphi2(1) = (- 3*dphi2(1) + 4*dphi2(2) - dphi2(3))/(2*h);
for i = 2:k-1
    ddphi2(i) = (dphi2(i+1) - dphi2(i-1))/(2*h);
end
ddphi2(k) = (dphi2(k-2) - 4*dphi2(k-1) + 3*dphi2(k))/(2*h);
% -----------------------

% =========================================
M1 = ddphi2;
M2 = dphi2;
N  = - g*sin(phi2);

sysGf2 = tf(1,[Tf2 1])^3;

Mf1 = lsim(sysGf2,M1,t);
Mf2 = lsim(sysGf2,M2,t);
Nf  = lsim(sysGf2,N ,t);

% =========================================
t_data = t;

phi2_data   = phi2;
dphi2_data  = dphi2;
ddphi2_data = ddphi2;

M1_data = M1;
M2_data = M2;
N_data  = N;

Mf1_data = Mf1;
Mf2_data = Mf2;
Nf_data  = Nf;

clear t phi2 dphi2 ddphi2
clear M1 M2 N Mf1 Mf2 Nf
% -----------------------
k = 0;
for i = 1:length(t_data)
    if t_data(i) >= tmin & t_data(i) <= tmax
        k = k + 1;
        t(k,1)   = t_data(i) - tmin;

        phi2(k,1)   = phi2_data(i);
        dphi2(k,1)  = dphi2_data(i);
        ddphi2(k,1) = ddphi2_data(i);

        Mf1(k,1) = Mf1_data(i);
        Mf2(k,1) = Mf2_data(i);
        Nf(k,1)  = Nf_data(i);

        M1(k,1) = M1_data(i);
        M2(k,1) = M2_data(i);
        N(k,1)  = N_data(i);
    end
end

% =======================
figure(3)
set(gcf,'Position',[200 200 800 300]) % [x0 y0 width height]
subplot('Position',[0.11 0.2 0.85 0.75])

p1 = plot(t+tmin,M1,'Color',myred,'LineWidth',1.5);
hold on
p2 = plot(t+tmin,Mf1,'Color',myskyblue,'LineWidth',1.5);
hold off

xlim([0 15])
ylim([-200 200])

set(gca,'XTick',0:1:15)
set(gca,'YTick',-200:50:200)

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$${M}_{1}(t)$$ and $${M}_{\rm f1}(t)$$', 'interpreter', 'latex','FontSize',18)

legend([p1 p2],{'$${M}_{1}(t)=\ddot{\phi}_{2}(t)$$\quad','$${M}_{\rm f1}(t)$$'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',18)

grid on

% -----------------------
figure(4)
set(gcf,'Position',[240 240 800 300]) % [x0 y0 width height]
subplot('Position',[0.11 0.2 0.85 0.75])

p1 = plot(t+tmin,M1,'Color',myred,'LineWidth',1.5);
hold on
p2 = plot(t+tmin,Mf1,'Color',myskyblue,'LineWidth',1.5);
hold off

xlim([15 30])
ylim([-200 200])

set(gca,'XTick',15:1:30)
set(gca,'YTick',-200:50:200)

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$${M}_{1}(t)$$ and $${M}_{\rm f1}(t)$$', 'interpreter', 'latex','FontSize',18)

legend([p1 p2],{'$${M}_{1}(t)=\ddot{\phi}_{2}(t)$$\quad','$${M}_{\rm f1}(t)$$'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',18)

grid on

% =======================
figure(5)
set(gcf,'Position',[320 320 800 300]) % [x0 y0 width height]
subplot('Position',[0.11 0.2 0.85 0.75])

p1 = plot(t+tmin,M2,'Color',myred,'LineWidth',1.5);
hold on
p2 = plot(t+tmin,Mf2,'Color',myskyblue,'LineWidth',1.5);
hold off

xlim([0 15])
ylim([-15 15])

set(gca,'XTick',0:1:15)
set(gca,'YTick',-15:5:15)

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$${M}_{2}(t)$$ and $${M}_{\rm f2}(t)$$', 'interpreter', 'latex','FontSize',18)

legend([p1 p2],{'$${M}_{2}(t)=\dot{\phi}_{2}(t)$$\quad','$${M}_{\rm f2}(t)$$'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',18)

grid on

% -----------------------
figure(6)
set(gcf,'Position',[360 360 800 300]) % [x0 y0 width height]
subplot('Position',[0.11 0.2 0.85 0.75])

p1 = plot(t+tmin,M2,'Color',myred,'LineWidth',1.5);
hold on
p2 = plot(t+tmin,Mf2,'Color',myskyblue,'LineWidth',1.5);
hold off

xlim([15 30])
ylim([-15 15])

set(gca,'XTick',15:1:30)
set(gca,'YTick',-15:5:15)

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$${M}_{2}(t)$$ and $${M}_{\rm f2}(t)$$', 'interpreter', 'latex','FontSize',18)

legend([p1 p2],{'$${M}_{2}(t)=\dot{\phi}_{2}(t)$$\quad','$${M}_{\rm f2}(t)$$'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',18)

grid on

% =======================
figure(7)
set(gcf,'Position',[440 440 800 300]) % [x0 y0 width height]
subplot('Position',[0.11 0.2 0.85 0.75])

p1 = plot(t+tmin,N,'Color',myred,'LineWidth',1.5);
hold on
p2 = plot(t+tmin,Nf,'Color',myskyblue,'LineWidth',1.5);
hold off

xlim([0 15])
ylim([-15 15])

set(gca,'XTick',0:1:15)
set(gca,'YTick',-15:5:15)

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$$N(t)$$ and $${N}_{\rm f}(t)$$', 'interpreter', 'latex','FontSize',18)

legend([p1 p2],{'$$N(t)=-g\sin{\phi}_{2}(t)$$\quad','$${N}_{\rm f}(t)$$'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',18)

grid on

% -----------------------
figure(8)
set(gcf,'Position',[480 480 800 300]) % [x0 y0 width height]
subplot('Position',[0.11 0.2 0.85 0.75])

p1 = plot(t+tmin,N,'Color',myred,'LineWidth',1.5);
hold on
p2 = plot(t+tmin,Nf,'Color',myskyblue,'LineWidth',1.5);
hold off

xlim([15 30])
ylim([-15 15])

set(gca,'XTick',15:1:30)
set(gca,'YTick',-15:5:15)

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$$N(t)$$ and $${N}_{\rm f}(t)$$', 'interpreter', 'latex','FontSize',18)

legend([p1 p2],{'$$N(t)=-g\sin{\phi}_{2}(t)$$\quad','$${N}_{\rm f}(t)$$'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',18)

grid on

% =========================================
Mf = [ Mf1  Mf2 ];

p2 = inv(Mf'*Mf)*Mf'*Nf;

a2 = p2(1);
b2 = p2(2);

fprintf('a2 = %3.2e;\n',a2)
fprintf('b2 = %3.2e;\n',b2)

% -----------------------
[phi20 i] = max(phi2);
t0 = t(i);

sim('sim_nonlinear_pend')

% ================================
figure(1)
hold on
q2 = plot(tsim+tmin,phi2sim*180/pi,'--','Color',mygreen,'LineWidth',2);
plot(tsim(1)+tmin,phi2sim(1)*180/pi,'o','Color',mydarkyellow,'LineWidth',1.5,'MarkerSize',8)
hold off

legend([q1 q2],{'Experiment','Simulation'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',18)

% ================================
figure(2)
hold on
r2 = plot(tsim+tmin,phi2sim*180/pi,'--','Color',mygreen,'LineWidth',2);
hold off

legend([r1 r2],{'Experiment','Simulation'},'Location','southeast','NumColumns',2)
set(legend, 'interpreter', 'latex','FontSize',18)

% ================================
figure(1); movegui('north')
figure(2); movegui('northeast')
figure(3); movegui('west')
figure(4); movegui('center')
figure(5); movegui('east')
figure(6); movegui('southwest')
figure(7); movegui('south')
figure(8); movegui('southeast')

% ================================
% figure(1); exportgraphics(gcf,'ident_pend1_sim.jpg','Resolution',1200)
% figure(2); exportgraphics(gcf,'ident_pend2_sim.jpg','Resolution',1200)
% figure(3); exportgraphics(gcf,'ident_pend1_M1.jpg','Resolution',1200)
% figure(4); exportgraphics(gcf,'ident_pend2_M1.jpg','Resolution',1200)
% figure(5); exportgraphics(gcf,'ident_pend1_M2.jpg','Resolution',1200)
% figure(6); exportgraphics(gcf,'ident_pend2_M2.jpg','Resolution',1200)
% figure(7); exportgraphics(gcf,'ident_pend1_N.jpg','Resolution',1200)
% figure(8); exportgraphics(gcf,'ident_pend2_N.jpg','Resolution',1200)

