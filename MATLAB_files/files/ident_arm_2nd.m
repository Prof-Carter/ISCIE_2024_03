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
load ident_arm_data

% -----------------------
rc = ref*pi/180;           % アーム角度の目標値
% -----------------------
t_data = t;                 % いったん，データを移動
theta1_data = theta1;
v_data = v;
clear t theta1 v
% -----------------------
k = 0;
for i = 1:length(t_data)    % 1 秒までのデータを捨てる
   if t_data(i) >= 1;
       k = k + 1;
       t(k) = t_data(i) - 1;
       theta1(k) = theta1_data(i);
       v(k) = v_data(i);
   end
end

% ++++++++++ Step 1 ++++++++++
[theta1_max imax] = max(theta1);    % ステップ応答の最大値の抽出

imax_last = imax;
for i = imax:length(t)
    if theta1(i) == theta1_max
        imax_last = i;
    end
    if theta1(i) < theta1_max
        break
    end
end

Tp    = (t(imax) + t(imax_last))/2; % 行き過ぎ時間
Amax  = theta1_max - rc;            % オーバーシュート

fprintf('Tp   = %3.2e;\n',Tp)
fprintf('Amax = %3.2e;\n',Amax)

% ++++++++++ Step 2 ++++++++++
xi   = - (1/Tp)*log(Amax/rc);
wn   = sqrt((pi/Tp)^2 + xi^2);      % 固有角周波数
zeta = xi/wn;                       % 減衰係数

fprintf('wn   = %3.2e;\n',wn)
fprintf('zeta = %3.2e;\n',zeta)

% ++++++++++ Step 3 ++++++++++
a1 = 2*zeta*wn;                     % a1 の同定
b1 = wn^2/kP;                       % b1 の同定

fprintf('a1   = %3.2e;\n',a1)
fprintf('b1   = %3.2e;\n',b1)

% -----------------------
t_sim = 0:0.001:1;                  % 同定されたパラメータを用いてシミュレーション
theta1_sim = rc*step(wn^2,[1 2*zeta*wn wn^2],t_sim);

figure(1)
movegui('north')
subplot('Position',[0.15 0.15 0.775 0.775])

plot([Tp Tp],[0 theta1_max]*180/pi,'Color',myred)
hold on
plot([0 Tp],[theta1_max theta1_max]*180/pi,'Color',myred)
plot([0 1],[rc rc]*180/pi,'Color',myred)

p1 = plot(t_sim,theta1_sim*180/pi,'--','Color',mygreen,'LineWidth',2);
p2 = stairs(t,theta1*180/pi,'Color',myblue,'LineWidth',1.5);

plot(Tp,theta1_max*180/pi,'o','LineWidth',1.5,'MarkerSize',8,'Color',myred)
hold off
xlim([0 1]);
ylim([0 120])

set(gca,'FontName','arial','FontSize',16)
xlabel('$$t$$ [s]', 'interpreter', 'latex','FontSize',18)
ylabel('$${\theta}_{1}(t)$$ [deg]', 'interpreter', 'latex','FontSize',18)

xtickangle(0)
set(gca,'XTick',0:0.1:1)
set(gca,'YTick',0:30:120)

legend([p2 p1],{'Experiment','Simulation'},'Location','southeast')
set(legend, 'interpreter','latex','FontSize',18)

grid on



