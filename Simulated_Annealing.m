clear all
close all
clc

%% Inlezen van excell file
M = xlsread('Design_space12.xlsx');
%% inputs
xa = M(:,1);
xb = M(:,2); % koud in
xc = M(:,3); % warm in
xd = M(:,4); % NaCl concentration
xe = M(:,5); % Airgap pressure
xf = M(:,6); % Flow
xg = M(:,7); % Q hot
xh = M(:,8); % Q membrane
xi = M(:,9); % Q airgap
xj = M(:,10); % Q permeat
xk = M(:,11); % Q cold
xl = M(:,12); % m membrane
xm = M(:,13); % m airgap
xbb = xb.*xb;
xcc = xc.*xc;
xdd = xd.*xd;
xee = xe.*xe;
xff = xf.*xf;
xgg = xg.*xg;
xhh = xh.*xh;
xii = xi.*xi;
xjj = xj.*xj;
xkk = xk.*xk;
xll = xl.*xl;
xmm = xm.*xm;
xbc = xb.*xc;
xbd = xb.*xd;
xbe = xb.*xe;
xbf = xb.*xf;
xbg = xb.*xg;
xbh = xb.*xh;
xbi = xb.*xi;
xbj = xb.*xj;
xbk = xb.*xk;
xbl = xb.*xl;
xbm = xb.*xm;
xcd = xc.*xd;
xce = xc.*xe;
xcf = xc.*xf;
xcg = xc.*xg;
xch = xc.*xh;
xci = xc.*xi;
xcj = xc.*xj;
xck = xc.*xk;
xcl = xc.*xl;
xcm = xc.*xm;
xde = xd.*xe;
xdf = xd.*xf;
xdg = xd.*xg;
xdh = xd.*xh;
xdi = xd.*xi;
xdj = xd.*xj;
xdk = xd.*xk;
xdl = xd.*xl;
xdm = xd.*xm;
xef = xe.*xf;
xeg = xe.*xg;
xeh = xe.*xh;
xei = xe.*xi;
xej = xe.*xj;
xek = xe.*xk;
xel = xe.*xl;
xem = xe.*xm;
xfg = xf.*xg;
xfh = xf.*xh;
xfi = xf.*xi;
xfj = xf.*xj;
xfk = xf.*xk;
xfl = xf.*xl;
xfm = xf.*xm;
xgh = xg.*xh;
xgi = xg.*xi;
xgj = xg.*xj;
xgk = xg.*xk;
xgl = xg.*xl;
xgm = xg.*xm;
xhi = xh.*xi;
xhj = xh.*xj;
xhk = xh.*xk;
xhl = xh.*xl;
xhm = xh.*xm;
xij = xi.*xj;
xik = xi.*xk;
xil = xi.*xl;
xim = xi.*xm;
xjk = xj.*xk;
xjl = xj.*xl;
xjm = xj.*xm;
xkl = xk.*xl;
xkm = xk.*xm;
xlm = xl.*xm;

Xf = horzcat(xa,xb,xc,xd,xe,xf,xg,xh,xi,xj,xk,xl,xm,xbb,xcc,xdd,xee,xff,xgg,xhh,xii,xjj,xkk,xll,xmm,xbc,xbd,xbe,xbf,xbg,xbh,xbi,xbj,xbk,xbl,xbm,xcd,xce,xcf,xcg,xch,xci,xcj,xck,xcl,xcm,xde,xdf,xdg,xdh,xdi,xdj,xdk,xdl,xdm,xef,xeg,xeh,xei,xej,xek,xel,xem,xfg,xfh,xfi,xfj,xfk,xfl,xfm,xgh,xgi,xgj,xgk,xgl,xgm,xhi,xhj,xhk,xhl,xhm,xij,xik,xil,xim,xjk,xjl,xjm,xkl,xkm,xlm);
[hf,pf] = size(Xf);
nf = length(Xf);

%% ouput
y1 = M(:,16); % flux
y2 = M(:,17); % gor
y3 = M(:,18); % Q heating
y4 = M(:,19); % Q cooling
y5 = M(:,20); % SPEC

%% polynomial regression with interaction terms of flux in coded variables
bf = inv(Xf.'*Xf)*Xf.'*y1; % = inv(Xf.'*Xf)*Xf.'*y1
bg = inv(Xf.'*Xf)*Xf.'*y2;
bh = inv(Xf.'*Xf)*Xf.'*y3;
bc = inv(Xf.'*Xf)*Xf.'*y4;
bs = inv(Xf.'*Xf)*Xf.'*y5;
%% formules voor flux en gor coded values
flux = @(xb,xc,xd,xe,xf,xg,xh,xi,xj,xk,xl,xm) bf(1)+bf(2)*xb+bf(3)*xc+bf(4)*xd+bf(5)*xe+bf(6)*xf+bf(7)*xg+bf(8)*xh+bf(9)*xi+bf(10)*xj+bf(11)*xk+bf(12)*xl+bf(13)*xm+bf(14)*xb.*xb+bf(15)*xc.*xc+bf(16)*xd.*xd+bf(17)*xe.*xe+bf(18)*xf.*xf+bf(19)*xg.*xg+bf(20)*xh.*xh+bf(21)*xi.*xi+bf(22)*xj.*xj+bf(23)*xk.*xk+bf(24)*xl.*xl+bf(25)*xm.*xm+bf(26)*xb.*xc+bf(27)*xb.*xd+bf(28)*xb.*xe+bf(29)*xb.*xf+bf(30)*xb.*xg+bf(31)*xb.*xh+bf(32)*xb.*xi+bf(33)*xb.*xj+bf(34)*xb.*xk+bf(35)*xb.*xl+bf(36)*xb.*xm+bf(37)*xc.*xd+bf(38)*xc.*xe+bf(39)*xc.*xf+bf(40)*xc.*xg+bf(41)*xc.*xh+bf(42)*xc.*xi+bf(43)*xc.*xj+bf(44)*xc.*xk+bf(45)*xc.*xl+bf(46)*xc.*xm+bf(47)*xd.*xe+bf(48)*xd.*xf+bf(49)*xd.*xg+bf(50)*xd.*xh+bf(51)*xd.*xi+bf(52)*xd.*xj+bf(53)*xd.*xk+bf(54)*xd.*xl+bf(55)*xd.*xm+bf(56)*xe.*xf+bf(57)*xe.*xg+bf(58)*xe.*xh+bf(59)*xe.*xi+bf(60)*xe.*xj+bf(61)*xe.*xk+bf(62)*xe.*xl+bf(63)*xe.*xm+bf(64)*xf.*xg+bf(65)*xf.*xh+bf(66)*xf.*xi+bf(67)*xf.*xj+bf(68)*xf.*xk+bf(69)*xf.*xl+bf(70)*xf.*xm+bf(71)*xg.*xh+bf(72)*xg.*xi+bf(73)*xg.*xj+bf(74)*xg.*xk+bf(75)*xg.*xl+bf(76)*xg.*xm+bf(77)*xh.*xi+bf(78)*xh.*xj+bf(79)*xh.*xk+bf(80)*xh.*xl+bf(81)*xh.*xm+bf(82)*xi.*xj+bf(83)*xi.*xk+bf(84)*xi.*xl+bf(85)*xi.*xm+bf(86)*xj.*xk+bf(87)*xj.*xl+bf(88)*xj.*xm+bf(89)*xk.*xl+bf(90)*xk.*xm+bf(91)*xl.*xm;

gor = @(xb,xc,xd,xe,xf,xg,xh,xi,xj,xk,xl,xm) bg(1)+bg(2)*xb+bg(3)*xc+bg(4)*xd+bg(5)*xe+bg(6)*xf+bg(7)*xg+bg(8)*xh+bg(9)*xi+bg(10)*xj+bg(11)*xk+bg(12)*xl+bg(13)*xm+bg(14)*xb.*xb+bg(15)*xc.*xc+bg(16)*xd.*xd+bg(17)*xe.*xe+bg(18)*xf.*xf+bg(19)*xg.*xg+bg(20)*xh.*xh+bg(21)*xi.*xi+bg(22)*xj.*xj+bg(23)*xk.*xk+bg(24)*xl.*xl+bg(25)*xm.*xm+bg(26)*xb.*xc+bg(27)*xb.*xd+bg(28)*xb.*xe+bg(29)*xb.*xf+bg(30)*xb.*xg+bg(31)*xb.*xh+bg(32)*xb.*xi+bg(33)*xb.*xj+bg(34)*xb.*xk+bg(35)*xb.*xl+bg(36)*xb.*xm+bg(37)*xc.*xd+bg(38)*xc.*xe+bg(39)*xc.*xf+bg(40)*xc.*xg+bg(41)*xc.*xh+bg(42)*xc.*xi+bg(43)*xc.*xj+bg(44)*xc.*xk+bg(45)*xc.*xl+bg(46)*xc.*xm+bg(47)*xd.*xe+bg(48)*xd.*xf+bg(49)*xd.*xg+bg(50)*xd.*xh+bg(51)*xd.*xi+bg(52)*xd.*xj+bg(53)*xd.*xk+bg(54)*xd.*xl+bg(55)*xd.*xm+bg(56)*xe.*xf+bg(57)*xe.*xg+bg(58)*xe.*xh+bg(59)*xe.*xi+bg(60)*xe.*xj+bg(61)*xe.*xk+bg(62)*xe.*xl+bg(63)*xe.*xm+bg(64)*xf.*xg+bg(65)*xf.*xh+bg(66)*xf.*xi+bg(67)*xf.*xj+bg(68)*xf.*xk+bg(69)*xf.*xl+bg(70)*xf.*xm+bg(71)*xg.*xh+bg(72)*xg.*xi+bg(73)*xg.*xj+bg(74)*xg.*xk+bg(75)*xg.*xl+bg(76)*xg.*xm+bg(77)*xh.*xi+bg(78)*xh.*xj+bg(79)*xh.*xk+bg(80)*xh.*xl+bg(81)*xh.*xm+bg(82)*xi.*xj+bg(83)*xi.*xk+bg(84)*xi.*xl+bg(85)*xi.*xm+bg(86)*xj.*xk+bg(87)*xj.*xl+bg(88)*xj.*xm+bg(89)*xk.*xl+bg(90)*xk.*xm+bg(91)*xl.*xm;

Qheat = @(xb,xc,xd,xe,xf,xg,xh,xi,xj,xk,xl,xm) bh(1)+bh(2)*xb+bh(3)*xc+bh(4)*xd+bh(5)*xe+bh(6)*xf+bh(7)*xg+bh(8)*xh+bh(9)*xi+bh(10)*xj+bh(11)*xk+bh(12)*xl+bh(13)*xm+bh(14)*xb.*xb+bh(15)*xc.*xc+bh(16)*xd.*xd+bh(17)*xe.*xe+bh(18)*xf.*xf+bh(19)*xg.*xg+bh(20)*xh.*xh+bh(21)*xi.*xi+bh(22)*xj.*xj+bh(23)*xk.*xk+bh(24)*xl.*xl+bh(25)*xm.*xm+bh(26)*xb.*xc+bh(27)*xb.*xd+bh(28)*xb.*xe+bh(29)*xb.*xf+bh(30)*xb.*xg+bh(31)*xb.*xh+bh(32)*xb.*xi+bh(33)*xb.*xj+bh(34)*xb.*xk+bh(35)*xb.*xl+bh(36)*xb.*xm+bh(37)*xc.*xd+bh(38)*xc.*xe+bh(39)*xc.*xf+bh(40)*xc.*xg+bh(41)*xc.*xh+bh(42)*xc.*xi+bh(43)*xc.*xj+bh(44)*xc.*xk+bh(45)*xc.*xl+bh(46)*xc.*xm+bh(47)*xd.*xe+bh(48)*xd.*xf+bh(49)*xd.*xg+bh(50)*xd.*xh+bh(51)*xd.*xi+bh(52)*xd.*xj+bh(53)*xd.*xk+bh(54)*xd.*xl+bh(55)*xd.*xm+bh(56)*xe.*xf+bh(57)*xe.*xg+bh(58)*xe.*xh+bh(59)*xe.*xi+bh(60)*xe.*xj+bh(61)*xe.*xk+bh(62)*xe.*xl+bh(63)*xe.*xm+bh(64)*xf.*xg+bh(65)*xf.*xh+bh(66)*xf.*xi+bh(67)*xf.*xj+bh(68)*xf.*xk+bh(69)*xf.*xl+bh(70)*xf.*xm+bh(71)*xg.*xh+bh(72)*xg.*xi+bh(73)*xg.*xj+bh(74)*xg.*xk+bh(75)*xg.*xl+bh(76)*xg.*xm+bh(77)*xh.*xi+bh(78)*xh.*xj+bh(79)*xh.*xk+bh(80)*xh.*xl+bh(81)*xh.*xm+bh(82)*xi.*xj+bh(83)*xi.*xk+bh(84)*xi.*xl+bh(85)*xi.*xm+bh(86)*xj.*xk+bh(87)*xj.*xl+bh(88)*xj.*xm+bh(89)*xk.*xl+bh(90)*xk.*xm+bh(91)*xl.*xm;

Qcool = @(xb,xc,xd,xe,xf,xg,xh,xi,xj,xk,xl,xm) bc(1)+bc(2)*xb+bc(3)*xc+bc(4)*xd+bc(5)*xe+bc(6)*xf+bc(7)*xg+bc(8)*xh+bc(9)*xi+bc(10)*xj+bc(11)*xk+bc(12)*xl+bc(13)*xm+bc(14)*xb.*xb+bc(15)*xc.*xc+bc(16)*xd.*xd+bc(17)*xe.*xe+bc(18)*xf.*xf+bc(19)*xg.*xg+bc(20)*xh.*xh+bc(21)*xi.*xi+bc(22)*xj.*xj+bc(23)*xk.*xk+bc(24)*xl.*xl+bc(25)*xm.*xm+bc(26)*xb.*xc+bc(27)*xb.*xd+bc(28)*xb.*xe+bc(29)*xb.*xf+bc(30)*xb.*xg+bc(31)*xb.*xh+bc(32)*xb.*xi+bc(33)*xb.*xj+bc(34)*xb.*xk+bc(35)*xb.*xl+bc(36)*xb.*xm+bc(37)*xc.*xd+bc(38)*xc.*xe+bc(39)*xc.*xf+bc(40)*xc.*xg+bc(41)*xc.*xh+bc(42)*xc.*xi+bc(43)*xc.*xj+bc(44)*xc.*xk+bc(45)*xc.*xl+bc(46)*xc.*xm+bc(47)*xd.*xe+bc(48)*xd.*xf+bc(49)*xd.*xg+bc(50)*xd.*xh+bc(51)*xd.*xi+bc(52)*xd.*xj+bc(53)*xd.*xk+bc(54)*xd.*xl+bc(55)*xd.*xm+bc(56)*xe.*xf+bc(57)*xe.*xg+bc(58)*xe.*xh+bc(59)*xe.*xi+bc(60)*xe.*xj+bc(61)*xe.*xk+bc(62)*xe.*xl+bc(63)*xe.*xm+bc(64)*xf.*xg+bc(65)*xf.*xh+bc(66)*xf.*xi+bc(67)*xf.*xj+bc(68)*xf.*xk+bc(69)*xf.*xl+bc(70)*xf.*xm+bc(71)*xg.*xh+bc(72)*xg.*xi+bc(73)*xg.*xj+bc(74)*xg.*xk+bc(75)*xg.*xl+bc(76)*xg.*xm+bc(77)*xh.*xi+bc(78)*xh.*xj+bc(79)*xh.*xk+bc(80)*xh.*xl+bc(81)*xh.*xm+bc(82)*xi.*xj+bc(83)*xi.*xk+bc(84)*xi.*xl+bc(85)*xi.*xm+bc(86)*xj.*xk+bc(87)*xj.*xl+bc(88)*xj.*xm+bc(89)*xk.*xl+bc(90)*xk.*xm+bc(91)*xl.*xm;

%% Model adequacy checking for flux
fluxcheck = zeros(hf,1);

for i = 1:hf
   fluxcheck(i)=flux(xb(i),xc(i),xd(i),xe(i),xf(i),xg(i),xh(i),xi(i),xj(i),xk(i),xl(i),xm(i));
end

residualf = zeros(hf,1);

for i = 1:hf
   residualf(i) = y1(i)-fluxcheck(i);
end

% Normal probability plot of residuals
figure('position', [100, 100, 700, 300])
subplot(1,2,1)
plot1 = normplot(residualf)
set(plot1, 'Color', 'k');
xlabel('Residual [l/m²/h]')
ylabel('Probability [-]')
title('Normal probability plot of the flux')

% Plot of residuals versus predicted response
subplot(1,2,2)
scatter(fluxcheck,residualf,'k')
hold on 
yline(0);
xlabel('Predicted response [l/m²/h]')
ylabel('Residuals [l/m²/h]')
title('Residuals plot of the flux')
% ylim([-2 2])
% xlim([0 20])

x = linspace(0,50);
y = x;

figure('position', [100, 100, 500, 450])
scatter(y1,fluxcheck,'k')
xlabel('J_{sim} [l/m²/h]','LineWidth',1.5,'FontSize',12)
ylabel('J_{reg} [l/m²/h]','LineWidth',1.5,'FontSize',12)
hold on
plot(x,y,'k')
title('Actual versus prediction plot of the flux')
legend('J_{sim} vs J_{reg}','J_{sim} = J_{reg}','Location','northwest','FontSize',12)
axis([0 20 0 20])

%% Model adequacy checking for gor
gorcheck = zeros(hf,1);

for i =1:hf
   gorcheck(i)=gor(xb(i),xc(i),xd(i),xe(i),xf(i),xg(i),xh(i),xi(i),xj(i),xk(i),xl(i),xm(i));
end

residualg = zeros(hf,1);

for i =1:hf
   residualg(i) = y2(i)-gorcheck(i);
end

% Normal probability plot of residuals
figure('position', [100, 100, 700, 300])
subplot(1,2,1)
plot2 = normplot(residualg)
set(plot2, 'Color', 'k');
xlabel('Residual [-]')
ylabel('Probability [-]')
title('Normal probability plot of the GOR')

% Plot of residuals versus predicted response
subplot(1,2,2)
scatter(gorcheck,residualg,'k')
hold on 
yline(0);
xlabel('Predicted response [-]')
ylabel('Residuals [-]')
title('Residuals plot of the GOR')
% ylim([-2 2])
% xlim([0 20])

x = linspace(0,300);
y = x;

figure('position', [100, 100, 500, 450])
scatter(y2,gorcheck,'k')
xlabel('GOR_{sim} [-]')
ylabel('GOR_{reg} [-]')
hold on
plot(x,y,'k')
title('Actual versus prediction plot of the GOR')
legend('GOR_{sim} vs GOR_{reg}','GOR_{sim} = GOR_{reg}','Location','northwest','FontSize',12)
axis([0 16 0 16])

%% Model adequacy checking for Q heating
Qheatcheck = zeros(hf,1);

for i = 1:hf
   Qheatcheck(i)=Qheat(xb(i),xc(i),xd(i),xe(i),xf(i),xg(i),xh(i),xi(i),xj(i),xk(i),xl(i),xm(i));
end

residualh = zeros(hf,1);

for i = 1:hf
   residualh(i) = y3(i)-Qheatcheck(i);
end

% Normal probability plot of residuals
figure('position', [100, 100, 700, 300])
subplot(1,2,1)
plot3 = normplot(residualh)
set(plot3, 'Color', 'k');
xlabel('Residual [kW]')
ylabel('Probability [-]')
title('Normal probability plot of Q_{heating}')

% Plot of residuals versus predicted response
subplot(1,2,2)
scatter(Qheatcheck,residualh,'k')
hold on 
yline(0);
xlabel('Predicted response [kW]')
ylabel('Residuals [kW]')
title('Residuals plot of Q_{heating}')
% ylim([-2 2])
% xlim([0 20])

x = linspace(0,70);
y = x;

figure('position', [100, 100, 500, 450])
scatter(y3,Qheatcheck,'k')
xlabel('Q_{in_{sim}} [kW]','LineWidth',1.5,'FontSize',12)
ylabel('Q_{in_{reg}} [kW]','LineWidth',1.5,'FontSize',12)
hold on
plot(x,y,'k')
title('Actual versus prediction plot of Q_{heating}')
legend('Q_{in_{sim}} vs Q_{in_{reg}}','Q_{in_{sim}} = Q_{in_{reg}}','Location','northwest','FontSize',12)
axis([0 70 0 70]);

%% Model adequacy checking for Q cooling
Qcoolcheck = zeros(hf,1);

for i = 1:hf
   Qcoolcheck(i)=Qcool(xb(i),xc(i),xd(i),xe(i),xf(i),xg(i),xh(i),xi(i),xj(i),xk(i),xl(i),xm(i));
end

residualc = zeros(hf,1);

for i = 1:hf
   residualc(i) = y4(i)-Qcoolcheck(i);
end

% Normal probability plot of residuals
figure('position', [100, 100, 700, 300])
subplot(1,2,1)
plot4 = normplot(residualc)
set(plot4, 'Color', 'k');
xlabel('Residual [kW]')
ylabel('Probability [-]')
title('Normal probability plot of Q_{cooling}')

% Plot of residuals versus predicted response
subplot(1,2,2)
scatter(Qcoolcheck,residualc,'k')
hold on 
yline(0);
xlabel('Predicted response [kW]')
ylabel('Residuals [kW]')
title('Residuals plot of Q_{cooling}')
% ylim([-2 2])
% xlim([0 20])

x = linspace(0,70);
y = x;

figure('position', [100, 100, 500, 450])
scatter(y4,Qcoolcheck,'k')
xlabel('Q_{out_{sim}} [kW]','LineWidth',1.5,'FontSize',12)
ylabel('Q_{out_{reg}}  [kW]','LineWidth',1.5,'FontSize',12) 
hold on
plot(x,y,'k')
title('Actual versus prediction plot of Q_{cooling}')
legend('Q_{out_{sim}} vs Q_{out_{reg}}','Q_{out_{sim}} = Q_{out_{reg}}','Location','northwest','FontSize',12)
axis([0 70 0 70]);

%% Simulated annealing of flux
Tmaxf = 1000; %1000
nepochf = 5000; %5000
Cf = [random('Uniform',15,40),random('Uniform',60,90),random('Uniform',0,60),random('Uniform',10000,1e5),0,random('Uniform',0.5,2),random('Uniform',0.5,2),random('Uniform',0.5,2),random('Uniform',0.5,2),random('Uniform',0.5,2),random('Uniform',0.5,2),random('Uniform',0.5,2)];
solutionf = [];
sig = 0.01; %0.01
q5 = linspace(0.2,2,300);
Efc = zeros(length(q5),1);
Egc = zeros(length(q5),1);
Efn = zeros(length(q5),1);
Egn = zeros(length(q5),1);
h = 1;
aantalruns = [];

for i = 1:nepochf
    a = round(Tmaxf) + 500; %500
    for j =1:a
        for k = 1:length(q5)
            Efc(k) = flux(Cf(1),Cf(2),Cf(3),Cf(4),q5(k),Cf(6),Cf(7),Cf(8),Cf(9),Cf(10),Cf(11),Cf(12));
            Egc(k) = gor(Cf(1),Cf(2),Cf(3),Cf(4),q5(k),Cf(6),Cf(7),Cf(8),Cf(9),Cf(10),Cf(11),Cf(12));
        end
        oppc = trapz(Efc,Egc)-trapz(Egc,Efc);
        N = zeros(12,1);
        %N(1) = 15;
        N(1) = random('Normal',Cf(1),sig);
        if N(1) < 15
            N(1) = 15;
        elseif N(1) > 35
            N(1) = 35;
        end
        %N(2) = 90;
        N(2) = random('Normal',Cf(2),sig);
        if N(2) < 65
            N(2) = 65;
        elseif N(2) > 90
            N(2) = 90;
        end
        %N(3) = 0;
        N(3) = random('Normal',Cf(3),sig);
        if N(3) < 0
            N(3) = 0;
        elseif N(3) > 200
            N(3) = 200;
        end
        %N(4) = 10000;
        N(4) = random('Normal',Cf(4),sig*100000);
        if N(4) < 10000
            N(4) = 10000;
        elseif N(4) > 1e5
            N(4) = 1e5;
        end
        N(6) = random('Normal',Cf(6),sig);
        if N(6) < 0.5
            N(6) = 0.5;
        elseif N(6) > 2
            N(6) = 2;
        end
        N(7) = random('Normal',Cf(7),sig);
        if N(7) < 0.5
            N(7) = 0.5;
        elseif N(7) > 2
            N(7) = 2;
        end
        N(8) = random('Normal',Cf(8),sig);
        if N(8) < 0.5
            N(8) = 0.5;
        elseif N(8) > 2
            N(8) = 2;
        end
        N(9) = random('Normal',Cf(9),sig);
        if N(9) < 0.5
            N(9) = 0.5;
        elseif N(9) > 2
            N(9) = 2;   
        end
        N(10) = random('Normal',Cf(10),sig);
        if N(10) < 0.5
            N(10) = 0.5;
        elseif N(10) > 2
            N(10) = 2;
        end
        N(11) = random('Normal',Cf(11),sig);
        if N(11) < 0.5
            N(11) = 0.5;
        elseif N(11) > 2
            N(11) = 2;
        end
        N(12) = random('Normal',Cf(12),sig);
        if N(12) < 0.5
            N(12) = 0.5;
        elseif N(12) > 2
            N(12) = 2;
        end
        
        for l = 1:length(q5)
            Efn(l) = flux(N(1),N(2),N(3),N(4),q5(l),N(6),N(7),N(8),N(9),N(10),N(11),N(12));
            Egn(l) = gor(N(1),N(2),N(3),N(4),q5(l),N(6),N(7),N(8),N(9),N(10),N(11),N(12));
        end
        oppn = trapz(Efn,Egn)-trapz(Egn,Efn);
        
        dE = oppn-oppc;
        if dE>0
            Cf = N;
        elseif exp(dE/Tmaxf)>rand()
            Cf = N;
        end  
    Edf = flux(Cf(1),Cf(2),Cf(3),Cf(4),Cf(5),Cf(6),Cf(7),Cf(8),Cf(9),Cf(10),Cf(11),Cf(12));
    Edg = gor(Cf(1),Cf(2),Cf(3),Cf(4),Cf(5),Cf(6),Cf(7),Cf(8),Cf(9),Cf(10),Cf(11),Cf(12));
    solutionf(i,1) = i;
    solutionf(i,2) = Tmaxf; 
    solutionf(i,3) = Edf;
    solutionf(i,4) = Edg;
    solutionf(i,5) = oppn;
    solutionf(i,6) = h;
    h = h+1;
    end
        Tmaxf = 0.995*Tmaxf;
end
%%
iterationsf = solutionf(:,6);
fluxsolution = solutionf(:,3);
temperaturef = solutionf(:,2);
oppsol = solutionf(:,5);

figure
subplot(2,1,1)
scatter(iterationsf,oppsol,5,'k','.')
title('Simulated annealing','LineWidth',1.5,'FontSize',12)
xlabel('n_{Iterations}','LineWidth',1.5,'FontSize',12)
ylabel('Surface','LineWidth',1.5,'FontSize',12)

subplot(2,1,2)
plot(iterationsf,temperaturef,'Color',[0.5 0.5 0.5],'LineWidth',1.25)
title('Temperature of annealing','LineWidth',1.5,'FontSize',12)
xlabel('n_{Iterations}','LineWidth',1.5,'FontSize',12)
ylabel(['Annealing' newline 'temperature'],'LineWidth',1.5,'FontSize',12)
