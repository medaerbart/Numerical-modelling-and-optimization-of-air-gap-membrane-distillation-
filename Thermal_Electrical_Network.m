clear all

N = 100; % Number of cells/sections
%% Mass matrix
E = zeros(8*N+4,8*N+4);
for i = 1:4*N+2
    E(i,i)=1;
end

%% Reading data from excel file
M = xlsread('Cal_data_el.xlsx');

Tfin = M(:,1); Tpin = M(:,2); Flowin = M(:,3); Sin=M(:,4); P = M(:,5); Lmd = M(:,6); Nchan = M(:,7);
J = zeros(length(Tfin),1);
Tpout = zeros(length(Tfin),1);
Tfout = zeros(length(Tfin),1);
Tfoutbeter = zeros(length(Tfin),1);
Tpoutbeter = zeros(length(Tfin),1);

for q = [ 6 14 21 ]

    %% initial conditions
    initial = 4500*ones(8*N+4,1);

    for i = 1:2*N
       initial(i) = 10000; 
    end
    
    f = linspace(Tfin(q),Tpin(q)+15,N);
    h = 1;
    for i = 2:2:2*N
        initial(i) = f(h);
        h = h+1;
    end
    h = 2;
    for i = 4*N+5:5*N+4
       initial(i) = initial(h) - 3 ; 
       h = h +2;
    end

    for i = 5*N+5:6*N+4
       initial(i) = initial(i-N) - 3; 
    end

    for i = 6*N+5:7*N+4
        initial(i) = initial(i-N) - 3;
    end

    for i = 7*N+5:8*N+4
       initial(i) = initial(i-N) - 3; 
    end
    f = linspace(Tfin(q)-15,Tpin(q),N);
    h = 1;
    for i = 2*N+3:2:4*N+1
         initial(i) = f(h);
         h = h+1;
    end
    
    initial(4*N+4) = Tpin(q)+Tfin(q)/2;
    initial(4*N+3) = Tfin(q)-Tpin(q)/1.25;
    %% ode15s
    tmax = 30;
    tspan = [0 tmax];
    options = odeset('Mass',E,'MassSingular','yes','RelTol',1e-1);
    [t,y] = ode15s(@(t,y) versie2_cal(t,y,N,Tpin(q),Tfin(q),Flowin(q),P(q),Sin(q),Nchan(q),Lmd(q)),tspan,initial,options);
    
    [m,n] = size(y);
    
    %% Getting the temperatures out of the model
    Tpout(q) = y(m,4*N+4);
    Tfout(q) = y(m,4*N+3);
    
    Tfoutbeter(q) = y(m,2*N);
    Tpoutbeter(q) = y(m,2*N+3);
    %% Getting the flux out of the model
    f = zeros(8*N+4,1);
    for i = 1:8*N+4
        f(i)  = y(m,i);
    end
    
    J(q) = sum(Data_cal(t,f,N,Tpin(q),Tfin(q),Flowin(q),P(q),Sin(q),Nchan(q),Lmd(q)))*3600;
end
