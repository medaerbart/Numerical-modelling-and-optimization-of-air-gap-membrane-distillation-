function ydot = pde_7cal(t,y,N,Tpin,Tfin,Flowin,P,Sin,Nchan,Lmd)
%%% This fill will not work when you run it because the parameters are not
%%% included in this file. This file is only to show how the code is
%%% structured.

%%% Variable input conditions
% Always after t>5 so the initial conditions are right.
%%% Process conditions
%Tpin = 21.1289047798268;                       % Cold temperature in          [°C]
%Tfin = 75.7121942039893;                       % Hot temperature in           [°C]
%Flowin = 0.4594;                               % Volume flow [m^3/h]
%P = 35000;                                     % Environmental pressure [Pa]

S = Sin*ones(N,1);
Sm = S;
aw = ones(N,1);
Nub = 14*ones(N,1);
    %% Calculations
for o = 1:30
    %%% Density
    rhob = zeros(N,1);
    for i = 1:N
        a1 = 999.9; a2 =0.02034; a3 = -6.162e-3; a4 = 2.261e-5; a5 = -4.657e-8;
        b1 = 802; b2=-2.001; b3 =0.01677; b4 = -3.060e-5; b5 = -1.613e-5;
        rhob(i) = a1+a2*y(i)+a3*y(i)^2+a4*y(i)^3+a5*y(i)^4+b1*S(i)/1000+b2*S(i)/1000*y(i)+b3*S(i)/1000*y(i)^2+b4*S(i)/1000*y(i)^3+b5*(S(i)/1000)^2*y(i);
    end
    rhod = zeros(N,1);
    for i = 1:N
        Sd = S(1);
        a1 = 999.9; a2 =0.02034; a3 = -6.162e-3; a4 = 2.261e-5; a5 = -4.657e-8;
        b1 = 802; b2=-2.001; b3 =0.01677; b4 = -3.060e-5; b5 = -1.613e-5;
        rhod(i) = a1+a2*y(5*N+i)+a3*y(5*N+i)^2+a4*y(5*N+i)^3+a5*y(5*N+i)^4+b1*Sd/1000+b2*Sd/1000*y(5*N+i)+b3*Sd/1000*y(5*N+i)^2+b4*Sd/1000*y(5*N+i)^3+b5*(Sd/1000)^2*y(5*N+i);
    end
    
        % Specific heat
    cpb = zeros(N,1);
    for i = 1:N
        A = 5.328 - 9.76e-2*S(i) + 4.04e-4*S(i)^2;
        B = -6.913e-3 + 7.351e-4*S(i) - 3.15*10^-6 * S(i)^2;
        C = 9.6*10^-6 - 1.927*10^-6*S(i) + 8.23*10^-9 * S(i)^2;
        D = 2.5*10^-9 + 1.666*10^-9*S(i) - 7.125*10^-12 * S(i)^2;
        cpb(i) =1000*( A + B*(y(i)+273.15) + C*(y(i)+273.15)^2 + D*(y(i)+273.15)^3);
    end

    cpd = zeros(N,1);
    for i = 1:N
        Sd = S(1);
        A = 5.328 - 9.76e-2*Sd + 4.04e-4*Sd^2;
        B = -6.913e-3 + 7.351e-4*Sd - 3.15*10^-6 * Sd^2;
        C = 9.6*10^-6 - 1.927*10^-6*Sd + 8.23*10^-9 * Sd^2;
        D = 2.5*10^-9 + 1.666*10^-9*Sd - 7.125*10^-12 * Sd^2;
        cpd(i) =1000*( A + B*(y(5*N+i)+273.15) + C*(y(5*N+i)+273.15)^2 + D*(y(5*N+i)+273.15)^3);
    end
    
    cpdes = zeros(N,1);
    for i = 1:N
        cpdes(i) =1000*( 5.328 -6.913e-3*((y(3*N+i)+y(4*N+i))/2+273.15) + 9.6*10^-6*((y(3*N+i)+y(4*N+i))/2+273.15)^2 + 2.5*10^-9*((y(3*N+i)+y(4*N+i))/2+273.15)^3);
    end
    
    mfin = (Flowin/Nchan)*rhob(1)/3600;  % Massflow in the hot channel of one evelope      [kg/s]  
    mp = mfin;                           % Flow in the cold channel of one evelope     [kg/s]
    Mbf = zeros(N,1);
    for i = 1:N
        Mbf(i) = Wmd*dz*tchan*rhob(i)*sp_por + (1-sp_por)*Wmd*dz*tchan*946*1920/cpb(i);    % Mass in hot channel          [kg]
    end
    Mbp = zeros(N,1);
    for i = 1:N
        Mbp(i) = Wmd*dz*tchan*rhod(i)*sp_por + (1-sp_por)*Wmd*dz*tchan*946*1920/cpd(i);    % Mass in cold channel         [kg]
    end                      

    %%% Calculation flux
    Pca = zeros(N,1);
    h = 1;
    for i = 3*N+1:4*N
       Pca(h) = exp(23.1964-3816.44/(y(i)+227.02));
       h = h+1;
    end

    %%% Water activity
    m = zeros(N,1);
    for i = 1:N
        m(i) = Sm(i)/(MNaCl*1000); %molality of solution
        aw(i) = 1-0.03112*m(i)-0.001482*m(i)^2;
    end

    Pmf = zeros(N,1);
    h = 1;
    for i = N+1:2*N
        Pmf(h) = exp(23.1964-3816.44/(y(i)+227.02))*aw(h);
        h = h+1;
    end

    Pma = zeros(N,1);
    h = 1;
    for i = 2*N+1:3*N
        Pma(h) = exp(23.1964-3816.44/(y(i)+227.02));
        h = h+1;
    end

    %%% Calculations thermophysical properties
    % Viscosity
    etab = zeros(N,1);
    for i = 1:N
        A = 1.474e-3+1.5e-5*y(i)-3.927e-8*y(i)^2;
        B = 1.073e-5-8.5e-8*y(i)+2.23e-10*y(i)^2;
        etab(i) = exp(-10.7019+604.129/(139.18+y(i)))*(1+A*S(i)+B*S(i)^2);
    end

    etad = zeros(N,1);
    for i = 1:N
        Sd = S(1);
        A = 1.474e-3+1.5e-5*y(5*N+i)-3.927e-8*y(5*N+i)^2;
        B = 1.073e-5-8.5e-8*y(5*N+i)+2.23e-10*y(5*N+i)^2;
        etad(i) = exp(-10.7019+604.129/(139.18+y(5*N+i)))*(1+A*Sd+B*Sd^2);
    end

    % Conductivity of air
    k_airm = zeros(N,1);
    for i = 1:N
        k_airm(i) = 2.72*10^-3+7.77*10^-5*((y(N+i)+y(2*N+i))/2+273.15);
    end

    k_aira = zeros(N,1);
    for i = 1:N
        k_aira(i) = 2.72*10^-3+7.77*10^-5*((y(2*N+i)+y(3*N+i))/2+273.15);
    end
    % Conductivity of water
    k_watb = zeros(N,1);
    for i = 1:N
        k_watb(i) = 10^(log10(240+0.0002*S(i))+0.434*(2.3-(345.5+0.037*S(i))/(y(i)+273.15))*(1-(y(i)+273.15)/(647+0.03*S(i)))^0.333)/1000;
    end

    k_watd = zeros(N,1);
    for i = 1:N
        Sd = S(1);
        k_watd(i) = 10^(log10(240+0.0002*Sd)+0.434*(2.3-(345.5+0.037*Sd)/(y(5*N+i)+273.15))*(1-(y(5*N+i)+273.15)/(647+0.03*Sd))^0.333)/1000;
    end

    k_watj = zeros(N,1);
    for i = 1:N
        k_watj(i) = 10^(log10(240)+0.434*(2.3-345.5/((y(2*N+i)+y(3*N+i))/2+273.15))*(1-((y(2*N+i)+y(3*N+i))/2+273.15)/647)^0.333)/1000;
    end

    %%% Calculations membrane
    % heat transfer
    beta = zeros(N,1);
    for i = 1:N
       beta(i) = (k_mem_mat-k_airm(i))/((k_mem_mat+2*k_airm(i))); 
    end

    k_mem = zeros(N,1); hm = zeros(N,1);
    for i = 1:N
        k_mem(i) = 0.93*k_airm(i)*(1+2*beta(i)*(1-por_mem))/(1-beta(i)*(1-por_mem));      % Conductivity of the membrane
        hm(i) = k_mem(i)/t_mem;                                  % Heat transfer membrane
    end

    % mass transfer
    Dk = zeros(N,1);
    for i = 1:N
        Dk(i) = 2/3*r_mem*(8*R*((y(N+i)+y(2*N+i))/2+273.15)/(pi*M))^0.5;
    end

    %%% Partial pressures
    Pemf = zeros(N,1);
    yaf = zeros(N,1);
    for i = 1:N
       yaf(i) = (P-Pmf(i))/P;
       Pemf(i) = P;
       if yaf(i)<0
           Pemf(i) = Pmf(i)+1;
           yaf(i) = (Pemf(i)-Pmf(i))/Pemf(i);
       end
    end
    
    Pema = zeros(N,1);
    yap = ones(N,1);
    for i = 1:N
       yap(i) = (P-Pma(i))/P;
       Pema(i) = P;
       if yap(i)<0
           Pema(i) = Pma(i)+1;
           yap(i) = (Pema(i)-Pma(i))/Pema(i);
       end
    end
    
    Peca = zeros(N,1);
    for i = 1:N
        Peca(i) = P; 
        if  P<Pca(i)
            Peca(i) = Pca(i)+1;
        end
    end

    %%% Calculation permeability
    K0 = 2*por_mem*r_mem/(3*tor_mem);
    K1 = por_mem/tor_mem;
    
    Dwa = zeros(N,1);
    for i = 1:N
        Dwa(i) = K1*1.895e-5*((y(N+i)+y(2*N+i))/2+273.15)^2.072;  % Molecular diffusion
    end

    vmol = zeros(N,1);
    for i = 1:N
       vmol(i) = (8*R*((y(N+i)+y(2*N+i))/2+273.15)/(pi*M))^0.5; 
    end
    Dek = zeros(N,1);
    for i = 1:N
       Dek(i)= K0*vmol(i); 
    end

    mefrepat = zeros(N,1);
    for i = 1:N
       mefrepat(i) = kb*((y(N+i)+y(2*N+i))/2+273.15)/(pi*P*de^2*2^0.5);
    end

    Kn = zeros(N,1);
    for i = 1:N
       Kn(i) = mefrepat(i)/(2*r_mem); 
    end

    Cm = zeros(N,1);
    for i = 1:N
       Cm(i) = Dwa(i)*(1+Kn(i))/(t_mem*R*((y(N+i)+y(2*N+i))/2+273.15))*log((Dek(i)*yap(i)+Dwa(i)*(1+Kn(i)))/(Dek(i)*yaf(i)+Dwa(i)*(1+Kn(i))))*1/(Pmf(i)-Pma(i));
    end
    
    Dag = zeros(N,1);
    for i = 1:N
        Dag(i) = 4.46e-6*((y(2*N+i)+y(3*N+i))/2+273.15)^2.334;
    end
    
    Ca = zeros(N,1);
    for i = 1:N
       Tav = (y(2*N+i)+y(3*N+i))/2;
       Ca(i) = Dag(i)/(t_airgap*R*(Tav+273.15))*log((P-Pca(i))/(P-Pma(i)))*1/(Pma(i)-Pca(i));
    end
    
    Cnonflud = zeros(N,1);
    for i = 1:N
       Cnonflud(i) = (Ca(i)^-1+Cm(i)^-1)^-1;
    end
    
    %%% Calculation flux
    wat = 0.3;
    for j = 1:50
        J = zeros(N,1);
        for i = 1:N
            J(i) = (wat*Cm(i)*(Pmf(i)-Pma(i))+(1-wat)*Cnonflud(i)*(Pmf(i)-Pca(i)));
        end
        wat = 0.1676*log(sum(J*3600))^0.571+0.02811;
        if wat<0
            wat = 0;
        elseif imag(wat)>0
            wat = 0;
        end
    end
    
    % effective diameter
    dh = 4*sp_por/(2/tchan+((1-sp_por)*4/(tchan/2)));
    
        %%% Calculation flow reduction
    mf = zeros(N,1);
    mf(1) = mfin;
    for i = 2:N
        mf(i) = mf(i-1)-J(i-1)*Am*2;
    end

    %%% Calculations hot channel
    % effective speed
    v_bef = zeros(N,1);
    for i = 1:N
       v_bef(i) = mf(i)/(rhob(i)*tchan*Wmd*sp_por); 
    end

    % Reynolds nuber
    Reb = zeros(N,1);
    for i = 1:1:N
        Reb(i) = rhob(i)*v_bef(i)*dh/etab(i);
    end
    % Nusselts number
    Cnu = 0.19;
    Bnu = 0.68;
    for i = 1:N
        Nub(i) = Cnu*Reb(i)^Bnu;
    end
    
    Dnaclb = zeros(N,1);
    Dnaclm = zeros(N,1);
    etam = zeros(N,1);
    Scb = zeros(N,1); %Schimdt number
    Scm = zeros(N,1); %Schimdt number
    Sh = zeros(N,1);
    K = zeros(N,1);
    for i = 1:N
        A = 1.474e-3+1.5e-5*y(i)-3.927e-8*y(i)^2;
        B = 1.073e-5-8.5e-8*y(i)+2.23e-10*y(i)^2;
        etam(i) = exp(-10.7019+604.129/(139.18+y(i+N)))*(1+A*Sm(i)+B*Sm(i)^2);
        
        Dnaclb(i) = 117.3e-18*(2.26*M*1000)^0.5*(y(i)+273.15)/(etab(i)*(0.04838)^0.6);
        Dnaclm(i) = 117.3e-18*(2.26*M*1000)^0.5*(y(i+N)+273.15)/(etam(i)*(0.04838)^0.6);
        
        Scb(i) = etab(i)/Dnaclb(i);
        Scm(i)= etam(i)/Dnaclm(i);
        Sh(i) = Nub(i)*Scb(i)^0.13*(Scb(i)/Scm(i))^0.25;
        K(i) = Sh(i)*((Dnaclb(i)+Dnaclm(i))/2)/dh;
    end
    
    for i = 1:N
       Sm(i) = S(i)*exp(J(i)/(rhod(i)*K(i))); 
    end
    
    hbf = zeros(N,1);
    for i = 1:N
        hbf(i) = Nub(i)*k_watb(i)/dh;
    end

    %%% Calculations air gag
    % heat transfer
    ha = zeros(N,1);
    for i = 1:N
        ha(i) = (sa_por*(k_aira(i)*(1-wat)+wat*k_watj(i))+(1-sa_por)*k_sp)/(t_airgap-t_cond);
    end
    
    hcon = zeros(N,1);
    for i = 1:N
        hcon(i) = k_watj(i)/t_cond;
    end
    
    hceff = zeros(N,1);
    for i = 1:N
        hceff(i) = (hcon(i)^-1+hc^-1)^-1; 
    end
    
    %%% Calculations cold channel
    v_pef = Flowin/(3600*Nchan*tchan*Wmd*sp_por);

    % Reynolds nuber
    Rep = zeros(N,1);
    for i = 1:N
        Rep(i) = rhod(i)*v_pef*dh/etad(i);
    end

    Nup = zeros(N,1);
    for i = 1:N
        Nup(i) = Cnu*Rep(i)^Bnu;
    end
    hbpf = zeros(N,1);
    for i = 1:N
        hbpf(i) = Nup(i)*k_watd(i)/dh;
    end
    %% D formula's
    D1 = zeros(N,N);
    for i = 1:1:N
       D1(i,i) = -Lmd*mf(i)/(Mbf(i)*dz)-2*Lmd*Wmd/(Mbf(i)*cpb(i))*(hbf(i)+J(i)*cpdes(i)); 
    end
    
    h = 2;
    for i = 2:1:N
       D1(i,i-1) = Lmd*mf(h)/(Mbf(h)*dz);
       h = h+1;
    end
    
    D2 = zeros(N,N);
    for i =1:N
       D2(i,i) = 2*Lmd*Wmd/(Mbf(i)*cpb(i))*(hbf(i)+J(i)*cpdes(i));
    end

    D3 = zeros(N,N);
    for i = 1:1:N
        D3(i,i) = 2*Lmd*Wmd*hbpf(i)/(Mbp(i)*cpd(i));
    end

    D4 = eye(N,N);
    for i = 1:1:N
        D4(i,i) = -Lmd*mp/(Mbp(i)*dz)-2*Lmd*Wmd*hbpf(i)/(Mbp(i)*cpd(i));
    end

    h=2;
    for i = 2:1:N
       D4(i-1,i) = Lmd*mp/(Mbp(h)*dz); 
       h = h+1;
    end
    %% Z formula's
    Z1 = zeros(N,N);
    Z2 = zeros(N,N);
    Z3 = zeros(N,N);
    Z4 = zeros(N,N);
    Z5 = zeros(N,N);
    Z6 = zeros(N,N);
    Z7 = zeros(N,N);
    Z8 = zeros(N,N);
    Z9 = zeros(N,N);
    Z10 = zeros(N,N);
    Z11 = zeros(N,N);
    Z12 = zeros(N,N);
    for i = 1:N
        Z1(i,i) = hbf(i)+J(i)*cpdes(i)*Am;
        Z2(i,i) = -hbf(i)-2*J(i)*cpdes(i)*Am-hm(i);
        Z3(i,i) = hm(i)+J(i)*cpdes(i)*Am;

        Z4(i,i) = hm(i)+J(i)*cpdes(i)*Am;
        Z5(i,i) = -hm(i)-2*J(i)*cpdes(i)*Am-ha(i);
        Z6(i,i) = ha(i)+J(i)*cpdes(i);

        Z7(i,i) = ha(i)+J(i)*cpdes(i)*Am;
        Z8(i,i) = -ha(i)-J(i)*cpdes(i)*Am-hceff(i);
        Z9(i,i) = hceff(i);
        Z10(i,i) = hceff(i);
        Z11(i,i) = -hceff(i)-hbpf(i);
        Z12(i,i) = hbpf(i);
    end

    %% F matrix
    F = [D1,D2,zeros(N,N),zeros(N,N),zeros(N,N),zeros(N,N);
         Z1,Z2,Z3,zeros(N,N),zeros(N,N),zeros(N,N);
         zeros(N,N),Z4,Z5,Z6,zeros(N,N),zeros(N,N);
         zeros(N,N),zeros(N,N),Z7,Z8,Z9,zeros(N,N);
         zeros(N,N),zeros(N,N),zeros(N,N),Z10,Z11,Z12;
         zeros(N,N),zeros(N,N),zeros(N,N),zeros(N,N),D3,D4];
    %% H matrix Heat of evaporation
    Hv1 = zeros(N,1);
    h = 1;
    for i = 2*N+1:3*N
        Hv1(h) = 1000*(2501.897149-2.407064037*y(i)+1.192217e-3*y(i)^2-1.5863e-5*y(i)^3);
        h = h+1;
    end

    Hv2 = zeros(N,1);
    h = 1;
    for i = 3*N+1:4*N
       Hv2(h)= wat*1000*(2501.897149-2.407064037*y(i)+1.192217e-3*y(i)^2-1.5863e-5*y(i)^3);
       h = h+1;
    end
    
    Hv3 = zeros(N,1);
    h = 1;
    for i = 4*N+1:5*N
       Hv3(h) = (1-wat)*1000*(2501.897149-2.407064037*y(i)+1.192217e-3*y(i)^2-1.5863e-5*y(i)^3);
       h = h+1;
    end

    H = zeros(6*N,1);
    h = 1;
    for i = 2*N+1:3*N
       H(i) = -J(h)*Hv1(h)*Am;
       h = h+1;
    end

    h = 1;
    for i = 3*N+1:4*N
       H(i) = J(h)*Hv2(h)*Am; 
    end
    
    h = 1;
    for i = 4*N+1:5*N
       H(i) = J(h)*Hv3(h)*Am;
       h = h+1;
    end

    for i = 2:N
       S(i) = S(i-1)* mf(i-1)/mf(i);
    end
end
%% q0 boundary conditions matrix

q0 = zeros(6*N,1);
q0(1) = Lmd*mf(1)/(Mbf(1)*dz)*Tfin;
q0(6*N) = Lmd*mp/(Mbp(N)*dz)*Tpin;

%% Formule van ydot
ydot = F*y+H+q0;
