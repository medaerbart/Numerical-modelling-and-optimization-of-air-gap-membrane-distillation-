function ydot = versie2_cal(t,y,N,Tpin,Tfin,Flowin,P,Sin,Nchan,Lmd)
%%% This fill will not work when you run it because the parameters are not
%%% included in this file. This file is only to show how the code is
%%% structured.

S = Sin*ones(N,1);
Sm = S;
aw = zeros(N,1);
m = zeros(N,1);
for i = 1:N
    m(i) = Sm(i)/(MNaCl*1000); %molality of solution
    aw(i) = 1-0.03112*m(i)-0.001482*m(i)^2;
end

Nub = 14*ones(N,1);

%% Calculations
for o = 1:30

    %%% Density
    rhob = zeros(N,1);
    h = 1;
    for i = 2:2:2*N
        a1 = 999.9; a2 =0.02034; a3 = -6.162e-3; a4 = 2.261e-5; a5 = -4.657e-8;
        b1 = 802; b2=-2.001; b3 =0.01677; b4 = -3.060e-5; b5 = -1.613e-5;
        rhob(h) = a1+a2*y(i)+a3*y(i)^2+a4*y(i)^3+a5*y(i)^4+b1*S(h)/1000+b2*S(h)/1000*y(i)+b3*S(h)/1000*y(i)^2+b4*S(h)/1000*y(i)^3+b5*(S(h)/1000)^2*y(i);
        h = h+1;
    end
    rhod = zeros(N,1);
    h = 1;
    for i = 2*N+3:2:4*N+1
        Sd = S(1);
        a1 = 999.9; a2 =0.02034; a3 = -6.162e-3; a4 = 2.261e-5; a5 = -4.657e-8;
        b1 = 802; b2=-2.001; b3 =0.01677; b4 = -3.060e-5; b5 = -1.613e-5;
        rhod(h) = a1+a2*y(i)+a3*y(i)^2+a4*y(i)^3+a5*y(i)^4+b1*Sd/1000+b2*Sd/1000*y(i)+b3*Sd/1000*y(i)^2+b4*Sd/1000*y(i)^3+b5*(Sd/1000)^2*y(i);
        h = h+1;
    end

    cpdes = zeros(N,1);
    for i = 1:N
        cpdes(i) =1000*( 5.328 -6.913e-3*((y(5*N+4+i)+y(6*N+4+i))/2+273.15) + 9.6*10^-6*((y(5*N+4+i)+y(6*N+4+i))/2+273.15)^2 + 2.5*10^-9*((y(5*N+4+i)+y(6*N+4+i))/2+273.15)^3);
    end
    
    mfin = (Flowin/Nchan)*rhob(1)/3600;  % Massflow in the hot channel      [kg/s]  
    mp = mfin;                           % Flow in the cold channel     [kg/s]

    Mbf = zeros(N,1);
    for i = 1:N
        Mbf(i) = Lmd*Wmd/N*tchan*rhob(i)*sp_por;    % Mass in hot channel          [kg]
    end
    Mbp = zeros(N,1);
    for i = 1:N
        Mbp(i) = Lmd*Wmd/N*tchan*rhod(i)*sp_por;% Mass in cold channel         [kg]
    end  
    
    %%% Calculation flux
    Pca = zeros(N,1);
    h = 1;
    for i = 6*N+5:7*N+4
       Pca(h) = exp(23.1964-3816.44/(y(i)+227.02));
       h = h+1;
    end
    
    %%% Water activity
    for i = 1:N
        m(i) = Sm(i)/(MNaCl*1000); %molality of solution
        aw(i) = 1-0.03112*m(i)-0.001482*m(i)^2;
    end
    
    Pmf = zeros(N,1);
    h = 1;
    for i = 4*N+5:5*N+4
        Pmf(h) = exp(23.1964-3816.44/(y(i)+227.02))*aw(h);
        h = h+1;
    end
    
    Pma = zeros(N,1);
    h = 1;
    for i = 5*N+5:6*N+4
        Pma(h) = exp(23.1964-3816.44/(y(i)+227.02));
        h = h+1;
    end
    
    % Specific heat
    cpb = zeros(N,1);
    h = 1;
    for i = 2:2:2*N
        A = 5.328 - 9.76e-2*S(h) + 4.04e-4*S(h)^2;
        B = -6.913e-3 + 7.351e-4*S(h) - 3.15*10^-6 * S(h)^2;
        C = 9.6*10^-6 - 1.927*10^-6*S(h) + 8.23*10^-9 * S(h)^2;
        D = 2.5*10^-9 + 1.666*10^-9*S(h) - 7.125*10^-12 * S(h)^2;
        cpb(h) =1000*( A + B*(y(i)+273.15) + C*(y(i)+273.15)^2 + D*(y(i)+273.15)^3);
        h = h+1;
    end

    cpd = zeros(N,1);
    h = 1;
    for i = 2*N+3:2:4*N+1
        Sd = S(1);
        A = 5.328 - 9.76e-2*Sd + 4.04e-4*Sd^2;
        B = -6.913e-3 + 7.351e-4*Sd - 3.15*10^-6 * Sd^2;
        C = 9.6*10^-6 - 1.927*10^-6*Sd + 8.23*10^-9 * Sd^2;
        D = 2.5*10^-9 + 1.666*10^-9*Sd - 7.125*10^-12 * Sd^2;
        cpd(h) =1000*( A + B*(y(i)+273.15) + C*(y(i)+273.15)^2 + D*(y(i)+273.15)^3);
        h = h+1;
    end
    
    % Viscosity
    etab = zeros(N,1);
    h = 1;
    for i = 2:2:2*N
        A = 1.474e-3+1.5e-5*y(i)-3.927e-8*y(i)^2;
        B = 1.073e-5-8.5e-8*y(i)+2.23e-10*y(i)^2;
        etab(h) = exp(-10.7019+604.129/(139.18+y(i)))*(1+A*S(h)+B*S(h)^2);
        h = h+1;
    end

    h = 1;
    etad = zeros(N,1);
    for i = 2*N+3:2:4*N+1
        Sd = S(1);
        A = 1.474e-3+1.5e-5*y(i)-3.927e-8*y(i)^2;
        B = 1.073e-5-8.5e-8*y(i)+2.23e-10*y(i)^2;
        etad(h) = exp(-10.7019+604.129/(139.18+y(i)))*(1+A*Sd+B*Sd^2);
        h = h+1;
    end
    
    % Conductivity of air
    k_airm = zeros(N,1);
    for i = 1:N
        k_airm(i) = 2.72*10^-3+7.77*10^-5*((y(4*N+4+i)+y(5*N+4+i))/2+273.15);
    end
    
    k_aira = zeros(N,1);
    for i = 1:N
        k_aira(i) = 2.72*10^-3+7.77*10^-5*((y(5*N+4+i)+y(6*N+4+i))/2+273.15);
    end
    
    % Conductivity of water
    k_watb = zeros(N,1);
    h = 1;
    for i = 2:2:2*N
        k_watb(h) = 10^(log10(240+0.0002*S(h))+0.434*(2.3-(345.5+0.037*S(h))/(y(i)+273.15))*(1-(y(i)+273.15)/(647+0.03*S(h)))^0.333)/1000;
        h = h+1;
    end
    
    k_watd = zeros(N,1);
    h = 1;
    for i = 2*N+3:2:4*N+1
        Sd = S(1);
        k_watd(h) = 10^(log10(240+0.0002*Sd)+0.434*(2.3-(345.5+0.037*Sd)/(y(i)+273.15))*(1-(y(i)+273.15)/(647+0.03*Sd))^0.333)/1000;
        h = h+1;
    end
    
    k_watj = zeros(N,1);
    for i = 1:N
        k_watj(i) = 10^(log10(240)+0.434*(2.3-345.5/((y(5*N+4+i)+y(6*N+4+i))/2+273.15))*(1-((y(5*N+4+i)+y(6*N+4+i))/2+273.15)/647)^0.333)/1000;
    end
    
    %%% Calculations membrane
    % heat transfer
    beta = zeros(N,1);
    for i = 1:N
       beta(i) = (k_mem_mat-k_airm(i))/((k_mem_mat+2*k_airm(i))); 
    end
    
    k_mem = zeros(N,1); hm = zeros(N,1); Rm = zeros(N,1);
    for i = 1:N
        k_mem(i) = 0.93*k_airm(i)*(1+2*beta(i)*(1-por_mem))/(1-beta(i)*(1-por_mem));      % Conductivity of the membrane
        hm(i) = k_mem(i)/t_mem;  % Heat transfer membrane
        Rm(i) = hm(i)^-1;        % Resistance membrane
    end
    
    % mass transfer
    Dk = zeros(N,1);
    for i = 1:N
        Dk(i) = 2/3*r_mem*(8*R*((y(4*N+4+i)+y(5*N+4+i))/2+273.15)/(pi()*M))^0.5;
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
        Dwa(i) = K1*1.895e-5*((y(4*N+4+i)+y(5*N+4+i))/2+273.15)^2.072;                                % Molecular diffusion
    end
    
    vmol = zeros(N,1);
    for i = 1:N
       vmol(i) = (8*R*((y(4*N+4+i)+y(5*N+4+i))/2+273.15)/(pi*M))^0.5; 
    end
    Dek = zeros(N,1);
    for i = 1:N
       Dek(i)= K0*vmol(i); 
    end
    
    mefrepat = zeros(N,1);
    for i = 1:N
       mefrepat(i) =  kb*((y(4*N+4+i)+y(5*N+4+i))/2+273.15)/(pi*P*de^2*2^0.5);
    end
    
    Kn = zeros(N,1);
    for i = 1:N
       Kn(i) = mefrepat(i)/r_mem; 
    end
    
    Cm = zeros(N,1);
    for i = 1:N
       Cm(i) = Dwa(i)*(1+Kn(i))/(t_mem*R*((y(4*N+4+i)+y(5*N+4+i))/2+273.15))*log((Dek(i)*yap(i)+Dwa(i)*(1+Kn(i)))/(Dek(i)*yaf(i)+Dwa(i)*(1+Kn(i))))*1/(Pmf(i)-Pma(i));
    end
    
    Dag = zeros(N,1);
    for i = 1:N
        Dag(i) = 4.46e-6*((y(5*N+4+i)+y(6*N+4+i))/2+273.15)^2.334;
    end
    
    Ca = zeros(N,1);
    for i = 1:N
       Tav = (y(5*N+4+i)+y(6*N+4+i))/2;
       Ca(i) = Dag(i)/(t_airgap*R*(Tav+273))*log((P-Pca(i))/(P-Pma(i)))*1/(Pma(i)-Pca(i));
    end
    
    Cnonflud = zeros(N,1);
    for i = 1:N
       Cnonflud(i) = (Ca(i)^-1+Cm(i)^-1)^-1;
    end
    
    Dnaclb = zeros(N,1);
    Dnaclm = zeros(N,1);
    etam = zeros(N,1);
    Scb = zeros(N,1); %Schimdt number
    Scm = zeros(N,1); %Schimdt number
    Sh = zeros(N,1);
    K = zeros(N,1);   
    mf = zeros(N,1);
    mf(1) = mfin;
    v_bef = zeros(N,1);
    Reb = zeros(N,1);
    
    % Nusselts number
    Cnu = 0.19;
    Bnu = 0.68;
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
        for i = 1:N
            A = 1.474e-3+1.5e-5*y(2*i)-3.927e-8*y(2*i)^2;
            B = 1.073e-5-8.5e-8*y(2*i)+2.23e-10*y(2*i)^2;
            etam(i) = exp(-10.7019+604.129/(139.18+y(4*N+4+i)))*(1+A*Sm(i)+B*Sm(i)^2);
        
            Dnaclb(i) = 117.3e-18*(2.26*M*1000)^0.5*(y(2*i)+273.15)/(etab(i)*(0.04838)^0.6);
            Dnaclm(i) = 117.3e-18*(2.26*M*1000)^0.5*(y(4*N+4+i)+273.15)/(etam(i)*(0.04838)^0.6);
        
            Scb(i) = etab(i)/Dnaclb(i);
            Scm(i)= etam(i)/Dnaclm(i);
            Sh(i) = Nub(i)*Scb(i)^0.13*(Scb(i)/Scm(i))^0.25;
            K(i) = Sh(i)*((Dnaclb(i)+Dnaclm(i))/2)/dh;
        end
        for i = 1:N
            Sm(i) = S(i)*exp(J(i)/(rhod(i)*K(i))); 
        end
        for i = 2:N
            mf(i) = mf(i-1)-J(i-1)*Am*2;
        end
        for i = 1:N
            v_bef(i) = mf(i)/rhob(i)/(tchan*Wmd*sp_por); 
        end
        % effective speed
        for i = 2:N
            S(i) = S(i-1)* v_bef(i-1)/v_bef(i);
        end
    
        % Reynolds nuber
        for i = 1:1:N
            Reb(i) = rhob(i)*v_bef(i)*dh/etab(i);
        end
        for i = 1:N
            Nub(i) = Cnu*Reb(i)^Bnu;
        end
    end
end
    
    hf = zeros(N,1); Rf = zeros(N,1);
    for i = 1:N
        hf(i) = Nub(i)*k_watb(i)/dh;
        Rf(i) = (Am*hf(i))^-1;
    end
    
    % heat transfer
    ha = zeros(N,1); Ra = zeros(N,1);
    for i = 1:N
        ha(i) = (sa_por*(k_aira(i)*(1-wat)+wat*k_watj(i))+(1-sa_por)*k_sp)/(t_airgap-t_cond);
        Ra(i) = (Am*ha(i))^-1;
    end
    
    hcon = zeros(N,1);
    for i = 1:N
        hcon(i) = k_watj(i)/t_cond;
    end
    
    hceff = zeros(N,1);
    Rc = zeros(N,1);
    for i = 1:N
        hceff(i) = (hcon(i)^-1+hc^-1)^-1;
        Rc(i) = 1/(Am*hceff(i));
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
    hp = zeros(N,1); Rp = zeros(N,1);
    for i = 1:N
        hp(i) = Nup(i)*k_watd(i)/dh;
        Rp(i) = (Am*hp(i))^-1;
    end  
    
    Cbf = zeros(N,1); Cbp = zeros(N,1);
    for i = 1:N
       Cbf(i) = Mbf(i)*cpb(i)*sp_por+(1-sp_por)*Am*tchan*946*1920; 
       Cbp(i) = Mbp(i)*cpd(i)*sp_por+(1-sp_por)*Am*tchan*946*1920;
    end

    Rfz = zeros(N,1);Lf = zeros(N,1);Rpz = zeros(N,1);Lp = zeros(N,1);
    for i = 1:1:N
       Rfz(i) = 1/(mf(i)^2*cpb(i)^2*(Rf(i)/2+Rm(i)/1.5+Ra(i)/1.75+1.5*Rc(i)+Rp(i))*3); 
       Rpz(i) = 1/(mp^2*cpd(i)^2*(1.5*Rf(i)+Rm(i)/1.1+Ra(i)+Rc(i)+Rp(i))); %Heeft invloed op sinus en op de uiteindelijke waarde, hoe kleiner hoe lager

       Lf(i)  = Rfz(i)^2*Cbf(i)/4;
       Lp(i)  = Rpz(i)^2*Cbp(i)/4;
    end

    Rfz = [Rfz(1);Rfz];
    Rpz = [Rpz(1);Rpz];
    Lf = [Lf(1);Lf];
    Lp = [Lp(1);Lp];

    %% Alle matrixen
    Z1 = zeros(N,N); Z2 = zeros(N,N);Z4 = zeros(N,N); Z5 = zeros(N,N);Z6 = zeros(N,N);Z7 = zeros(N,N);Z8 = zeros(N,N); Z9 = zeros(N,N);Z10 = zeros(N,N);Z11 = zeros(N,N);
    for i = 1:N
       Z1(i,i) = -1/Rc(i)-1/Rp(i); 
       Z2(i,i) = 1/Rc(i);
       Z4(i,i) = 1/Rc(i);
       Z5(i,i) = -1/Ra(i)-1/Rc(i)-J(i)*Am*cpb(i);
       Z6(i,i) = 1/Ra(i)+ J(i)*Am*cpb(i);
       Z7(i,i) = 1/Ra(i) + J(i)*Am*cpb(i);
       Z8(i,i) = -1/Rm(i)-1/Ra(i);
       Z9(i,i) = 1/Rm(i);
       Z10(i,i) = 1/Rm(i)+J(i)*Am*cpb(i);
       Z11(i,i) = -1/Rf(i)-1/Rm(i)- 2*J(i)*Am*cpb(i);
    end

    Z3 = zeros(N,2*N+1);
    i = 1;
    for k = 2:2:2*N+1
        Z3(i,k) = 1/Rp(i);
        i = i+1;
    end

    Z12 = zeros(N,2*N+1);
    i = 1;
    for k = 2:2:2*N+1
        Z12(i,k) = 1/Rf(i)+J(i)*Am*cpb(i);
        i = i+1;
    end

    Zf1 = zeros(2*N+1,2);
    Zf1(2*N+1,1) = -1/Lf(1);

    I = eye(2,2);

    Af = zeros(2*N+1);
    g = 1;
    h = 1;
    for i = 1:2:2*N+1
       Af(i,g) = -Rfz(h)/Lf(h);
       g = g+2;
       h = h+1;
    end
    g = 2;
    h = 1;
    for i = 1:2:2*N
       Af(i,g) = -1/Lf(h);
       g = g+2;
       h = h+1;
    end
    g = 1;
    h = 1;
    for i = 2:2:2*N+1
       Af(i,g) = 1/Cbf(h);
       h = h+1;
       g = g+2; 
    end
    g = 2;
    h = 1;
    for i = 2:2:2*N+1
       Af(i,g) = -2/(Cbf(h)*Rf(h)) - 2*J(h)*Am*cpb(h)/Cbf(h);
       g = g+2;
       h = h+1;
    end
    g = 3;
    h = 1;
    for i = 2:2:2*N+1
       Af(i,g) = -1/Cbf(h);
       g = g+2;
       h = h+1;
    end
    g = 2;
    h = 1;
    for i = 3:2:2*N+1
       Af(i,g) = 1/Lf(h);
       g = g+2;
       h = h+1;
    end
    Af(1,2) = -1/(Lf(1)/2);
    Af(2*N+1,2*N) = 1/(Lf(1)/2);

    Zf2 = zeros(2*N+1,N);
    i = 1;
    for k = 2:2:2*N+1
        Zf2(k,i) = 2/(Cbf(i)*Rf(i));
        i = i+1;
    end

    Ap = zeros(2*N+1,2*N+1);

    g = 1;
    h = 1;
    for i = 1:2:2*N+1
       Ap(i,g) = -Rpz(h)/Lp(h);
       g = g+2;
       h = h+1;
    end
    g = 2;
    h = 1;
    for i = 1:2:2*N
       Ap(i,g) = -1/Lp(h);
       h = h+1;
       g = g+2;
    end
    g = 1;
    h = 1;
    for i = 2:2:2*N+1
       Ap(i,g) = 1/Cbp(h);
       g = g+2;
       h = h+1;
    end
    g = 2;
    h = 1;
    for i = 2:2:2*N+1
       Ap(i,g) = -2/(Cbp(h)*Rp(h));
       g = g+2;
       h = h+1;
    end
    g = 3;
    h = 1;
    for i = 2:2:2*N+1
       Ap(i,g) = -1/Cbp(h);
       g = g+2;
       h = h+1;
    end
    g = 2;
    h = 1;
    for i = 3:2:2*N+1
       Ap(i,g) = 1/Lp(h);
       g = g+2;
       h = h+1;
    end
    Ap(1,2) = -1/(Lp(1)/2);
    Ap(2*N+1,2*N) = 1/(Lp(N)/2);

    Zp1 = zeros(2*N+1,2);
    Zp1(1,2) = 1/(Lp(1)/2);

    Zp2 = zeros(2*N+1,N);
    i = 1;
    for k = 2:2:2*N+1
        Zp2(k,i) = 2/(Cbp(i)*Rp(i));
        i = i+1;
    end

    Tf0 = zeros(2,2*N+1);
    Tf0(1,2*N+1) = -1/(mf(1)*cpb(1));

    Tp0 = zeros(2,2*N+1);
    Tp0(2,1) = 1/(mp*cpd(1));

    %% F matrix

    F = [Af,zeros(2*N+1),Zf1,Zf2,zeros(2*N+1,N),zeros(2*N+1,N),zeros(2*N+1,N);
         zeros(2*N+1),Ap,Zp1,zeros(2*N+1,N),zeros(2*N+1,N),zeros(2*N+1,N),Zp2;
         Tf0,Tp0,I,zeros(2,N),zeros(2,N),zeros(2,N),zeros(2,N);
         Z12,zeros(N,2*N+1),zeros(N,2),Z11,Z10,zeros(N,N),zeros(N,N);
         zeros(N,2*N+1),zeros(N,2*N+1),zeros(N,2),Z9,Z8,Z7,zeros(N,N);
         zeros(N,2*N+1),zeros(N,2*N+1),zeros(N,2),zeros(N,N),Z6,Z5,Z4;
         zeros(N,2*N+1),Z3,zeros(N,2),zeros(N,N),zeros(N,N),Z2,Z1];

    %% B matrix
    B = zeros(8*N+4,2);
    B(1,1) = 1/(Lf(1)/2);
    B(4*N+2,2) = -1/(Lp(N)/2);
    B(4*N+3,2) = -1;
    B(4*N+4,1) = -1;
    Bu = B*[Tfin;Tpin];

    %% H matrix voor Hv
    H = zeros(8*N+4,1);
    g = 1;
    for i = (4*N+5):(5*N+4) 
        H(i) = -Am*J(g)*(2501.897149-2.407064037*y(i)+1.192217e-3*y(i)^2-1.5863e-5*y(i)^3)*1000;
        g = g+1;
    end
    g = 1;
    for i = (5*N+5):(6*N+4)
       H(i) =  wat*Am*J(g)*(2501.897149-2.407064037*y(i)+1.192217e-3*y(i)^2-1.5863e-5*y(i)^3)*1000;
       g = g+1;
    end
    
    g = 1;
    for i = (6*N+5):(7*N+4)
       H(i) =  (1-wat)*Am*J(g)*(2501.897149-2.407064037*y(i)+1.192217e-3*y(i)^2-1.5863e-5*y(i)^3)*1000;
       g = g+1;
    end

%% Formule van ydot
ydot = F*y+Bu+H;
