%%% 1D gas-phase combustion sub-model

% solves the ODE system with the functions using Strang splitting: 
% y1' = f1(t,y1)    -- the transport part in equations 15 & 16 in methods
% y2' = f2(t,y2)	-- the reaction part in equations 15 & 16 in methods


% load pyrolysis model solution

load ('pyrolysis_data.mat');

% load species definitions
load ('solid_kinetics_data.mat');

%%%%%%%%%%%%%%%% create the gas object %%%%%%%%%%%%%%%%%%%%%%%%
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties
%
% The gas phase will be taken from the definition of phase 'gas' in
% input file 'gas_kinetics_data.cti,' which contains the kinetics model for
% gas-phase combustion

gas = Solution('gas_kinetics_data.cti');

% mesh setup

Mesh.Jnodes = 5;  % grid size
domain_height = .013; % [m]
Mesh.dz = domain_height/Mesh.Jnodes;
Mesh.a = 10e-2^2;  % sample cross-sectional area [m2]
Mesh.dv = Mesh.a * Mesh.dz;
Mesh.dx = 10e-2;

% solver options
dt =.001;
nstep = 15000;
options = odeset('RelTol',1.e-4,'AbsTol',1e-7,'NonNegative',1,'MaxStep',.01);

%set initial conditions
time = 0;
t1 = zeros(nstep+1,1); % time

nsp = nSpecies(gas);
MW = molecularWeights(gas);
species(47)={'N2'};

T0 = zeros(Mesh.Jnodes,1)+300; T0(end)=1273;
yj0 = zeros(nsp,Mesh.Jnodes);

set(gas,'Temperature',300,'Pressure',101325.0,'MassFractions', 'Ar:1,O2:20,N2:78');
cpa = cp_mass(gas); da = thermalConductivity(gas)/(cpa*density(gas)); 
ya = massFractions(gas);rhoa = density(gas);
for i=1:Mesh.Jnodes
      yj0(:,i)= ya;
end
y0 = transpose([yj0(:); T0(:)]);

% solution array
y = zeros(nstep+1,length(y0)); 

t2 = transpose(0:.01:nstep*dt);
y(1,:) = y0;

s = length(t2);
x = zeros(nsp,s);
g_index = [3 4 5 6 7 8 9 10 11 12 13 14 20 21 29 30 31 47];
ye1 = ye(1:s,[1:12 14:15 17:19 23]); x1=ye1; j01=j0; 
Ts = Ts(1:s,:); j01 = j01(1:s,:);


for i=2:s(1)
   x1(i,:) = ye1(i,:)./sum(ye1(i,:));
   j01(i) = j0(i)*sum(ye(i,[1:12 14:15 17:19 23]));
	if j0(i)-j0(i-1) > .01
		j01(i) = j0(i-1);
	end
end

species1 = string(species(g_index));
species1(4)='C6H6O3'; species1(16)='C9H10O2'; species1(17)='C6H5OH'; 
for j=1:s(1)
    for i=1:length(g_index)
        x(speciesIndex(gas,char(species1(i))),j) = x1(j,i);
    end
end


% time integration 
for i=1:nstep
    Ts1 = interp1(t2,Ts,time,'linear');
    j011 = interp1(t2,j01,time,'linear');
    yin = interpy(time,t2,x);

    [~,a] = ode15s(@(t,y1) yprime2(t,y1,Mesh,Ts1,j011,yin, i),[time time+dt/2],y(i,:),options);
	temp = a(end,:);
    temp(temp<0)=1e-30;

    [~,b] = ode15s(@(t,y2) yprime1(t,y2,Mesh,Ts1,j011,yin),[time time+dt],temp,options);
	temp = b(end,:);
    temp(temp<0)=1e-30;

    [~,c] = ode15s(@(t,y1) yprime2(t,y1,Mesh,Ts1,j011,yin, i),[time+dt/2 time+dt],temp,options);
	temp = c(end,:);
    temp(temp<0)=1e-30;
    y(i+1,:) = temp;
        
    time = time + dt
    t1(i+1) = time;
    y(i+1,end)
	save results.mat y t1
end


global MW cpa da ya nsp gas rhoa 

%ODE function y1'

function [dydt] = yprime1(t,y,Mesh,Ts,j0,yin)

global cpa da ya nsp gas rhoa 
    
    T = y(nsp*Mesh.Jnodes+1:end);
    yj = zeros(nsp,Mesh.Jnodes);
    
    cp = zeros(Mesh.Jnodes,1);
    cps = zeros(Mesh.Jnodes,1);
    
    Tprime = zeros(Mesh.Jnodes,1);
    d = zeros(Mesh.Jnodes,1);
    ds = zeros(Mesh.Jnodes,1);
    
    dyjdt = zeros(nsp,Mesh.Jnodes);
    
    
    for i=1:Mesh.Jnodes
        yj(:,i)=y(nsp*(i-1)+1:nsp*(i-1)+nsp);
        
        set(gas,'Temperature',T(i),'Pressure',101325.0,'MassFractions', yj(:,i));
        
        cp(i) = cp_mass(gas);
        
        cps(i) = (cp(i)+cpa)/2;
        d(i) = thermalConductivity(gas)/(rhoa*cp(i));
        ds(i) = d(i).*da/(da+(d(i)-da)/2);
        
    end
    
    if Mesh.Jnodes ==1
        jxw = -rhoa*ds(1).*(yj(:,1)-ya)./(Mesh.dx/2);
        jxe = -rhoa*ds(1).*(ya-yj(:,1))./(Mesh.dx/2);
        Tprime(1) = ( - j0*cp(1)*(T(1)-Ts)/Mesh.dz...
               +4*rhoa*ds(1)*cps(1)*(300-T(1))/(Mesh.dx^2))/(rhoa*cp(1));

        dyjdt(:,1) = ( - j0*(yj(:,1)-yin)./Mesh.dz -(jxe-jxw)/Mesh.dx)/rhoa;
    end
    
    if Mesh.Jnodes ==2
        
        
        jxw = -rhoa*ds(1).*(yj(:,1)-ya)./(Mesh.dx/2);
        jxe = -rhoa*ds(1).*(ya-yj(:,1))./(Mesh.dx/2);
       
        
        Tprime(1) = ( - j0*cp(1)*(T(1)-Ts)/Mesh.dz...
                   +4*rhoa*ds(1)*cps(1)*(300-T(1))/(Mesh.dx^2))/(rhoa*cp(1));
        dyjdt(:,1) = ( - j0*(yj(:,1)-yin)./Mesh.dz -(jxe-jxw)/Mesh.dx)/rhoa;
        jxw = -rhoa*ds(end).*(yj(:,end)-ya)./(Mesh.dx/2);
        jxe = -rhoa*ds(end).*(ya-yj(:,end))./(Mesh.dx/2);
        
        Tprime(end) = ( - j0*cp(end)*(T(end)-T(1))/(Mesh.dz)...
                   +4*rhoa*ds(end)*cps(end)*(300-T(end))/(Mesh.dx^2))/(rhoa*cp(end));
        dyjdt(:,end) = ( - j0*(yj(:,end)-yj(:,end-1))./Mesh.dz -(jxe-jxw)/Mesh.dx)/rhoa;
    end
        

    if Mesh.Jnodes>2
        
        jxw = -rhoa*ds(1).*(yj(:,1)-ya)./(Mesh.dx);
        jxe = -rhoa*ds(1).*(ya-yj(:,1))./(Mesh.dx);
        
        Tprime(1) = ( - j0*cp(1)*(T(1)-Ts)/Mesh.dz...
                   +2*rhoa*ds(1)*cps(1)*(300-T(1))/(Mesh.dx^2))/(rhoa*cp(1));
        dyjdt(:,1) = ( - j0*(yj(:,1)-yin)./Mesh.dz -(jxe-jxw)/Mesh.dx)/rhoa;
        
        for i=2:Mesh.Jnodes-1
            jxw = -rhoa*ds(i).*(yj(:,i)-ya)./(Mesh.dx);
            jxe = -rhoa*ds(i).*(ya-yj(:,i))./(Mesh.dx);
            
            Tprime(i) = ( - j0*cp(i)*(T(i)-T(i-1))/Mesh.dz...
                   +2*rhoa*ds(i)*cps(i)*(300-T(i))/(Mesh.dx^2))/(rhoa*cp(i));
            dyjdt(:,i) = ( - j0*(yj(:,i)-yj(:,i-1))./Mesh.dz -(jxe-jxw)/Mesh.dx)/rhoa;
        end
        jxw = -rhoa*ds(end).*(yj(:,end)-ya)./(Mesh.dx);
        jxe = -rhoa*ds(end).*(ya-yj(:,end))./(Mesh.dx);
        Tprime(end) = ( - j0*cp(end)*(T(end)-T(end-1))/Mesh.dz...
                   +2*rhoa*ds(end)*cps(end)*(300-T(end))/(Mesh.dx^2))/(rhoa*cp(end));
        dyjdt(:,end) = ( - j0*(yj(:,end)-yj(:,end-1))./Mesh.dz -(jxe-jxw)/Mesh.dx)/rhoa;
        
    end
             
    Tprime(end)=0;    
    dydt = [dyjdt(:); Tprime(:)];
    
end


%ODE function y2'

function [dydt] = yprime2(t,y,Mesh,Ts,j0,yin, k)

global MW cpa nsp gas rhoa 
    
    T = y(nsp*Mesh.Jnodes+1:end);
    yj = zeros(nsp,Mesh.Jnodes);
    wdot = zeros(nsp,Mesh.Jnodes);
    cp = zeros(Mesh.Jnodes,1);
    cps = zeros(Mesh.Jnodes,1);
    Tprime = zeros(Mesh.Jnodes,1);
    dh = zeros(nsp,Mesh.Jnodes);
    dyjdt = zeros(nsp,Mesh.Jnodes);
        
    for i=1:Mesh.Jnodes
        yj(:,i)=y(nsp*(i-1)+1:nsp*(i-1)+nsp);
        set(gas,'Temperature',T(i),'Pressure',101325.0,'MassFractions', yj(:,i));

        if T(i)>= 600
            dh(:,i) = enthalpies_RT(gas).*gasconstant.*T(i);
            wdot(:,i) = netProdRates(gas); %kmol/m3
        end
        cp(i) = cp_mass(gas);
        
        
        cps(i) = (cp(i)+cpa)/2;
        
    end
      
    Tprime(1) = (-sum(wdot(:,1).*dh(:,1)) )/(rhoa*cp(1));

    dyjdt(:,1) = (wdot(:,1).*MW )/rhoa;

    if Mesh.Jnodes>1
        for i=2:Mesh.Jnodes
            Tprime(i) = (-sum(wdot(:,i).*dh(:,i)))/(rhoa*cp(i));
            dyjdt(:,i) = (wdot(:,i).*MW )/rhoa;
        end
    end        
        
    dydt = [dyjdt(:); Tprime(:)];
    
end

% function to evaluate mass fractions at intlet at desired time

function yinter = interpy(tt,xx,yy)
    global nsp
    yinter = zeros(nsp,1);
    for j=1:nsp
        yinter(j) = interp1(xx,yy(j,:),tt,'nearest');
    end
end
   
     
