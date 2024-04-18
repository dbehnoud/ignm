%%% 1D pyrolysis model

% solves the two ODE system with the functions: 
% y1' = f1(t,y1)    -- equations 1 & 2 in methods
% y2' = f2(t,y2)	-- equations 10 & 11 in methods

% solution vectors y1 and y2 store the system variables: 

% y1 = [mass_species1_cell_1, mass_species2_cell_1, ... , mass_speciesM_cell_1, ...
% , mass_solid_species1_cell_2, mass_solid_species2_cell_2, ... , mass_solid_speciesM_cell_2, ... 
% , mass_solid_speciesM_cell_N, T_cell_1, T_cell_2, ... , T_cell_N, rho_cell_1, ...
% , solid_rho_cell_2, ... , solid_rho_cell_N]

% y2 = [gas_rho*g*phi_cell_1, ... , gas_rho*g*phi_cell_N, gas_rho*g*phi*y_species1_cell_1, ...
% , gas_rho*g*phi*y_gas_species2_cell_1, ..., gas_rho*g*phi*y_gas_speciesM_cell_1, ...
% , gas_rho*g*phi*y_gas_species1_cell_2, ... , gas_rho*g*phi*y_gas_speciesM_cell_N]



% load reaction rates parameters 
load ('kinetics_data.mat');


% mesh set-up

Mesh.Jnodes = 10; % number of cells
sample_height = 3.8e-2;
Mesh.dz = sample_height/(Mesh.Jnodes); 
Mesh.a = 3.8e-2^2; % cross-sectional area of each cell 
Mesh.dv = Mesh.a * Mesh.dz;

%%%% initial conditions  %%%%%%%%%%%%

% 'nsp' is total number of species
% 'gsp' is number of gas-phase species
% 'g_index' stores the indices of gas-phase species
% 's_index' stores the indices of solid-phase species
% MW stores moecular weight of species (Kg/mol)

T0 = zeros(Mesh.Jnodes,1) + 300;  %temperature (K)
mass0 = zeros(nsp,Mesh.Jnodes); % mass of solid species

%white pine
mass0(1,:) = 0.4254; % CELL
mass0(17,:) = 0.1927; % HCE
mass0(24,:) = 0.0998; % LIGH
mass0(25,:) = 0.0482 % LIGO
mass0(23,:) = 0.1658; % LIGC
mass0(38,:) = 0.0326; %TGL
mass0(37,:) = 0.0354; %CTANN
mass0(39,:) = 0.05; %moisture

rhos0 = zeros(Mesh.Jnodes,1) + 380*(1+mass0(39,1)); % initial solid density (Kg/m3)
sample_mass = Mesh.a*sample_height*rhos0(1);
mass0 = mass0./sum(mass0(s_index,1))*sample_mass./Mesh.Jnodes;
y10 = [mass0(:); T0(:); rhos0(:)]; % initial solution vector y1
yi0 = mass0(s_index,1)./sum(mass0(s_index,1)); % initial solid mass fraction


p0 = 1.013e5; %initial pressure (Pa)
yj0 = zeros(gsp,1);  % initial gas mass fraction
yj0(end) = 1;  % only N2 present
M = 1/sum(yj0./MW(g_index)); % average MW
R = 8.314; % gas constant
rhog0 = zeros(Mesh.Jnodes,1) + (p0)*M/(R*T0(1)); % initial gas density (Kg/m3)
rhogphi0 = rhog0*phii(yi0,rhos0(1)); % gas_rho*g*phi
rgpy0 = zeros(gsp,Mesh.Jnodes) + rhogphi0(1).*yj0; %gas_rho*g*phi*y_gas_species
y20 = [rhogphi0(:); rgpy0(:)]; % initial solution vector y2

% input radiative heat flux (W/m2)
qs = 40000; 


%%% variable initialization  %%%%%%%%%%%%%

dt =.01; % time step size
nstep = 20000; % number of time steps
time = 0;
t = zeros(nstep+1,1); 
t(1)= 0;

yy = zeros(nstep+1,length(y20));
yy1 = zeros(nstep+1,length(y10));
yy(1,:) = y20;
yy1(1,:) = y10;



%%%% time integration %%%%%%%%%%%%%

for i=1:nstep
    tspan = [t(i) t(i)+dt]
    [~,b] = ode113(@(t,y)yprime1(time,y,Mesh),tspan,yy1(i,:),options1); % equation 1
    [~,a] = ode113(@(t,y)yprime(time,y,Mesh,yy1(i,:)),tspan,yy(i,:),options2); % equation 2
	
	% this step ensures the mass fration values are non-negative
    temp = a(end,:);
    temp(temp<0)=1e-30;
	
    yy(i+1,:) = temp;
    yy1(i+1,:) = b(end,:);
    
    t(i+1) = t(i) + dt;
end
 

global ycoeff afac nfac ea istart qs g_index s_index  MW gsp nsp p0 yj0


% ODE function y2'

function [dydt] = yprime(t,yy,Mesh,yy1)

global ycoeff afac nfac ea istart s_index g_index MW gsp nsp p0 yj0

    wdot_mass = zeros(nsp,Mesh.Jnodes); % species mass production rate
    k = zeros(28,Mesh.Jnodes); % reaction rate coefficient
    m = zeros(nsp,Mesh.Jnodes); % mass
    phi = zeros(Mesh.Jnodes,1); % porosity
    kb = zeros(Mesh.Jnodes,1); % %thermal conductivity 
    e = zeros(Mesh.Jnodes,1); % emissivity
    rhos = zeros(Mesh.Jnodes,1); 
    drhosdt = zeros(Mesh.Jnodes,1); 
    mprime = zeros(nsp,Mesh.Jnodes);
    p = zeros(Mesh.Jnodes,1); % non-staggered pressure points
    flux = zeros(Mesh.Jnodes,1); %flux based on staggered vertical velocity points
    yj = zeros(gsp,Mesh.Jnodes); % gas-phase species mass fraction
    j = zeros(gsp,Mesh.Jnodes); % diffusive flux
    D3 = zeros(Mesh.Jnodes,1); %  diffusivity
    rhogphi = zeros(Mesh.Jnodes,1);
    drgpydt = zeros(gsp,Mesh.Jnodes);
    drhogphidt = zeros(Mesh.Jnodes,1);
    
    for i=1:Mesh.Jnodes
        m(:,i)=yy1(nsp*(i-1)+1:nsp*(i-1)+nsp);
        m(:,i)= m(:,i)./MW;
    end
    
    yi = zeros(length(s_index),Mesh.Jnodes);
    T = yy1(nsp*Mesh.Jnodes+1:(nsp+1)*Mesh.Jnodes);
    rhos(:) = yy1((nsp+1)*Mesh.Jnodes+1:(nsp+2)*Mesh.Jnodes);
    rhogphi(:) = yy(1:Mesh.Jnodes);
    rgpy = reshape(yy(Mesh.Jnodes+1:end),gsp,Mesh.Jnodes);
    for i=1:gsp
        yj(i,:)=rgpy(i,:)./transpose(rhogphi(:));
    end
    
    R = 8.314; % gas constant
	Kd = 1e-10; % porous solid permeability
    
    for i=1:Mesh.Jnodes
          
        % ycoeff: reaction stoichiometry
		% afac: Z_A parameter in eqn 9 in methods
		% nfac: n parameter in eqn 9 in methods
		% ea: E_A parameter in eqn 9 in methods
		% istart: index of the reacting species for each reaction
		
        yi(:,i) = m(s_index,i).*MW(s_index)./sum(m(s_index,i).*MW(s_index));
        phi(i) = phii(yi(:,i),rhos(i));
        k(:,i) = afac .*((T(i)).^nfac).* exp(-ea ./(R*T(i)));
        mprime(:,i) = ycoeff*(k(:,i).*m(istart,i)).*MW; 
        wdot_mass(:,i) = mprime(:,i)./ Mesh.dv;
        kb(i)= kba(T(i),yi(:,i), phi(i),rhos(i)); 
        e(i) = epsilon(yi(:,i),rhos(i),phi(i));
        M = 1/sum(yj(:,i)./MW(g_index)); 
        p(i) = rhogphi(i)/phi(i)*R*abs(T(i))/M-p0;
        D3(i) = .018829*sqrt(T(i)^3*(1/32+1/28))/((p(i)+p0)*5.061^2*.93);

    end 
    
     
     for i=1:Mesh.Jnodes-1
         
         flux(i) = -Kd*1/D3(i)*((p(i+1)-p(i))/Mesh.dz- rhogphi(i)/phi(i)*10*0);
        if flux(i)<0
            flux(i)=0;
        end
		
        for k=1:gsp
            D = D3(i+1)*D3(i)/(D3(i+1)+(D3(i)-D3(i+1))/2);
            rhophi = (rhogphi(i)+rhogphi(i+1))/2;
            j(k,i) = -D/MW(k)*rhophi*(yj(k,i+1)-yj(k,i))/Mesh.dz;
        end
     end
    
    flux(Mesh.Jnodes) = -Kd/D3(end)*((0-p(Mesh.Jnodes))/(Mesh.dz/2)...
        - rhogphi(Mesh.Jnodes)/phi(Mesh.Jnodes)*10*0);
    if flux(Mesh.Jnodes)<0
        flux(Mesh.Jnodes)=0;
    end
	
	yja = zeros(gsp,1); % surrounding gas mass fraction 
    yja(end)= 1; % only N2 present
	
    for ii=1:gsp
        j(ii,Mesh.Jnodes) = 0-.01*(yja(ii)-yj(ii,end));
        if j(ii,Mesh.Jnodes)<0
            j(ii,Mesh.Jnodes)=0;
        end
    end

    for i=2:Mesh.Jnodes-1
        drhogphidt(i) = sum(wdot_mass(g_index,i)) - (flux(i)-flux(i-1))/Mesh.dz;
        yfi = yj(:,i).*flux(i);
        yfii = yj(:,i-1).*flux(i-1);
        drgpydt(:,i) = wdot_mass(g_index,i) - (yfi-yfii)./Mesh.dz -(j(:,i)-j(:,i-1))/Mesh.dz;    
    end
    
	% bottom boundary
    drhogphidt(1) = sum(wdot_mass(g_index,1)) - flux(1)/Mesh.dz;
    drgpydt(:,1) = wdot_mass(g_index,1) - (yj(:,1)*flux(1))./Mesh.dz - j(:,1)./Mesh.dz;

	% top boundary
    drhogphidt(end) = sum(wdot_mass(g_index,end))-(flux(end)-flux(end-1))/Mesh.dz;
    yfi = yj(:,end)*flux(end);
    yfii = yj(:,end-1)*flux(end-1);
    drgpydt(:,end) = wdot_mass(g_index,end)-(yfi-yfii)./Mesh.dz -(j(:,end)-j(:,end-1))/Mesh.dz;
    
	for i=1:Mesh.Jnodes
        drhosdt(i) = - sum(wdot_mass(g_index,i));
    end

    dydt = [drhogphidt(:); drgpydt(:)];  
end


% ODE function y1'

function [dydt] = yprime1(t,yy,Mesh)

global ycoeff afac nfac ea istart s_index g_index qs MW deltah gsp nsp tempflux p0

    wdot_mass = zeros(nsp,Mesh.Jnodes);
    k = zeros(28,Mesh.Jnodes);
    m = zeros(nsp,Mesh.Jnodes);
    phi = zeros(Mesh.Jnodes,1);
    kb = zeros(Mesh.Jnodes,1);
    e = zeros(Mesh.Jnodes,1);
    rho_s_mass = zeros(Mesh.Jnodes,1);
    drhosdt = zeros(Mesh.Jnodes,1);
    mprime = zeros(nsp,Mesh.Jnodes);
    Tprime = zeros(Mesh.Jnodes,1);
	yi = zeros(length(s_index),Mesh.Jnodes);
    
    for i=1:Mesh.Jnodes
        m(:,i)=yy(nsp*(i-1)+1:nsp*(i-1)+nsp);
        m(:,i)= m(:,i)./MW;
    end
    
    T = yy(nsp*Mesh.Jnodes+1:(nsp+1)*Mesh.Jnodes);
    rho_s_mass(:) = yy((nsp+1)*Mesh.Jnodes+1:(nsp+2)*Mesh.Jnodes);
    
    R = 8.314; 
    sigma = 5.670374419e-8; 
	h =10; % heat transfer coefficient at top boundary
	tr=0;
    
    
    for i=1:Mesh.Jnodes
          
        yi(:,i) = m(s_index,i).*MW(s_index)./sum(m(s_index,i).*MW(s_index));
        phi(i) = phii(yi(:,i),rho_s_mass(i));
        k(:,i) = afac .*((T(i)).^nfac).* exp(-ea ./(R*T(i)));
        mprime(:,i) = ycoeff*(k(:,i).*m(istart,i)).*MW; 
        wdot_mass(:,i) = mprime(:,i)./ Mesh.dv;
        kb(i)= kba(T(i),yi(:,i), phi(i),rho_s_mass(i)); 
        e(i) = epsilon(yi(:,i),rho_s_mass(i),phi(i));
    end 
    
    for j=2:Mesh.Jnodes-1
       ddd = cp(T(j));
       deltah = q_srxns(T(end));
       Tprime(j) = (1/(Mesh.dz^2)*((kb(j)+kb(j+1))/2*(T(j+1)-T(j))+(kb(j)+kb(j-1))/2*(T(j-1)-T(j)))...
           +e(j)*tr/Mesh.Jnodes*qs/Mesh.dz+sum(abs(wdot_mass(istart,j)).*q_srxns(T(j))))...
           /(rho_s_mass(j).*sum(ddd(s_index).*yi(:,j))); 
    end
    
	% top boundary
	de = cp(T(end));
    c = sum(de(s_index).*yi(:,end));
	
    Tprime(Mesh.Jnodes)= (Mesh.a*(e(end)*qs*(1-tr)-h*(T(end)-300)-e(end)*sigma*(T(end)^4-300^4))...
        -Mesh.a*(kb(end)+kb(end-1))/2*(T(end)-T(end-1))/Mesh.dz...
     +Mesh.dv*sum(abs(wdot_mass(istart,end)).*q_srxns(T(end))))/(Mesh.dv*rho_s_mass(end)*c);
	 
	% bottom boundary
    d1 = cp(T(1));
	
    Tprime(1)=(Mesh.a*kb(1)/(Mesh.dz)*(T(2)-T(1))+Mesh.a*e(1)*tr/Mesh.Jnodes*qs...
        +Mesh.dv*sum(abs(wdot_mass(istart,1)).*q_srxns(T(1))))...
        /(Mesh.dv*rho_s_mass(1)*sum(d1(s_index).*yi(:,1)));
     
     
    for i=1:Mesh.Jnodes
        drhosdt(i) = - sum(wdot_mass(g_index,i));
    end

      
    dydt = [mprime(:); Tprime(:); drhosdt(:)];  
end


% function defining heat conductivity [W/m/K]
function kb = kba(T,yi,phi,rho_s_mass)
    
    global s_index MW
    k = zeros(length(s_index),1)+.17*(T/300)^.594;
    k(19)=.6;
    k(3)=.065*(T/300)^.435+5.670374419e-8*3.3e-3*(T)^3;
   
    s_density = [9.37;9.37;25;11.5;11.5;11.5;5.88;3.48;3.59;5.88;4;7.29;5.76;7.22;...
	5;1.67;55;0.0037;0.0058;0.0054;0.0807;0.01014;0.0051;0.0058].*MW(s_index)*1000; 
    yi(18:end)=0;
    kb = rho_s_mass*sum(yi.*k./(s_density))/(1-phi);
    
end

% function defining emissivity 
function e = epsilon(yi,rho_s_mass,phi)
    global s_index MW
    e = zeros(length(s_index),1)+0.757;
    e(3)=0.957; % char
    e(19)=.95; %H2O

    s_density = [9.37;9.37;25;11.5;11.5;11.5;5.88;3.48;3.59;5.88;4;7.29;5.76;7.22;5;1.67;...
        55;0.0037;0.0058;0.0054;0.0807;0.01014;0.0051;0.0058].*MW(s_index)*1000; 
    yi(18:end)=0;
    e = rho_s_mass*sum(yi.*e./(s_density))/(1-phi);
end 

% function defining heat capacity [J/kg/K]
function cp = cp(T)
global nsp
 
    cp = zeros(nsp,1)+(1.5+.001*T)*1000;
    cp(15)= (.7+.0035*T)*1000; %char
    cp(39) = 4188; %H2O

end


% function defining heat of reactions [J/kg of reactant]
function q_srxns = q_srxns(T)

    global ycoeff MW istart
    
    deltah = [-1300; 27100; 23200; -62700; -5000; -500; -42400; 17900; 12000;...
	-10300; 30700; 26000; -31100; -26100; 46200; -21100; -83600; 1300; 1300;...
	10100; -29100; -13400; 48600; 0; 0; 0; 0; 0]*4.184;
    q_srxns = deltah./MW(istart);
    q_srxns(28) = -2.41e6;
end

% function defining porosity 
function phi = phii(yi,rho_s_mass)
    
    global MW s_index
    
    s_density = [9.37;9.37;25;11.5;11.5;11.5;5.88;3.48;3.59;5.88;4;7.29;5.76;7.22;...
	5;1.67;55;0.0037;0.0058;0.0054;0.0807;0.01014;0.0051;0.0058]*1000;
    yi(18:end)=0;
    phi = 1-sum(yi./(s_density.*MW(s_index)))*rho_s_mass;
end
