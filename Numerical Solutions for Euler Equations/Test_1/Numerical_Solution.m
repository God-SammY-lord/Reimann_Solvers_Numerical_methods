%% 1D Shock Tube using Lax-Friedrich and Adjustable Time Stepping
close all;
clear;
clc;
%% Initialization of Parameters
N       = 100 ;             % Number of grid points
gamma   = 1.4 ;
endTime = 0.2 ;             % Physical End Time
CFL     = 0.5 ;
%% Grid Initialization
x       = linspace(0,1,N+1) ;
xc      = 0.5 * ( x(1:N) + x(2:N+1) ) ;
xc(1)   = 0 ;               % Xmin
xc(N)   = 1 ;               % Xmax
time    = 0 ;
%% Initial Conidtions
denistyR    = 0.125 ;       densityL    = 1.0 ;
pressureR   = 0.1   ;       pressureL   = 1.0 ;

rho     = zeros(N,1) ;
p       = zeros(N,1) ;
u       = zeros(N,1) ;

for i =  1:N
    if i<=N/2
        rho(i)  = densityL  ;
        p(i)    = pressureL ;        
    else
        rho(i)  = denistyR  ;
        p(i)    = pressureR ;
    end
end
e   = p/(gamma-1) + 0.5*rho.*u.*u ;
%%
new_rho = rho ;
new_u   = u   ;
new_e   = e   ;

while time <= endTime
    
    for i=2:N-1
        p(i)    = (gamma-1)*(e(i) - 0.5*rho(i)*u(i)*u(i)) ;
    end
    a       = sqrt(gamma*p./rho) ;
    lamda   = max(a) ;
    max_vel = max(u) ;
    
    dt      = CFL/N/(max_vel+lamda) ;  % adjustable Time Step
    time    = time + dt ;
    
    for i=2:N-1
        dx          = xc(i) - xc(i-1) ;
        
        mom_R       = rho(i+1)*u(i+1) ;     rho_R = rho(i+1) ;      u_R = u(i+1) ;      p_R = p(i+1) ;
        mom_P       = rho(i)*u(i)     ; 	rho_P = rho(i)   ;      u_P = u(i)   ;      p_P = p(i)   ;
        mom_L       = rho(i-1)*u(i-1) ;    	rho_L = rho(i-1) ;      u_L = u(i-1) ;      p_L = p(i-1) ;
        
        vel_flux_R  = rho_R*u_R*u_R +p_R ;    e_R = e(i+1) ;
        vel_flux_P  = rho_P*u_P*u_P +p_P ;    e_P = e(i)   ;
        vel_flux_L  = rho_L*u_L*u_L +p_L ;    e_L = e(i-1) ;
        
        energy_flux_R   = u_R * ( e_R + p_R ) ;
        energy_flux_P   = u_P * ( e_P + p_P ) ;
        energy_flux_L   = u_L * ( e_L + p_L ) ;
        
        rho_fluxR   = 0.5*( mom_P + mom_R ) -0.5*lamda*( rho_R - rho_P ) ;
        rho_fluxL   = 0.5*( mom_P + mom_L ) -0.5*lamda*( rho_P - rho_L ) ;
        
        vel_fluxR   = 0.5*( vel_flux_P + vel_flux_R ) -0.5*lamda*( mom_R - mom_P ) ;
        vel_fluxL   = 0.5*( vel_flux_P + vel_flux_L ) -0.5*lamda*( mom_P - mom_L ) ;
        
        energy_fluxR    = 0.5*( energy_flux_P + energy_flux_R ) -0.5*lamda*( e_R - e_P );
        energy_fluxL    = 0.5*( energy_flux_P + energy_flux_L ) -0.5*lamda*( e_P - e_L ) ;
        
        new_rho(i)  = rho_P - dt/dx * ( rho_fluxR - rho_fluxL ) ;
        vel_flux    = mom_P - dt/dx * ( vel_fluxR - vel_fluxL ) ;
        
        new_u(i)    = vel_flux/new_rho(i) ;
        new_e(i)    = e_P - dt/dx * ( energy_fluxR - energy_fluxL ) ;
        
    end
    
    rho     = new_rho ;
    u       = new_u ;
    e       = new_e ;
    
end

%pressure = dlmread('pressure.dat') ;
%density  = dlmread('density.dat')  ;
%velocity = dlmread('velocity.dat') ;



%Gudunovs stuff
%Gadunov scheme
L = 1;
dx_gud = 0.001;
x_gud = (0:dx_gud:L);
N_gud = (L/dx_gud)+1;
T_gud = 0.2;
dt_gud = 0.0001;
n_gud = (T_gud/dt_gud);
C = dt_gud/dx_gud;
gamma = 1.4;

%Initial conditions
for i = 1:N_gud
    if x_gud(i) <= 0.5
        rho_gud(i) = 1;
        u_gud(i) = 0;
        p_gud(i) = 2.5*(gamma - 1);
    else
        rho_gud(i) = 0.125;
        u_gud(i) = 0;
        p_gud(i) = 0.25*(gamma - 1);
    end
end
ET(:) = ((1/2).*rho_gud(:).*(u_gud(:).^2)) + (p_gud(:)/(gamma - 1));



for i = 1:n_gud

    %1) Construct Q and E vectors

    Q_gud(1,:) = rho_gud(:);
    Q_gud(2,:) = rho_gud(:).*u_gud(:);
    Q_gud(3,:) = ET(:);

    E_gud(1,:) = rho_gud(:).*u_gud(:);
    E_gud(2,:) = E_gud(1,:).*u_gud(1,:) + p_gud(1,:);
    E_gud(3,:) = ET(:).*u_gud(:) + p_gud(:).*u_gud(:);

    %2) Calculate Q(n+1) using Gadunov method
    Q1 = Q_gud;
    flux = zeros(3,N_gud);
    for k = 1:N_gud-1
        c = sqrt(gamma*p_gud(k)/rho_gud(k));
        A = [abs(Q_gud(2,k)) abs(Q_gud(2,k)+c) abs(Q_gud(2,k)-c)];
        al = max(A);
        flux(:,k) = 0.5*(E_gud(:,k) + E_gud(:,k+1)) - 0.5*(al)*(Q1(:,k+1) - Q1(:,k));
    end
    for k = 2:N_gud-1
        Q_gud(:,k) = Q1(:,k) - C*(flux(:,k) - flux(:,k-1));
    end
    
    %3) Update solution (rho,u,p)
    rho_gud(:) = Q_gud(1,:);
    u_gud(1,:) = Q_gud(2,:)./rho_gud(1,:);
    ET(:) = Q_gud(3,:);
    p_gud(1,:) = (ET(1,:) - ((1/2).*rho_gud(1,:).*(u_gud(1,:).^2)))*(gamma - 1);
    
end






time = 0.2;
data = Test_1(time);


figure(1)
hold on
%plot(density(:,1),density(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(data.x,data.rho,'b','LineWidth',2);
plot(xc,rho,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
plot(x_gud, rho_gud, 'g', 'LineWidth', 2);
xlabel(' x ','FontSize',18,'FontWeight','bold');
ylabel(' Density ','FontSize',18,'FontWeight','bold');
legend('Exact','Lax Friedrich','Gudunov','Location','northeast','FontSize',16,'FontWeight','bold');
%print(gcf,'Density.jpg','-dpng','-r300');

figure(2)
hold on
%plot(pressure(:,1),pressure(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(data.x,data.P,'b','LineWidth',2);
plot(xc,p,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
plot(x_gud, p_gud, 'g', 'LineWidth', 2);
xlabel(' x ','FontSize',18,'FontWeight','bold');
ylabel(' Pressure ','FontSize',18,'FontWeight','bold');
legend('Exact','Lax Friedrich','Gudunov','Location','northeast','FontSize',16,'FontWeight','bold');
%print(gcf,'Pressure.jpg','-dpng','-r300');

figure(3)
hold on
%plot(velocity(:,1),velocity(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(data.x,data.u,'b','LineWidth',2);
plot(xc,u,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
plot(x_gud, u_gud, 'g', 'LineWidth', 2);
xlabel(' x ','FontSize',18,'FontWeight','bold');
ylabel(' Velocity ','FontSize',18,'FontWeight','bold');
legend('Exact','Lax Friedrich','Gudunov','Location','south','FontSize',16,'FontWeight','bold');
%print(gcf,'Velocity.jpg','-dpng','-r300');





