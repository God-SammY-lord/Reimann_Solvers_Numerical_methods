function [data] = Test_5(t)

if nargin < 1
    %set default value
    t = 0.25;
end
%Initial conditions
x0 = 0;
rho_l = 5.99924;
P_l = 460.894;
u_l = 19.5975;


rho_r =  5.99242;
P_r = 46.0950;
u_r = -6.19633;

gamma = 1.4;
mu = sqrt( (gamma-1)/(gamma+1) );

%speed of sound
a_l= power( (gamma*P_l/rho_l),0.5);
a_r = power( (gamma*P_r/rho_r),0.5);

A_R = 2/((gamma+1)*rho_r);
B_R = mu*mu*P_r;
A_L = 2/((gamma+1)*rho_l);
B_L = mu*mu*P_l;


%Test 1
P_star = fzero('Function_set',2000);
v_star = u_l - (P_star-P_l)*power((A_L/(P_star+B_L)),0.5);
v_shock_r = u_r + a_r*power((((gamma+1)/(2*gamma))*(P_star/P_r)) + ((gamma-1)/(2*gamma)) ,0.5);
v_shock_l = u_l - a_l*power((((gamma+1)/(2*gamma))*(P_star/P_l)) + ((gamma-1)/(2*gamma)) ,0.5);
rho_star_R = rho_r*(( (P_star/P_r) + mu^2 )/(1 + mu*mu*(P_star/P_r)));
rho_star_L = rho_l*(( (P_star/P_l) + mu^2 )/(1 + mu*mu*(P_star/P_l)));
%v_shock_r = v_star*((rho_star_R/rho_r)/( (rho_star_R/rho_r) - 1));
%v_shock_l = v_star*((rho_star_L/rho_l)/( (rho_star_L/rho_l) - 1));

%Key Positions
x1 = x0 - v_shock_l*t;
x2 = x0 + v_star*t;
x3 = x0 + v_shock_r*t;
%determining x2
%c_2 = a_l - ((gamma - 1)/2)*v_star;
%x2 = x0 + (v_star - c_2)*t;

disp("p*")
disp(P_star);
disp("u*")
disp(v_star);
disp("rho*L")
disp(rho_star_L);
disp("rho*R")
disp(rho_star_R);
%disp(v_shock);



%start setting values
n_points = 10000;    %set by user
%boundaries (can be set)
x_min = -0.5;
x_max = 0.5;

x = linspace(x_min,x_max,n_points);
data.x = x';
data.rho = zeros(n_points,1);   %density
data.P = zeros(n_points,1); %pressure
data.u = zeros(n_points,1); %velocity
data.e = zeros(n_points,1); %internal energy

for index = 1:n_points
    if data.x(index) < x1
        %Solution b4 x1 (b4 left shock)
        data.rho(index) = rho_l;
        data.P(index) = P_l;
        data.u(index) = u_l;
    elseif (x1 <= data.x(index) && data.x(index) <= x2)
        %Solution b/w x1 and x2 (btw left shock and contact)
        data.rho(index) = rho_star_L;
        data.P(index) = P_star;
        data.u(index) = v_star;
    elseif (x2 <= data.x(index) && data.x(index) <= x3)
        %Solution b/w x2 and x3 (btw contact and right shock)
        data.rho(index) = rho_star_R;
        data.P(index) = P_star;
        data.u(index) = v_star;
    elseif x3 < data.x(index)
        %Solution after x3 (after right shock)
        data.rho(index) = rho_r;
        data.P(index) = P_r;
        data.u(index) = u_r;
    end
    data.e(index) = data.P(index)/((gamma - 1)*data.rho(index));
end
end
