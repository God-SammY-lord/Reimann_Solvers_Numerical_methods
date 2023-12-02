function [data] = Test_4(t)

if nargin < 1
    %set default value
    t = 0.25;
end
%Initial conditions
x0 = 0;
rho_l = 1;
P_l = 0.01;
u_l = 0;


rho_r = 1;
P_r = 100;
u_r = 0;

gamma = 1.4;
mu = sqrt( (gamma-1)/(gamma+1) );

%speed of sound
a_l= power( (gamma*P_l/rho_l),0.5);
a_r = power( (gamma*P_r/rho_r),0.5);

%Test 1
P_star = fzero('Function_set',50);
%v_star = 2*(sqrt(gamma)/(gamma - 1))*(1 - power(P_star, (gamma - 1)/(2*gamma)));
v_star = u_r - 2*(a_r/(gamma-1))*(1-power((P_star/P_r),((gamma-1)/(2*gamma))));
rho_star_R = (rho_r)*power((P_star/P_r),1/gamma);
rho_star_L = rho_l*(( (P_star/P_l) + mu^2 )/(1 + mu*mu*(P_star/P_l)));
v_shock = v_star*((rho_star_L/rho_l)/( (rho_star_L/rho_l) - 1));

%Key Positions
x1 = x0 + v_shock*t;
x2 = x0 + v_star*t;
x4 = x0 + a_r*t;
%determining x3
c_3 = a_r + ((gamma - 1)/2)*v_star;
x3 = x0 + (v_star + c_3)*t;

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
n_points = 1000000;    %set by user
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
        %Solution b4 x1 (b4 shock)
        data.rho(index) = rho_l;
        data.P(index) = P_l;
        data.u(index) = u_l;
    elseif (x1 <= data.x(index) && data.x(index) <= x2)
        %Solution b/w x1 and x2 (shock to contact)
        data.rho(index) = rho_star_L;
        data.P(index) = P_star;
        data.u(index) = v_star;
    elseif (x2 <= data.x(index) && data.x(index) <= x3)
        %Solution b/w x2 and x3 (contact to fan)
        data.rho(index) = rho_star_R;
        data.P(index) = P_star;
        data.u(index) = v_star;
    elseif (x3 <= data.x(index) && data.x(index) <= x4)
        %Solution b/w x3 and x4 (fan)
        c = -mu*mu*(u_r- ((data.x(index))/t)) + (1 - mu*mu)*a_r; 
        data.rho(index) = rho_r*power((c/a_r),2/(gamma - 1));
        data.P(index) = P_r*power((data.rho(index)/rho_r),gamma);
        data.u(index) = (1 - mu*mu)*( ((data.x(index))/t) - a_r+((gamma-1)/2)*u_r);
    elseif x4 < data.x(index)
        %Solution after x4 (after fan)
        data.rho(index) = rho_r;
        data.P(index) = P_r;
        data.u(index) = u_r;
    end
    data.e(index) = data.P(index)/((gamma - 1)*data.rho(index));
end
end
