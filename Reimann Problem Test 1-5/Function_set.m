function y = Function_set(P)
    %Initial conditions
rho_l = 5.99924;
P_l = 460.894;
u_l = 19.5975;


rho_r =  5.99242;
P_r = 46.0950;
u_r = -6.19633;

    gamma = 1.4;
    mu = sqrt( (gamma-1)/(gamma+1) );
    A_R = 2/((gamma+1)*rho_r);
    B_R = (mu*mu*P_r);
    A_L = 2/((gamma+1)*rho_l);
    B_L = (mu*mu*P_l);
    a_l= power( (gamma*P_l/rho_l),0.5);
    a_r = power( (gamma*P_r/rho_r),0.5);


       
    % Test 1 Sod Problem
    %y = (P - P_r)*((1-mu*mu)*(rho_r*(P + mu*mu*P_r))^-1)^(0.5)-(power(P_l , (gamma-1)/(2*gamma))-power(P , (gamma-1)/(2*gamma)))*(((1-mu*mu*mu*mu)*P_l^(1/gamma))*(mu*mu*mu*mu*rho_l)^-1)^(0.5);

    %Test 2 123 Problem
    %y = ((2*a_l)/(gamma-1))*(-1+power((P/P_l),((gamma-1)/(2*gamma)))) + ((2*a_r)/(gamma-1))*(-1+power((P/P_r),((gamma-1)/(2*gamma)))) +u_r-u_l;

    %Test 3 
    %y = (P - P_r)*((1-mu*mu)*(rho_r*(P + mu*mu*P_r))^-1)^(0.5)-(power(P_l , (gamma-1)/(2*gamma))-power(P , (gamma-1)/(2*gamma)))*(((1-mu*mu*mu*mu)*P_l^(1/gamma))*(mu*mu*mu*mu*rho_l)^-1)^(0.5);

    %Test 4
    %y = (P-P_l)*power((A_L/(P+B_L)),0.5) + ((2*a_r)/(gamma-1))*(-1+power((P/P_r),((gamma-1)/(2*gamma)))) +u_r-u_l;

    %Test 5
    y = (P - P_r)*power(((A_R/(B_R+P))),0.5) + (P - P_l)*power(((A_L/(B_L+P))),0.5) + u_r-u_l;
end
