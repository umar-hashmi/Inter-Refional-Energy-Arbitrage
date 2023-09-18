%% Multi-regional Energy Arbitrage with linear programming
% Code Author: Md Umar Hashmi
% Date: 7 August 2023
% Location: KU Leuven and EnergyVille, Belgium

clear
close all
clc

tic
load('BE_UK_interconFlow_2019.mat');

e_ch=0.95;                          % Charging efficiency
e_dis =0.95;                        % Discharging efficiency
e_con = 0.95;
e_net_ch = e_ch*e_con;
e_net_dis = e_dis*e_con;

% battery parameters
del_max = 0.5;                      % Maximum charging rate
del_min = -del_max;                 % Minimum discharging rate
b_0 = 0.5;                          % Initial battery capacity
b_max = 1;                          % Maximum battery capacity
b_min = 0.1;                        % Minimum permissible battery capacity
h=1;

n_days = length(UKBE19)/24;

profit_mat =zeros(n_days,3);
b_cap_mat = zeros(24,n_days);
x1_mat = zeros(24,n_days);
x2_mat = zeros(24,n_days);
x3_mat = zeros(48,n_days);

del_min_adj_mat = zeros(24,n_days);
del_max_adj_mat = zeros(24,n_days);

int_rent = 2;                       %euros per MWh

int_efficiency = 1-0.025;

samp = 24;                          % total samples

epsilon = 1e-3;

for i = 1:n_days
    st = (i-1)*samp+1;
    ed = samp*i;
    pri_be = max(UKBE19(st:ed,1),epsilon);
    pri_uk = max(UKBE19(st:ed,2),epsilon);
    int_flow = UKBE19(st:ed,3);
    
    P_B_A_mat = pri_be/e_net_ch;
    P_S_A_mat = pri_be*e_net_dis;
    
    P_B_B_mat = (pri_uk + int_rent)/(int_efficiency*e_net_ch);
    P_S_B_mat = (pri_uk - int_rent)*(int_efficiency*e_net_dis);
    
    C1 = diag(P_B_A_mat);
    C2 = diag(P_S_A_mat);
    C3 = diag(P_B_B_mat);
    C4 = diag(P_S_B_mat);
    Ze = zeros(samp);
    Oe = ones(samp);
    Tr = tril(ones(samp));
    Ey = eye(samp);
    
    del_min_adj = zeros(samp,1);
    del_max_adj = zeros(samp,1);
    
    for kj =1:samp
        if int_flow(kj) >= 0
            del_min_adj(kj,1) = min(0, max(del_min, -1000 + int_flow(kj)));
            del_max_adj(kj,1) = del_max;
        else
            del_min_adj(kj,1) = del_min;
            del_max_adj(kj,1) = max(0, min(del_max, 1000 + int_flow(kj)));
        end
    end
        
    A = [C1,    Ze,     -1*Ey,      Ze,     Ze,     Ze;...
        C2,     Ze,     -1*Ey,      Ze,     Ze,     Ze;...
        Ze,     C3,     Ze,         -1*Ey,  Ze,     Ze;...
        Ze,     C4,     Ze,         -1*Ey,  Ze,     Ze;...
        Tr,     Tr,     Ze,         Ze,     Ze,     Ze;...
        -1*Tr,  -1*Tr,  Ze,         Ze,     Ze,     Ze;...
        Ey,     Ey,     Ze,         Ze,     Ze,     Ze;...
        -1*Ey, -1*Ey,   Ze,         Ze,     Ze,     Ze;...
        -1*Ey,  Ze,     Ze,         Ze,     del_min*Ey, Ze;...
        Ey,     Ze,     Ze,         Ze,     Ze,     -1*del_max*Ey;...
        Ze,     -1*Ey,  Ze,         Ze,     diag(del_min_adj), Ze;...
        Ze,     Ey,     Ze,         Ze,     Ze,     -1*diag(del_max_adj);...
        Ze,     Ze,     Ze,         Ze,     Ey,     Ey];
    
    b = [zeros(4*samp,1);...
        (b_max-b_0)*ones(samp,1);...
        (b_0-b_min)*ones(samp,1);...
        del_max*ones(samp,1);...
        -del_min*ones(samp,1);...
        zeros(4*samp,1);...
        ones(samp,1)];
    
    lb=[del_min*h*ones(samp,1);del_min*h*ones(samp,1); -1000000*ones(samp,1);-1000000*ones(samp,1);  zeros(2*samp,1)];
    ub=[del_max*h*ones(samp,1);del_max*h*ones(samp,1); 1000000*ones(samp,1);1000000*ones(samp,1);  ones(2*samp,1)];
    
    Aeq=[];
    beq=[];
    f=[zeros(samp,1);zeros(samp,1); ones(samp,1);ones(samp,1); zeros(2*samp,1)];
    intcon = (4*samp+1):6*samp;
    
    x_state = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
    
    X1 = x_state(1:samp);
    X2 = x_state(samp+1:2*samp);
    X3 = x_state(2*samp+1:4*samp);
    Profit1 = sum(P_B_A_mat'*subplus(X1) - P_S_A_mat'*subplus(-X1));
    Profit2 = sum(P_B_B_mat'*subplus(X2) - P_S_B_mat'*subplus(-X2));
    
    profit_mat(i,:) = [Profit1, Profit2, Profit1+Profit2];
    
    b_cap = [b_0 + cumsum(X1) + cumsum(X2)];
    
    b_cap_mat(:,i) = b_cap;
    x1_mat(:,i) = X1;
    x2_mat(:,i) = X2;
    x3_mat(:,i) = X3;
    
    del_min_adj_mat(:,i) = del_min_adj;
    del_max_adj_mat(:,i) = del_max_adj;
    
    b_0 = b_cap(end);
end

toc

figure
plot(-1*profit_mat)

tot=-1*sum(profit_mat);

[total_cycle_100] = calculate_cycles(x1_mat+x2_mat);

euros_per_cycle = -1*profit_mat(:,3)./total_cycle_100;

%%
% id =304;
id = 23;
figure
plot(x1_mat(:,id))
hold on
plot(x2_mat(:,id))
hold on
plot(x2_mat(:,id)+x1_mat(:,id), '+')
hold on
plot(del_min_adj_mat(:,id), '*')
hold on
plot(del_max_adj_mat(:,id), 'o')

figure
plot(-1*profit_mat)
