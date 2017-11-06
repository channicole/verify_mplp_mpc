clear;clc;
%%
%set up the initial conditions: total time; time step(optional); initial
%states; perturbation of the initial state;
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\SeDuMi_1_3'));
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\yalmip'));
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\SDPT3'))
%%
% load ../Trajectories/Twelve_Traj.mat
% Twelve_Traj = Twelve_Traj([1:2000],:);

ops = sdpsettings('verbose',0);
enable_disp = 0;
dt = 0.0001;
t_dur = 0:dt:10;
inner_step = 20;
initial_state = [1.505,1.505,0.005,0.005]';
dimension = 4;

%Initial_P = eye(dimension);
Initial_P = [2.5063   -0.0276    0.9436    0.0280 ; -0.0276    2.4722    0.0013    0.9537; 0.9436    0.0013    2.2441    0.0151; 0.0280    0.9537 0.0151    2.2845];
Initial_l = (0.005)^2;

delta_initial_up = zeros(4,1);
delta_initial_up(1) = sqrt(Initial_l);
%Lipschitz constant of te system, user input
Lip_constant = 0;
[t_precise,x_precise] = ode45(@robot,t_dur,initial_state);
% [t_true,x_true] = ode45(@Twelve,t_dur,initial_state+delta_initial_up);
%[t_true_low,x_true_low] = ode45(@Engine,t_dur,initial_state+delta_initial_low);
% t_precise = Twelve_Traj(:,1);
% x_precise = Twelve_Traj(:,[2:13]);

%% Get index of variables
% INDEX_CH = zeros(2^8,8);
% for i = 0 : 2^8-1
%     INDEX_CH(i+1,:) = de2bi(i,8);
% end
% Jacbian_Function_index = [3,2;3,3;4,3;4,4;5,4;5,5;7,2;7,7];

%%
if enable_disp
    figure;
    scatter(x_precise(:,1),x_precise(:,2));
    xlim = ylim;
    hold on;
    %scatter(x_true(:,1),x_true(:,2),'r');
    scatter(initial_state(1),initial_state(2),180,'x','b','linewidth',2);
    %scatter(initial_state(1)+delta_initial_up(1),initial_state(2)+delta_initial_up(2),180,'x','r','linewidth',2);
    xlabel('x1','FontSize', 24);
    ylabel('x2','FontSize', 24);
    if enable_disp
        plot_ellipsoid(x_precise(1,[1,2])',Initial_P([1,2],[1,2])/Initial_l,'k');
    end
end



index = 1:inner_step:length(t_precise);
t = t_precise(index);
x = x_precise(index,:);

% figure
% for i = 1:size(t,1)
%     lambda = max(real(eig(Twelve_Jac(x(i,:)))));
%     tmp(i) = lambda;
%     [gamma,Q] = QLF_restricted(Twelve_Jac(x(i,:)));
%     disp(gamma);
%     disp(max(eig(Q)));
% end
% plot(tmp)

TELASPE = tic;
%%
n = length(x);   
%dimention of the system

Reach_info = zeros(n,dimension^2+1);
Reach_info(1,:) = [reshape(Initial_P,1,dimension^2),Initial_l];

for i = 1:length(t)-1
%for i = 1:1
    if enable_disp
        disp('----------------');
        disp(['Computing ', num2str(t(i)), ' time interval ']);
    end
    M = reshape(Reach_info(i,1:dimension^2),dimension,dimension);
    level = Reach_info(i,dimension^2+1);

    Average_value = mean(x_precise(index(i):index(i+1)-1,:));
    A_c = robot_jac(Average_value);
    shrink = max(real(eig(A_c)));
    %shrink = -0.01;
    
    %Local Lipschitz constant
    Lip = norm(A_c);   
   %delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_max = (max(x(i+1,:),x(i,:))+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_min = (min(x(i+1,:),x(i,:))-sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    if max(eig(A_c'*M+M*A_c-shrink*M)) > 0 
        check_ch_flag = 0;
    else
        check_ch_flag = 1;
    end  
    
    if check_ch_flag == 1 %&& check_all_vertices==1
%         Q = M;
%         if check_all_vertices == 1
%             disp('Pass all vertex check for old shrink')      
%             gamma = shrink;
%         else
        gamma = shrink;
        Delta = abs(robot_jac(x_relative_min)-robot_jac(x_relative_max));
        SDP_Delta_matrix = Delta'*Q+Q*Delta;
        epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
        %disp(epsilon_matrix);
        if enable_disp
            disp('Estimation error!')
            disp(gamma*min(eig(Q)) + epsilon_matrix);
        end
        gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q))/2;
        shrink = gamma;
        %end
        level_next = level*exp(shrink*(t(i+1)-t(i)));
        if enable_disp
            disp('remain P in last round');
            disp('expontial change rate:');
            disp(shrink);
        end
        P_level_original = level;
    else
        if enable_disp
            disp('Compute new P')
        end
        %Lyapunov-like function
        [gamma,Q] = QLF_restricted(A_c);
%         check_rest_vertices = 1;
%         [Current_A_min,Current_A_max] = LL_Jac_Vertice(x_relative_min,x_relative_max);
%         Current_Vertex = Current_A_min;
%         for cnt_ver = 1:2^8
%             for cnt_inx = 1:8
%                 if INDEX_CH(cnt_ver,cnt_inx) == 1
%                     Current_Vertex(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2)) = Current_A_max(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2));
%                 end
%             end
%             while max(eig(Current_Vertex'*Q+Q*Current_Vertex - gamma*Q)) > 0
%                 if gamma <=0.3
%                     gamma = gamma + 0.1;
%                 else
%                     disp('fail using new computed Q and gamma for rest check');
%                     check_rest_vertices = 0;
%                     break;
%                 end
%             end
%             if check_rest_vertices == 0
%                 break
%             end
%         end
%        %disp(gamma);
%         if check_rest_vertices == 1
%             disp('pass the rest vertex check')
%             disp(gamma);
%         else
%             [gamma,Q] = QLF(A_c);
        delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
        Delta_matrix = abs(robot_jac(Average_value+delta_theta)-robot_jac(Average_value-delta_theta))/2;
        %disp((Delta_matrix));
        SDP_Delta_matrix = Delta_matrix'*Q+Q*Delta_matrix;
        epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
        if enable_disp
            disp('Esitimation Error!')
            disp(gamma*min(eig(Q)) + epsilon_matrix);
        end
        gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q))/2;
        if enable_disp
            disp('using interval matrix disturbance to over-approximate \gamma')
            disp('expontial change rate:');
            disp(gamma);
        end
%         else
%             if enable_disp
%                 disp('Computed P satifies convex hull constrains')
%                 disp('expontial change rate:');
%                 disp(gamma);
%             end
%    end
    %     disp('eigenvalue of P');disp(eig(Q));
        P_level_original = min_level_set4(Q,M,level);
    end
    level_next = P_level_original*(exp(gamma*(t(i+1)-t(i))));
    for cnt = index(i):index(i)+inner_step-1
        if enable_disp
            plot_ellipsoid(x_precise(cnt,[1,2])',Q([1,2],[1,2])/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))),'k');
        end
        fp_checker_inner = fix_point_checker(x_precise(cnt+1,:),x_precise(cnt,:),Q/(P_level_original*(exp(gamma*(t_precise(cnt+1)-t(i))))),Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))));
        if fp_checker_inner == 1
            if enable_disp
                disp('find invariant set in between');
            end
            break;
        end
    end
    if enable_disp
       plot_ellipsoid(x(i+1,[1,2])',Q([1,2],[1,2])/level_next,'k');
    end
    Reach_info(i+1,:) = [reshape(Q,1,dimension^2),level_next];
    fp_checker = fix_point_checker(x(i+1,:),x(i,:),Q/level_next,M/level);
    if fp_checker == 1 || fp_checker_inner == 1
        disp('find invariant set');
        break;
    end
end

toc(TELASPE);

Reach_final = reshape(Reach_info(i,1:dimension^2),dimension,dimension)/Reach_info(i,dimension^2+1);
dia_final = sqrt(eig(inv(Reach_final)));
disp(prod(dia_final)/(Initial_l)^6);

dim = 4;


Reach_initial = reshape(Reach_info(1,1:dim^2),dim,dim)/Reach_info(1,dim^2+1);
dia_initial = sqrt(eig(inv(Reach_initial)));
Reach_final = reshape(Reach_info(i,1:dim^2),dim,dim)/Reach_info(i,dim^2+1);
dia_final = sqrt(eig(inv(Reach_final)));
disp('F/I:');
disp(prod(dia_final)/prod(dia_initial));

sum_volume = 0;
for i = 1:length(t)
    Reach_ellip = reshape(Reach_info(i,1:dim^2),dim,dim)/Reach_info(i,dim^2+1);
    dia_now = sqrt(eig(inv(Reach_ellip)));
    current_volume = prod(dia_now);
    sum_volume = sum_volume + current_volume;
end
sum_volume = sum_volume/length(t);
disp('Average/Initial:')
disp(sum_volume/prod(dia_initial));
