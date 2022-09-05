clear all;

simulate = true; 

t0 = 0;       % Initial time
dt = 1e-4;    % Fixed time step
tf = 500000;       % Final time INCREASE IT IF NEEDED

tlen = tf/dt;
print_counter = tlen/20000; % for the histogram
%%% for D=1e-3; dt=1e-4, tf = 500000 (output sample 20000)
%%% for D=1e-4; dt=1e-2, tf=10000000
%%% for D=1e-5; dt=1e-2, tf=10000000

%%%% THIS IS WHAT WE SOLVE
%%%% dX_t = \beta -k*X_t - sqrt(D)*X_t*dW
%%%% \beta = g/(1+rho*R)
%%%% model parameters %%%%%

g=0.1;
k=5e-4;
rho=10;
R=100;
D=1e-5; % noise intensity

beta = g/(1+rho*R);
beta/(k+D)  %%% this is theoretical prediction of the stationary mode

f_func = @(t,x) beta - k*x(1);      % Drift function
g_func = @(t,x)-sqrt(D)*x;          % Diffusion function

if simulate==1
    x = zeros(tlen/print_counter,1); % Allocate output
    x_new = 0.0;
    x_prev = 0.0;
    x0 = beta/(k+D);   % set initial conditin to the expected stationary mode
    
    seed = 10;  % Seed value
    rng(seed); % Always seed random number generator
    
    % Euler-Heun for Stratonovich interpretation
    counter=1;
    for i = 1:tlen-1
        dW = randn(1,1)*sqrt(2*dt);
        tspan_i = t0+i*dt;
        gnbar = g_func(tspan_i,x_prev+g_func(tspan_i,x_prev)*dW);
        x_new = x_prev + f_func(tspan_i,x_prev)*dt + 0.5*(g_func(tspan_i,x_prev)+gnbar).*dW;
        x_prev = x_new;
        if(mod(i-1,print_counter)==0)
            x(counter,:) = x_new;
            counter = counter+1;
        end
    end
else
    %%% load from file
    x = importdata('C++ Solver/Outputs/output_D1e-5.txt');
end

x_mesh = linspace(0.001,0.5,1000);
%%% theoretical stationary distribution from the paper:
p = (beta./D).^(k./D).*x_mesh.^(-1-k./D).*exp(-beta./D./x_mesh)./gamma(k./D);

%%%%% PLOT

figure(2)
clf;
%%% for D=1e-3 'NumBins' 500000
%%% for D=1e-5 or 1e-4 'NumBins' around 100-300
h=histogram(x(:),'Normalization','pdf','NumBins',200);
xlim([0 0.5])
hold on
plot(x_mesh,p,'Color','green','LineWidth',3)

%%% mode from numerical results (should coincide with theoretical
%%% prediction)
xmode = 0.5*(h.BinEdges(find(h.Values==max(h.Values)))+h.BinEdges(find(h.Values==max(h.Values))+1));
xline(xmode,'color','red','linewidth',3)
xmode

%%% if need to save data %%%%
%save('x_milstein_long.mat','x')
% writetable(table(x),'tmp.txt','WriteRowNames',false)
