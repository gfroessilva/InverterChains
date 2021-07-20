function plotResults(invChains,ChainToPlot,N,control,vars,wrap)
time = invChains(ChainToPlot).time; 
y = invChains(ChainToPlot).y;
u = invChains(ChainToPlot).u;
P = invChains(ChainToPlot).P;
Q = invChains(ChainToPlot).Q;
n =  invChains(ChainToPlot).sys.n;
S = diag(ones(1,n))-diag((ones(1,n-1)),-1); % to calculate ang differences

if strcmp(vars,'incremental') 
    S=S(1:n-1,1:n-1);
    Del = S'*y(:,1:n-1)'; % delta (angle differences)    
else
    Del = S'*y(:,1:n)'; % delta (angle differences)    
end


    
%% Angles plot
figure(1)
if strcmp(vars,'incremental')        
    if wrap
        plot(time,[wrapToPi(y(:,1:n-1)) zeros(length(time),1)],'LineWidth',2);
    else
        plot(time,[y(:,1:n-1) zeros(length(time),1)],'LineWidth',2);
    end
    title('Angles $\Theta_i$','interpreter','latex')
    legend(string(1:n),'location','best')
else %% for absolute
    if wrap
        plot(time,wrapToPi(y(:,1:n)),'LineWidth',2);
    else
        plot(time,(y(:,1:n)),'LineWidth',2);
    end    
    title('Angles $\delta_i$','interpreter','latex')
    legend(string(1:n),'location','best')
%     plot(time,(y(:,1:n)),'LineWidth',2)
end
set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% Angle differences plot
figure(2)
if wrap
    plot(time,wrapToPi(Del(1:n-1,:))','LineWidth',2)
else
    plot(time,Del(1:n-1,:)','LineWidth',2)
end
if strcmp(vars,'incremental')      
    title('Angle differences $\Theta_{ij}$','interpreter','latex')
else
    title('Angle differences $\delta_{ij}$','interpreter','latex')
end
legend(strcat(string(1:n-1),string((1:n-1)+1)),'location','best')
set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% Frequency plot
figure(3)
if strcmp(vars,'incremental')
    plot(time,y(:,n:2*n-1),'LineWidth',2) 
    title('Frequency $\tilde \omega_i$','interpreter','latex')
else %% for absolute
    plot(time,y(:,n+1:2*n),'LineWidth',2) 
    title('Frequency $\omega_i$','interpreter','latex')
end
legend(string(1:n),'location','best')

set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% Voltage plot
figure(4)
if strcmp(vars,'incremental')
    plot(time,y(:,2*n:3*n-1),'LineWidth',2) 
    title('Voltage $\tilde V_i$','interpreter','latex')
else %% for absolute
    plot(time,y(:,n+1:2*n),'LineWidth',2) 
    title('Voltage $\tilde V_i$','interpreter','latex')
end
legend(string(1:n),'location','best')

set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% Secondary Variable plot (freq)
figure(5)
if strcmp(vars,'incremental')
    plot(time,y(:,3*n:4*n-1),'LineWidth',2)
    title('Secondary variable $\tilde{\xi_i}$','interpreter','latex')
else %% for absolute    
    plot(time,y(:,2*n+1:end),'LineWidth',2)
    title('Secondary variable $\xi_i$','interpreter','latex')
end
legend(string(1:n),'location','best')
set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% Secondary Variable plot (volt)
figure(6)
if strcmp(vars,'incremental')
    plot(time,y(:,4*n:end),'LineWidth',2)
    title('Secondary variable $\tilde{\zeta_i}$','interpreter','latex')
else %% for absolute    
    plot(time,y(:,2*n+1:end),'LineWidth',2)
    title('Secondary variable $\zeta_i$','interpreter','latex')
end
legend(string(1:n),'location','best')
set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% "Control" effort
figure(7)
if strcmp(vars,'incremental')
    plot(time,u,'LineWidth',2);   
    title('Control effort','interpreter','latex')
else %% absolute
    plot(time,u,'LineWidth',2);   
    title('Control effort','interpreter','latex')
end
legend(string(1:n),'location','best')
set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% Active Power P
figure(8)
if strcmp(vars,'incremental')
    plot(time,P,'LineWidth',2);   
    title('Active Power P','interpreter','latex')
else %% absolute
    plot(time,P,'LineWidth',2);   
    title('Active Power P','interpreter','latex')
end
legend(string(1:n),'location','best')
set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% Reactive Power Q
figure(9)
if strcmp(vars,'incremental')
    plot(time,Q,'LineWidth',2);   
    title('Rective Power Q','interpreter','latex')
else %% absolute
    plot(time,Q,'LineWidth',2);   
    title('Rective Power Q','interpreter','latex')
end
legend(string(1:n),'location','best')
set(gca,'FontSize',12,'LineWidth',3,'Box','on')
grid on

%% Infty norms plot
if length(N)<2
    return
else
    figure(10)
end
if strcmp(vars,'incremental')
    stem(N,[invChains([invChains.control]==control).z_infnorm],'LineWidth',2) 
else %% for absolute 
    stem(N,[invChains([invChains.control]==control).z_infnorm],'LineWidth',2) 
end
axis([N(1)-10 N(end)+10 min([invChains([invChains.control]==control).z_infnorm])*0.95 max([invChains([invChains.control]==control).z_infnorm])*1.05])
xticks(N)
title('$\|\cdot\|_\infty$ for different strings','interpreter','latex')
set(gca,'FontSize',14,'LineWidth',1.5,'Box','on')
grid on
end