function plotResults(invChains,ChainToPlot,N,control,vars)
time = invChains(ChainToPlot).time; 
y = invChains(ChainToPlot).y;
u = invChains(ChainToPlot).u;
n =  invChains(ChainToPlot).n;
S = diag(ones(1,n))-diag((ones(1,n-1)),-1); % to calculate ang differences
    %S=S(1:n-1,1:n-1);
Del = S'*y(:,1:n)'; % delta (angle differences)
    
%% Angles plot
figure(1)
if strcmp(vars,'incremental')    
    plot(time,[zeros(length(time),1) wrapToPi(y(:,1:n-1))],'LineWidth',2);
%     plot(time,[zeros(length(time),1) y(:,1:n-1)],'LineWidth',2);
else %% for absolute
    plot(time,wrapToPi(y(:,1:n)),'LineWidth',2);
%     plot(time,(y(:,1:n)),'LineWidth',2)
end
title('Angles \delta_i')
set(gca,'FontSize',16,'LineWidth',3,'Box','on')
grid on

%% Angle differences plot
figure(2)
plot(time,wrapToPi(Del(1:n-1,:))','LineWidth',2)
% plot(time,Del(1:n-1,:)','LineWidth',2)
title('Angle differences \delta_{ij}')
set(gca,'FontSize',16,'LineWidth',3,'Box','on')
grid on

%% Frequency plot
figure(3)
if strcmp(vars,'incremental')
    plot(time,y(:,n:2*n-1),'LineWidth',2) 
else %% for absolute
    plot(time,y(:,n+1:2*n),'LineWidth',2) 
end
title('Frequency')
set(gca,'FontSize',16,'LineWidth',3,'Box','on')
grid on

%% Secondary Variable plot
figure(4)
if strcmp(vars,'incremental')
    plot(time,y(:,2*n:end),'LineWidth',2) 
else %% for absolute
    plot(time,y(:,2*n+1:end),'LineWidth',2) 
end
title('Secondary variable')
set(gca,'FontSize',16,'LineWidth',3,'Box','on')
grid on

%% "Control" effort
% if control=="DSS"
    figure(5)
    if strcmp(vars,'incremental')
        plot(time,u,'LineWidth',2);
    else %% absolute

    end
    title('Control effort')
    set(gca,'FontSize',16,'LineWidth',3,'Box','on')
    grid on
% end
%% Infty norms plot
figure(6)
if strcmp(vars,'incremental')
    stem(N,[invChains([invChains.control]==control).z_infnorm],'LineWidth',2) 
else %% for absolute 
    
end
axis([N(1)-10 N(end)+10 min([invChains([invChains.control]==control).z_infnorm])*0.95 max([invChains([invChains.control]==control).z_infnorm])*1.05])
xticks(N)
title('Inf norm for different strings')
set(gca,'FontSize',14,'LineWidth',1.5,'Box','on')
grid on
end