Ns = length(N)*length(controllers);
%%
% f = figure(5);
% clf
% 
% % leg = string(1:Ns);
% for x = 1:Ns
%     if mod(x,2)==1 % Controller C1
% %         leg(x) = strcat("C1 N=",);
%         colvec = [1 0 0];
%     else % Controller C2
% %         leg(x) = strcat("C2 N=",string(mod(x,length(controllers))));
%         colvec = [0 1 0];
%     end
%     if invChains(x).ctrl.control == "C1"
%         plot(invChains(x).time,invChains(x).norms.l2,...
%             '-','color',colvec.*x/Ns)    
%     elseif invChains(x).ctrl.control == "C2"
%         plot(invChains(x).time,invChains(x).norms.l2,...
%             '--','color',colvec.*x/Ns)    
%     end
%     leg(x) = strcat("N=",num2str(invChains(x).sys.n)," ",...
%         invChains(x).ctrl.control);
%     hold on
% end
% 
% xlabel('Time [s]')
% ylabel('$\sup_i\left|x_i-x_i^*\right|_2$')
% legend(leg)

%%

f = figure(6); clf; 
f.Name = 'linf_x';
set(groot,'DefaultAxesUnits','normalized')

if length(controllers)>1
    norms_ours_v = [invChains(1:2:end).norms];
    norms_naiv_v = [invChains(2:2:end).norms];
    stem(N(1:end)-.1,[norms_ours_v(1:end).linf],'-*','r','linewidth',2); hold on
    stem(N(1:end)+.1,[norms_naiv_v(1:end).linf],'--','g','linewidth',2);
    pause(.5)
    legend('DSS Controller $C_1$','Standard Controller $C_2$','location','best')
else    
    norms_ours_v = [invChains(1:end).norms];
%     norms_naiv_v = [invChains(2:2:end).norms];
    stem(N-.1,[norms_ours_v.linf],'--*','r','linewidth',2); hold on
%     stem(N+.1,[norms_naiv_v.linf],'--*','g','linewidth',2)    
end
xlabel('Number of agents $n$')
ylabel('$\sup_i\left|\left|x_i-x_i^*\right|\right|_2$')
xlim([N(1)-1 N(end)+1])
xticks(N)
f.Units='centimeters';
f.Position(4) = 7;

%% Max P Power
fx = figure(7); clf; 

set(groot,'DefaultAxesUnits','normalized')
maxP = zeros(length(N),1);
for idxN = 1:Ns    
    if plots.powerinj
        fx.Name = 'linf_Pinj';
        maxP(idxN) = max(max(abs(invChains(idxN).Pi)));
    else
        fx.Name = 'linf_P';
        maxP(idxN) = max(max(abs(invChains(idxN).P)));        
    end
%     maxP(idxN) = min(min((invChains(idxN).P)));
end
if length(controllers)>1
    stem(N(1:end)-.1,maxP(1:2:end),'-*','r','linewidth',2); hold on
    stem(N(1:end)+.1,maxP(2:2:end),'--','g','linewidth',2);
    legend('DSS Controller $C_1$','Standard Controller $C_2$','location','best')
else
    stem(N-.1,maxP(1:end),'--*','r','linewidth',2); hold on
end
xlabel('Number of agents $n$')
if plots.powerinj
    ylabel('$\sup_i\left|\left|P_{I,i}\right|\right|_2$')
else
    ylabel('$\sup_i\left|\left|P_i\right|\right|_2$')
end
xlim([N(1)-1 N(end)+1])
ylim([min(maxP)-10 max(maxP)+10])
xticks(N)
fx.Units='centimeters';
fx.Position(4) = 7;

%% Max Q Power
fx = figure(8); clf; 

set(groot,'DefaultAxesUnits','normalized')
maxQ = zeros(length(N),1);
for idxN = 1:Ns    
    if plots.powerinj
        fx.Name = 'linf_Qinj';
        maxQ(idxN) = max(max(abs(invChains(idxN).Qi)));
    else
        fx.Name = 'linf_Q';
        maxQ(idxN) = max(max(abs(invChains(idxN).Q)));        
    end
end
if length(controllers)>1
    stem(N(1:end)-.1,maxQ(1:2:end),'-*','r','linewidth',2); hold on
    stem(N(1:end)+.1,maxQ(2:2:end),'--','g','linewidth',2);
    legend('DSS Controller $C_1$','Standard Controller $C_2$','location','best')
else
    stem(N-.1,maxQ(1:end),'--*','r','linewidth',2); hold on
end
xlabel('Number of agents $n$')
if plots.powerinj
    ylabel('$\sup_i\left|\left|Q_i\right|\right|_2$')
else
    ylabel('$\sup_i\left|\left|Q_i\right|\right|_2$')
end
xlim([N(1)-1 N(end)+1])
ylim([min(maxQ)-10 max(maxQ)+10])
xticks(N)
fx.Units='centimeters';
fx.Position(4) = 7;