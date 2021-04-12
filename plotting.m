fprintf('Plotting...\n')

Nplot = 12;
% Nplot = 28;
Cplot = controllers(1);
idx_N = find(ismember(N,Nplot));

ylim1 = [-Inf Inf];
% xlim = [0 1000];
ylimP = [-Inf Inf];

subplots = 1;

callstack = dbstack('-completenames');
if strcmp(callstack(end).name,'plotting')
    close all
    plots.powerinj = 1;
    plots.delta = 1;
    plots.omega = 1;
    plots.xi = 0;
    plots.P = 1;
    plots.Q = 0;
    plots.Pfactor = 0;
end

if ~ismember(Nplot,N)
    idx_N = find(ismember(N(1),N));
    Nplot = N(idx_N);
end

if Cplot=="C2" && length(controllers)>1
    idx_N = idx_N+1 + 1*(idx_N-1);
elseif length(controllers)>1
    idx_N = idx_N + 1*(idx_N-1);
end

y = invChains(idx_N).y;
n = invChains(idx_N).sys.n;
time = invChains(idx_N).time;
P_final = invChains(idx_N).P;
Q_final = invChains(idx_N).Q;
Pi_final = invChains(idx_N).Pi;
Qi_final = invChains(idx_N).Qi;
    
colvect = invChains(idx_N).colvec;

set(groot,'DefaultAxesFontSize', 12)
set(groot,'DefaultAxesLineWidth',1.5)
set(groot,'DefaultAxesBox','on')
set(groot,'DefaultLineLinewidth',2.5)
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultAxesColorOrder',colvect)
set(groot,'DefaultAxesColorMap',colvect)
set(groot,'DefaultAxesYGrid','on')
set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesUnits','centimeters')

delta = (y(:,1:n)); % WRAP TO PI
omegastart = (conf.omegaast*time);
omega = y(:,n+1:2*n);
V = y(:,2*n+1:3*n);
xi = y(:,3*n+1:4*n);
zeta = y(:,4*n+1:5*n);

[~, idx_DistOn] = min(abs(time - flag.DistOnT));
[~, idx_DistOff] = min(abs(time - flag.DistOffT));

xlimm = [flag.DistOnT flag.DistOnT + 20];
xlimm2 = [flag.DistOffT-5 flag.DistOffT + 20];

[~, idx_xlim1] = min(abs(time - xlimm(2)));
[~, idx_xlim21] = min(abs(time - xlimm2(1)));
[~, idx_xlim22] = min(abs(time - xlimm2(2)));

xlim1 = 1:idx_xlim1+1;
xlim2 = idx_xlim21-1:idx_xlim22+1;

delta_s = zeros(size(delta));


if flag.DistOffT == sys.simtime
    delta_s(:) = repmat(delta(idx_DistOn-1,:) - omegastart(idx_DistOn-1),length(delta),1);
else
    delta_s(:) = repmat(delta(end-1,:) - omegastart(end-1),length(delta),1);
end

if flag.Dist==1
    delta_s(idx_DistOn:idx_DistOff,:) = repmat(delta(idx_DistOff-1,:) ...
        - omegastart(idx_DistOff-1),idx_DistOff-idx_DistOn+1,1);
end

%%
ylimm = [-90 90];
xlimmd = [flag.DistOnT flag.DistOnT + 300];
[~, idx_xlimd1] = min(abs(time - xlimmd(2)));
xlimd1 = 1:idx_xlimd1+1;
if subplots && plots.delta
    f = figure(1); clf; hold on
    f.Name = 'delta';
    fs1 = subplot(121);
%     plot(time(xlim1),(delta(xlim1,:) - omegastart(xlim1) - delta_s(xlim1,:)));
    plot(time(xlimd1),180/pi*(delta(xlimd1,:) - omegastart(xlimd1) - 0*delta_s(xlimd1,:)),'-');
    axis([xlimmd(1) xlimmd(2) ylimm])
    xlabel('$t$ [s]')
    ylabel('$\delta_{i} - \delta^*$ [deg]')
    fs2 = subplot(122);
    plot(time(xlim2),180/pi*(delta(xlim2,:) - omegastart(xlim2) - 0*delta_s(xlim2,:)));
    axis([xlimm2(1) xlimm2(2) ylimm])
    xlabel('$t$ [s]')
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
    f.Units = 'centimeters';
    f.Position(4) = 7.5;
    fs1.Position(3) = 4.75;
    fs2.Position(3) = 4.75;
    fs2.Position(4) = 5;
    fs1.Position(4) = 5;
    fs2.Position(1) = fs2.Position(1) - 1.25;
    f.Position(3) = 14.25;
    
elseif plots.delta
    f = figure(1); clf; hold on
    f.Name = 'delta';
    %     plot(time,(delta - omegastart - delta_s)); hold on
    plot(time,(delta - omegastart),'--'); hold on
    
    xlabel('$t$ [s]')
    ylabel('$\delta_{i} - \delta^* [deg]$')
    
    axis([-Inf Inf ylim1])
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
end

%%
ylimm = [49.94 50.06];
if subplots && plots.omega
    f = figure(2); clf;
    f.Name = 'omega';
    fs1 = subplot(121);
    plot(time(xlim1),omega(xlim1,:)/2/pi);
    axis([xlimm(1) xlimm(2) ylimm])
    xlabel('Time [s]')
    ylabel('$\omega_i$ [Hz]')
    fs2 = subplot(122);
    plot(time(xlim2),omega(xlim2,:)/2/pi);
    axis([xlimm2(1) xlimm2(2) ylimm])
    xlabel('$t$ [s]')
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
    f.Units = 'centimeters';
    f.Position(4) = 7.5;
    fs1.Position(3) = 4.75;
    fs2.Position(3) = 4.75;
    fs2.Position(4) = 5;
    fs1.Position(4) = 5;
    fs2.Position(1) = fs2.Position(1) - 1;
    f.Position(3) = 14.4;
    
    fs1.YTick = ylimm(1):(ylimm(2)-ylimm(1))/4:ylimm(2);
    fs2.YTick = ylimm(1):(ylimm(2)-ylimm(1))/4:ylimm(2);
    
elseif plots.omega
    f = figure(2); clf;
    f.Name = 'omega';
    plot(time,omega/2/pi);
    xlabel('Time [s]')
    ylabel('$\omega_i$ [Hz]')
    axis([-Inf Inf ylim1])
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
end

%%
ylimm = [-0.2 0.2];
if subplots && plots.xi
    f = figure(3);
    f.Name = 'xi';
    fs1 = subplot(121);
    plot(time(xlim1),xi(xlim1,:));
    axis([xlimm(1) xlimm(2) ylimm])
    xlabel('$t$ [s]')
    ylabel('$\xi_i$')
    fs2 = subplot(122);
    plot(time(xlim2),xi(xlim2,:));
    axis([xlimm2(1) xlimm2(2) ylimm])
    xlabel('$t$ [s]')
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
    f.Units = 'centimeters';
    fs1.Position(3) = 4.5;
    fs2.Position(3) = 4.5;
    fs2.Position(1) = fs2.Position(1) - 1;
    f.Position(3) = 13.5;
elseif plots.xi
    f = figure(3);
    f.Name = 'xi';
    plot(time,xi);
    axis([-Inf Inf ylim1])
    xlabel('$t$ [s]')
    ylabel('$\xi_i$')
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
end

% figure(3)
% plot(time,V);
% % legend('DG1','DG2','DG3','DG4','location','best');
% % title('Voltage $V$')
% xlabel('Time [s]')
% ylabel('Voltage [V]')
% % axis([-Inf Inf 100 400])
% axis([xlimm -Inf Inf])
% colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
%     'TickLabels',string(1:Nplot),...
%     'LineWidth',1,...
%     'Direction','reverse',...
%     'TickLength',0);
% colormap(colvect)

% figure(4)
% plot(time,zeta);
% % legend('DG1','DG2','DG3','DG4','location','best');
% % title('Secondary volt var $\zeta$')
% xlabel('Time [s]')
% ylabel('$\zeta$')
% axis([xlimm -Inf Inf])
% colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
%     'TickLabels',string(1:Nplot),...
%     'LineWidth',1,...
%     'Direction','reverse',...
%     'TickLength',0);
% colormap(colvect)

%% Power
% Plr = repmat(conf.Prat(1)*.9E3,n,length(time))';
% Plr(idx_DistOn:idx_DistOff,1:4:end) = Plr(idx_DistOn:idx_DistOff,1:4:end).*1.1;
if plots.powerinj
    ylimm = [-250 250];
%     Plr = repmat(sys.Plr,length(time),n);
else
    ylimm = [1000 1600];
%     Plr = 0*repmat(sys.Plr,length(time),n);
end


if subplots && plots.P
    f = figure(4); clf;
    f.Name = 'P';
    fs1 = subplot(121);
    if plots.powerinj
        plot(time(xlim1),Pi_final(xlim1,:));
    else
        plot(time(xlim1),P_final(xlim1,:));
    end
    axis([xlimm(1) xlimm(2) ylimm])
    %     axis([xlimm(1) xlimm(2) -Inf Inf])
    xlabel('t [s]')
    ylabel('$P_{I,i}$ [W]')
    fs2 = subplot(122);
    if plots.powerinj
        plot(time(xlim2),Pi_final(xlim2,:));
    else
        plot(time(xlim2),P_final(xlim2,:));
    end
    axis([xlimm2(1) xlimm2(2) ylimm])
    %     axis([xlimm2(1) xlimm2(2) -Inf Inf])
    xlabel('t [s]')
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
    f.Units = 'centimeters';
    f.Position(4) = 7.5;
    fs1.Position(3) = 4.75;
    fs2.Position(3) = 4.75;
    fs2.Position(4) = 5;
    fs1.Position(4) = 5;
    fs2.Position(1) = fs2.Position(1) - 1;
    f.Position(3) = 14.4;
    
    fs1.YTick = ylimm(1):(ylimm(2)-ylimm(1))/4:ylimm(2);
    fs2.YTick = ylimm(1):(ylimm(2)-ylimm(1))/4:ylimm(2);
elseif plots.P
    f = figure(4); clf;
    f.Name = 'P';
    if plots.powerinj
        plot(time,Pi_final);
    else
        plot(time,P_final);
    end
    axis([-Inf Inf ylimP])
    xlabel('t [s]')
    ylabel('$P$ [W]')
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
end

%% 
if plots.Q
    f = figure(8);clf
    f.Name = 'Q';
    plot(time,Q_final);
    % legend('DG1','DG2','DG3','DG4','location','best');
    % title('Reactive Power $Q$')
    xlabel('Time [s]')
    ylabel('$Q$ [VAr]')
    % ylabel('Power Factor')
    axis([-Inf Inf -Inf Inf])
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
end

%%
if plots.Pfactor
    f = figure(8);clf
    f.Name = 'Pfactor';
    
    set(groot,'DefaultAxesUnits','normalized')
    
    Ns = length(N)*length(controllers);
    leg = string(1:Ns);
    for x = 1:Ns
        if mod(x,2)==1 % Controller C1
            %         leg(x) = strcat("C1 N=",);
            colvec = [1 0 0];
        else % Controller C2
            leg(x) = strcat("C2 N=",string(mod(x,length(controllers))));
            colvec = [0 1 0];
        end
        S = sqrt(invChains(x).P.^2 + invChains(x).Q.^2);
        if invChains(x).ctrl.control == "C1"
            plot(invChains(x).time,invChains(x).P./S,...
                '-','color',colvec.*x/Ns);
        elseif invChains(x).ctrl.control == "C2"
            plot(invChains(x).time,invChains(x).P./S,...
                '--','color',colvec.*x/Ns)
        end
        leg(x) = strcat("N=",num2str(invChains(x).sys.n)," ",...
            invChains(x).ctrl.control);
        hold on
    end
    xlabel('Number of agents $n$')
    ylabel('Power Factor P/Q')
    
    axis([-Inf Inf -Inf Inf])
    colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
        'TickLabels',string(1:Nplot),...
        'LineWidth',1,...
        'Direction','reverse',...
        'TickLength',0);
    colormap(colvect)
    % legend(leg)
    
    % Max Power Factor P/S
    fx = figure(9); clf;
    fx.Name = 'PfactorN';
    set(groot,'DefaultAxesUnits','normalized')
    maxPF = zeros(length(N),1);
    for idxN = 1:Ns
        S = sqrt(invChains(idxN).P.^2 + invChains(idxN).Q.^2);
        maxPF(idxN) = min(min(invChains(idxN).P./S));
    end
    if length(controllers)>1
        stem(N(1:end)-.1,maxPF(1:2:end),'-*','r','linewidth',2); hold on
        stem(N(1:end)+.1,maxPF(2:2:end),'--','g','linewidth',2);
        legend('DSS Controller','Standard Controller','location','best')
    else
        stem(N-.1,maxPF(1:end),'--*','r','linewidth',2); hold on
    end
    ylim([0.8 0.86])
    xlabel('Number of agents $n$')
    ylabel('Power Factor P/Q')
    xlim([N(1)-1 N(end)+1])
    xticks(N)
end


pause(0.1)