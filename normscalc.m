if calcnorms == 1
    fprintf('Calculating norms...\n')
for idx = 1:length(invChains)
    fprintf('\t%d/%d\n',idx,length(N)*length(controllers))
    n = invChains(idx).sys.n;
    y = invChains(idx).y;
    P = invChains(idx).P;
%     y(:,1:n) = wrapTo2Pi(y(:,1:n));
    time = invChains(idx).time;    
    delta = y(:,1:n);
%     delta = wrapTo2Pi(y(:,1:n));
    omegastart = conf.omegaast*time;
%     omegastart = wrapTo2Pi(conf.omegaast*time);
    xi = y(:,3*n+1:4*n);
    l2norma = zeros(size(y,1),n);                   % ok
    z = reshape(y, size(y,1),n,(size(y,2))/n);      % ok
    [~, idx_DistOff] = min(abs(time - flag.DistOffT));
    for tt = time' % every time step
        idx_t = find(ismember(time,tt));
        for idx_a = 1:n % every agent
            ztemp = z(idx_t,idx_a,:);                  % ok
            ztemp = reshape(ztemp,1,(size(y,2))/n);

            ztemp(2) = ztemp(2) - conf.omegaast;
            ztemp(3) = ztemp(3) - conf.Vast;            

            % disturbance on?                
            if flag.Dist && tt > flag.DistOnT && tt < flag.DistOffT 
%                 [~, idx_DistOn] = min(abs(time - flag.DistOnT));
                ztemp(1) = ztemp(1) - omegastart(idx_t) - ...
                  (delta(idx_DistOff-1,idx_a) - omegastart(idx_DistOff-1));
                ztemp(1) = ztemp(1).*180/pi;
                ztemp(4) = ztemp(4) - xi(idx_DistOff-1,idx_a); 
            elseif flag.Dist && flag.DistOffT==sys.simtime
                ztemp(1) = ztemp(1) - omegastart(idx_t) - ...
                    (delta(idx_DistOn-1,idx_a) - omegastart(idx_DistOn-1));
                ztemp(1) = ztemp(1).*180/pi;
                ztemp(4) = ztemp(4) - xi(end-1,idx_a);
            else
                ztemp(1) = ztemp(1) - omegastart(idx_t) - ...
                    (delta(end-1,idx_a) - omegastart(end-1));
                ztemp(1) = ztemp(1).*180/pi;
                ztemp(4) = ztemp(4) - xi(end-1,idx_a);
            end
%             ztemp(1) = 0;
            ztemp(5) = 0;
            % l2 norm for each agent
            l2norma(idx_t,idx_a) = norm(ztemp);
        end
    end
    % max l2 norm among agents
    l2normt = max(l2norma,[],2);                
    % supreme max l2 
    [linfnorm, ~] = max(max(l2norma));

    invChains(idx).norms.l2 = l2normt;
    invChains(idx).norms.linf = linfnorm;
%         norms=0;
    %%% --- ----------------- --- %%%

    % y = y(find(time>10):end,:);
    % P_final = P_final(find(time>10):end,:);
    % Q_final = Q_final(find(time>10):end,:);
    % time = time(time>10);
end
end
