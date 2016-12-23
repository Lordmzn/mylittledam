
param = param(:); % it must be a vector
for useless = 1:n_experiments
    storage = NaN(2401,1);
    stepcosts = NaN(2401,4);
    storage(1) = 100;
    inflow = 40;
    D = 1;
    N = 3;
    phi = NaN(N,1);
    for t = 1:2400
        hh = storage(t) / 155;
        count = 1;
        for j = 1:N
            bf = 0;
            bf = bf + (hh - param(count))^2 / param(count+1)^2;
            phi(j) = exp(-bf);
            count = count + 2;
        end
        % weighted sum + bias
        decision = sum(param(count:end-1) .* phi) + param(end);
        decision = decision * 155;
        max_rel = storage(t);
        min_rel = max(storage(t) - 100, 0);
        release = min(max_rel, max(min_rel, decision));
        storage(t+1) = storage(t) + inflow - release;
        stepcosts(t+1, 1) = max(50 - release, 0);
        stepcosts(t+1, 2) = max(storage(t+1) - 50, 0);
        hp = 1 * 9.81 * 1000 / 3600000 * storage(t + 1) * max(release, 0);
        stepcosts(t+1, 3) = max(4.36 - hp, 0);
        stepcosts(t+1, 4) = max(release - 30, 0);
    end
    
    disp(mean(stepcosts(~isnan(stepcosts(:,1)), :)))
end