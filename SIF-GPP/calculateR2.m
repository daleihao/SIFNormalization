function R2 = calculateR2(obs, sim)
filters = ~isnan(obs) & ~isnan(sim) & obs>0 & sim>0;
obs = obs(filters);
sim = sim(filters);

if(sum(filters)<1)
    R2 = nan;
elseif(sum(filters)==1)
    R2 = 1;
else
R = corrcoef(obs, sim);
R = R(2,1);
R2 = R*R;
end