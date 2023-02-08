z%%
n =10;

N = rand(n,1);

%%
evol = [];
for gamma  = 0:0.001:0.4
    disp(gamma);
    cvx_begin quiet
    variable W(n,1)
    minimize ( 0.5*(W - N)'*(W - N) + gamma*sum(abs((W(2:end)-W(1:end-1)))) );
    cvx_end
    
    evol = [evol, W];
end

%%
gamma = 0;
s = sign(diff(N));  % this contains the difference in sign between consecutive entries
sjumps = s;         % this contains the difference in sign between consecutive blocks
R = ones(n,1);      % this will contain all the averages of each block
Istart = nan(n,1);  % for each node we will know the start and end of its block
Iend = nan(n,1);

N_ave = N;
%%
disp(gamma);

% compute averages from s
tmps = [s;0];
ix_start = 1;
ix_end = 1;
for ix = 1:n-1
    if (s(ix) == 0)
        ix_end = ix + 1;
    else
        N_ave(ix_start:ix_end) = mean(N(ix_start:ix_end));
        R(ix_start:ix_end) = ix_end - ix_start + 1;
        sjumps(ix_start:ix_end-1) = tmps(ix_end);
        Istart(ix_start:ix_end) = ix_start;
        Iend(ix_start:ix_end) = ix_end;
        ix_start = ix+1;
        ix_end = ix+1;
    end
end
N_ave(ix_start:ix_end) = mean(N(ix_start:ix_end));
R(ix_start:ix_end) = ix_end - ix_start + 1;
sjumps(ix_start:ix_end-1) = tmps(ix_end);
Istart(ix_start:ix_end) = ix_start;
Iend(ix_start:ix_end) = ix_end;


%% compute the new alg parameters

block_ix_start = unique(Istart);
block_ix_end = unique(Iend);
new_ix = -1;
gamma_new = inf;
for i = 2:length(block_ix_start)

    
    ix_p = block_ix_start(i);
    ix_m = block_ix_start(i-1);
    
    if (i > 2)
        ix_mm = block_ix_start(i-2);
        sjumps_mm = sjumps(ix_mm);
    else
        sjumps_mm = 0;
    end
    
    if (ix_p < n)
        sjumps_p = sjumps(ix_p);
    else
        sjumps_p = 0;
    end
        
    
    d2 = (sjumps_p - sjumps(ix_m)) / R(ix_p);
    d1 = (sjumps(ix_m) - sjumps_mm) / R(ix_m);
    
    gamma_pot = - (N_ave(ix_p) - N_ave(ix_m)) / (d2 - d1);
     
    if (gamma_pot > gamma && gamma_pot < gamma_new && ~isnan(gamma_pot) && gamma_pot > -inf)
        gamma_new = gamma_pot;
        new_ix = block_ix_end(i);
    end
     
end
%%
s(new_ix)=0;
gamma = gamma_new;

% 
% %%
% 
% sfR = diff([0;sjumps;0]);
% sfR = sfR./R;
% sfR = -diff(sfR);
% 
% pot_gam = diff(N_ave)./(sfR);
% disp(pot_gam');
% 
% pot_gam = max([pot_gam,gamma*ones(n-1,1)],[],2);
% pot_gam = pot_gam + 10000*(pot_gam == gamma);
% 
% [gamma, new_i] = min(pot_gam);
% 
% s(new_i) = 0;
% 
% 
% 
