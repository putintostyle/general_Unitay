global M N r
r = 4; % number of unitary matrices
itnumb = 100; % given iteration number
%% Generating size of unitary matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  we give U_i = m_i*n_i where m_i>=n_i
%%%  First determine the row number of U_i 
%%%  then the column number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_list = [6, 4, 3, 2]; % m_list = [m_1, m_2, ..., m_r] 3\leq m_i\leq 4
n_list = [2, 3, 2, 2];

%%%%%%%%%%%%%%%%%%%%%%%%
M = prod(m_list);
N = prod(n_list);
%% Given an initial state


A = 0.5-rand([M,N]);% initialize A

%% Set experience, 500 trails
trailResHistory = []; % store ultimate res
for trails = 1:500
    UList = {};
    U = 1; % initialize U
    for i = 1:r
        %%%%%%%%%%%%%%% initial guess %%%%%%%%%%%%%%%%%%%
        %%%% size(A0) = (m_i, n_i) and U = U_1\ot U_2 \ot ...\ot U_r
        [U0, S, V] = svds(rand(m_list(i)), n_list(i));
        U = kron(U, U0);
        %%%% append U_i to generating array
        UList{end+1} = U0;
    end
    
    %% Polar decomposition 
    
    %%% Store U_i^p at each iteration p
    UItes = {UList};
    iteP = 1;
    %%% Give U = U_1\ot...U_r
    UNext = U;
    % trailRes = [1/2*norm(A-UNext, 'fro')^2];
    while iteP < itnumb
    
        UTrail =  UItes{end};
        
        for idx = 1:r
            [hGrad_re, hGrad_im, hGrad] = compGrad(A, UTrail, m_list, n_list, idx); % compute the gradient dh/dU_i
            [U_polor, P_polor] = poldec(hGrad); % do polar decomp
            
            UTrail{idx} = U_polor; % renew U_i^p -> U_i^{p+1}
            UNext = 1;
            for j=1:r
                UNext = kron(UNext, UTrail{j}); % compute new U = U_1^{p+1}\ot... U_i^{p+1}\ot U_{i+1}^{p}\ot U_r^p
            end
            
    
        end
    
        UItes{end+1} = UTrail;
        % trailRes = [trailRes, 1/2*norm(A-UNext, 'fro')^2];
        iteP = iteP+1;        
    end
   
            
   trailResHistory = [trailResHistory; 1/2*norm(A-UNext, 'fro')^2] ;
end

%% PLOT       
% first figure: error norm        
figure(1)
scatter(1:500, trailResHistory, '.');

hold on
axscatter = gca; % current axes

[counts,bins] = hist(trailResHistory, 200); %# get counts and bin locations
% Bh = bar(bins,counts,'facecolor','');

h = barh(bins,counts, 'facecolor', 'black');
set(get(h,'Parent'),'xdir','r')



hold off

