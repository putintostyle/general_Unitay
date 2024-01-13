global M N r
r = 4; % number of unitary matrices
itnumb = 100; % iteration number
TRnum = 1000; % Repeat trail number
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
for trails = 1:TRnum
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
color = [205, 133, 63]./255; % define color in RGB space
%%%%%%% PLOT SCATTER %%%%%%%
figure(1)

scatter(1:TRnum, trailResHistory, '.');
axscatter = gca;

%%%%%%% PLOT HISTOGRAM %%%%%%%
figure(2)
[counts,bins] = hist(trailResHistory, 100); % get counts and bin locations

h = barh(bins,counts, 'facecolor', color);
set(get(h,'Parent'),'xdir','r') % reverse histogram in x direction
axhist = gca;

%%%%%%% PLOT HISTOGRAM %%%%%%%
figure(3)
hold on
scatter(1:TRnum, trailResHistory, '.');

[counts,bins] = hist(trailResHistory, 100);
h = barh(bins,counts, 'facecolor', color);
set(get(h,'Parent'),'xdir','r');
xticks = axscatter.XTick;

%%%% set xticks %%%
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
% Set the x-tick and y-tick  labels for the second axes
XTick = linspace(min(axhist.XTick), max(axhist.XTick), length(axscatter.XTick));
set(ax2, 'XTickLabel', XTick(end:-1:1),'YTickLabel',axhist.YTick, 'XColor', color, 'YColor', color);
set(ax1, 'XTickLabel', axscatter.XTick(end:-1:1));