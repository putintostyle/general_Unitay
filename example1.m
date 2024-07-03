global M N r
r = 3; % number of unitary matrices
itnumb = 10; % given iteration number
%% Generating size of unitary matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  we give U_i = m_i*n_i where m_i>=n_i
%%%  First determine the row number of U_i 
%%%  then the column number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m_list = randi([4,9],[1, r]); % m_list = [m_1, m_2, ..., m_r] 3\leq m_i\leq 4
% n_list = [];
% for i=1:r
%     n_list = [n_list, randi([3,m_list(i)])];
% end
m_list = [2, 2, 2]; % m_list = [m_1, m_2, ..., m_r] 3\leq m_i\leq 4
n_list = [2, 2, 2];
%%%%%%%%%%%%%%%%%%%%%%%%
M = prod(m_list);
N = prod(n_list);
%% Given an initial state
AList = {}; % generate ground state \rho = A_1\ot A_2...A_r
UList = {}; % generate initial guesses U_1\ot U_2...U_r
A = 1; % initialize A
U = 1; % initialize U
for i = 1:r
    %%%% size(A0) = (m_i, n_i) and A = A_1\ot A_2 \ot ...\ot A_r
    [A0, S, V] = svds(rand(m_list(i)), n_list(i));
    A = kron(A, A0);
    %%%% append A_i to generating array
    AList{end+1} = A0; 

    %%%%%%%%%%%%%%% initial guess %%%%%%%%%%%%%%%%%%%
    %%%% size(A0) = (m_i, n_i) and U = U_1\ot U_2 \ot ...\ot U_r
    [U0, S, V] = svds(rand(m_list(i)), n_list(i));
    U = kron(U, U0);
    %%%% append U_i to generating array
    UList{end+1} = U0;
end
A = [ 
     1     0     0     0     0     0     0     0;
     0     1     0     0     0     0     0     0;
     0     0     1     0     0     0     0     0;
     0     0     0     1     0     0     0     0;
     0     0     0     0     1     0     0     0;
     0     0     0     0     0     1     0     0;
     0     0     0     0     0     0     0     1;
     0     0     0     0     0     0     1     0];
%% Polar decomposition
P = 1; % The number of repeated trials.

%%% Store U_i^p at each iteration p
UItes = {UList};
iteP = 1;
%%% Give U = U_1\ot...U_r
UNext = U;
ResAll = [1/2*norm(A-UNext,'fro')^2];


while iteP < itnumb

    UTrail =  UItes{end};
    
    for idx = 1:r
        [hGrad_re, hGrad_im, hGrad] = compGrad(A, UTrail, m_list, n_list, idx); % compute the gradient dh/dU_i
        [U_polor, P_polor] = poldec_new(hGrad); % do polar decomp
        
        UTrail{idx} = U_polor; % renew U_i^p -> U_i^{p+1}
        UNext = 1;
        for j=1:r
            UNext = kron(UNext, UTrail{j}); % compute new U = U_1^{p+1}\ot... U_i^{p+1}\ot U_{i+1}^{p}\ot U_r^p
        end
        
        ResAll = [ResAll, 1/2*norm(A-UNext,'fro')^2];

    end

    UItes{end+1} = UTrail;
    
    iteP = iteP+1;        
end


%% PLOT       
% first figure: error norm        
figure(1)
semilogy(ResAll','-')
title('Evolution of residuals','Interpreter','latex','FontSize',15)
xlabel('Update steps','Interpreter','latex','FontSize',12)
ylabel('Residuals, $\frac{1}{2}||A-U_1\otimes\cdots\otimes U_3||_F^2$','Interpreter','latex','FontSize',12)
grid on
% exportgraphics(gcf,'example_1_closest.eps','Resolution',300);


