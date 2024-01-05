r = 3; 
m_list = randi([3,4],[1, r]); % m_list = [m_1, m_2, ..., m_r] 1\leq m_i\leq 6
n_list = [];
for i=1:r
    n_list = [n_list, randi([2,m_list(i)])];
end

 % n_list = [n_1, n_2, ..., n_r] 1\leq n_i\leq 6
M = prod(m_list);
N = prod(n_list);
%% Given an initial state
AList = {}; % generate ground state \rho = A_1\ot A_2...A_r
UList = {}; % generate initial guesses U_1\ot U_2...U_r
A = 1;
U = 1;
for i = 1:r
    [A0, S, V] = svd(rand(m_list(i)));
    A0 = A0(:, 1:n_list(i)); %size(Q) = (m_i, n_i)
    AList{end+1} = A0; %% append A_i to generating array
    A = kron(A, A0);

    %%%%%%%%%%%%%%% initial guess %%%%%%%%%%%%%%%%%%%
    [U0, S, V] = svd(rand(m_list(i)));
    U0 = U0(:, 1:n_list(i)); %size(Q) = (m_i, n_i)
    UList{end+1} = U0; %% append U_i to generating array
    U = kron(U, U0);
end

%% Calculate \nabla h
hGrad = {}; % store \partial h/\partial U
for i=2:2
    % reverse U_1 \ot U_2 \ot U_3 --> U_2 \ot U_3 \ot U_1 
    mBackSize = prod(m_list(1:i-1));
    nBackSize = prod(n_list(1:i-1));
    
    mForSize = prod(m_list(i:end));
    nForSize = prod(n_list(i:end));
    
    [BA, S1, S2] = reverse_kron(U, [mBackSize, nBackSize], [mForSize, nForSize]);

    
    
end

