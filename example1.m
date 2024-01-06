global M N r
r = 2;
itnumb = 100;
% m_list = randi([3,4],[1, r]); % m_list = [m_1, m_2, ..., m_r] 1\leq m_i\leq 6
% n_list = [];
% for i=1:r
%     n_list = [n_list, randi([2,m_list(i)])];
% end
m_list = [5, 4];
n_list = [3, 2];
 % n_list = [n_1, n_2, ..., n_r] 1\leq n_i\leq 6

M = prod(m_list);
N = prod(n_list);
%% Given an initial state
AList = {}; % generate ground state \rho = A_1\ot A_2...A_r
UList = {}; % generate initial guesses U_1\ot U_2...U_r
A = randn(M,N);
U = 1;
for i = 1:r
    % [A0, S, V] = svds(rand(m_list(i)));
    % A0 = A0(:, 1:n_list(i)); %size(Q) = (m_i, n_i)
    % AList{end+1} = A0; %% append A_i to generating array
    % A = kron(A, A0);

    %%%%%%%%%%%%%%% initial guess %%%%%%%%%%%%%%%%%%%
    [U0, S, V] = svds(rand(m_list(i)), n_list(i));
     %size(Q) = (m_i, n_i)
    UList{end+1} = U0; %% append U_i to generating array
    U = kron(U, U0);
end

%% Calculate \nabla h


%% Polar decomposition
P = 1; % The number of repeated trials.

UItes = {UList};
iteP = 1;
UNext = U;
ResAll = [1/2*norm(A-UNext,'fro')^2];

for idx = 1:r
    disp(UTrail{idx})
end
while iteP < itnumb
    UTrail =  UItes{end};
    

    for idx = 1:r
        [hGrad_re, hGrad_im, hGrad] = compGrad(A, UTrail, m_list, n_list, idx);
        [U_polor, P_polor] = poldec(hGrad);
        
        UTrail{idx} = U_polor;
        UNext = 1;
        for j=1:r
            UNext = kron(UNext, UTrail{j});
        end
        disp(UTrail{idx})
        ResAll = [ResAll, 1/2*norm(A-UNext,'fro')^2];

    end

    UItes{end+1} = UTrail;
    
    iteP = iteP+1;        
end


        
        
figure(1)
semilogy(ResAll','-')

