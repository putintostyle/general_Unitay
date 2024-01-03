r = 4; 
m_list = randi(6,[1, r]); % m_list = [m_1, m_2, ..., m_r] 1\leq m_i\leq 6
n_list = randi(6,[1, r]); % n_list = [n_1, n_2, ..., n_r] 1\leq n_i\leq 6
%% Given an initial state
QList = [];
rho = 1;
for i = 1:r
    [Q, S, V] = svd(rand(m_list(i)));
    Q = Q(:, 1:n_list(i)); %size(Q) = (m_i, n_i)
    QList = [QList, Q]; %這裡要修，要去找append 的攻勢
    rho = kron(rho, Q);
end

%% Calculate \nabla h

