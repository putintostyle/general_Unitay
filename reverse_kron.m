function [S_1, S_2] = reverse_kron( ASize, BSize)

%%
% B\kron A = S_1 A\kron B S_2'
% S_1 = M_m1xm2, S_2 = M_n1xn2
m1 = ASize(1);
n1 = ASize(2);

m2 = BSize(1);
n2 = BSize(2);


I_m = eye(m1*m2);
I_n = eye(n1*n2);
%% This can bbe improved
S_1 = [];
for i=1:m2
    S_1 = [S_1;I_m(i:m2:m1*m2, :)];
end

S_2 = [];
for i=1:n2
    S_2 = [S_2;I_n(i:n2:n1*n2, :)];
end

% BA = S_1*AB*S_2';