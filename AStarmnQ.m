function [RAQ] = AStarmnQ(RA,Q,m2,m1,n2,n1,ind)
% %
% % This code is to calculate 
% %      [<A_{ij}, Q2>] \in \mathbb{R}^{m1 \times n1} if ind == 2
% % or 
% %      [<tilde{A}_{ij}, Q1>] \in \mathbb{R}^{m2 \times n2} if ind == 1

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%% Save R(A) \in \mathbb{R}^{m1n1\times m2n2}
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % rowDist = m2*ones(1,m1);
% % colDist = n2*ones(1,n1);
% % C = mat2cell(A, rowDist, colDist); C = C(:)';
% % RA = reshape(cell2mat(C),m2*n2,[])';

if ind == 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate [<A_{ij}, Q2>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q2 = Q(:);
    RAQ = reshape(RA*Q2, m1, n1);
  
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate [<tilde{A}_{ij}_{ij}, Q1>], where R(tilde{A}) = R(A)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Q1 = Q(:);
    RAQ = reshape(RA'*Q1, m2, n2);
end

