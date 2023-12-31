function [U, H] = poldec(A)
% %
% %   Calculate the polar decomposition of a matrix.
% %   [U H] = poldec(A) factorizes a full-rank matrix A \in C^{m\times n}, m\geq n 
% %   such that A=U*H, 
% %       (Note that H = (A'*A)^{1/2} is uniquely determined even if A is singular.)
% %   where
% %    A = P*blkdiag(Sigma'; 0 )'*Q' is the singular value decomposition of A,
% %        Sigma \in C^{n\times n} and 0 \in C^{m-n\times n}
% %        P= [P1,P2]\in C^{m\times m} and Q \in C^{n\times n} are unitary,
% %        where P1 \in C^{m\times n} and P2 \in C^{m\times m-n},
% %              P1'*P1 = In,
% %    That is, A = P1*Sigma*Q'
% %    H = Q*Sigma*Q' \in C^{n\times n}, and   
% %    U = P1*Q'.
% %
C = A'*A;   % C = A'*A = Q*Sigma^2*Q'
[Q, lambda] = eig(C);
Sigma = sqrt(diag((lambda))); 

Hinv = repmat(1./Sigma',size(A,2),1).*Q*Q'; % Hinv = Q*(1./Sigma)*Q' 
U = A*Hinv;
H = U'*A;




