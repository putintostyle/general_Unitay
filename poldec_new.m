function [U, H] = poldec_new(A)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %   This code is to calculate the polar decomposition of a matrix
        % %   [U H] = poldec(A) factorizes a matrix A \in C^{m\times n}, m\geq n
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [P1,Sigma,Q] = svd(A,'econ');
        H = Q*Sigma*Q';
        U = P1*Q';

    end