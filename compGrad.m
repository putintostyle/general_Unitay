function [hGrad_re, hGrad_im, hGrad] = compGrad(A, UList, m_list, n_list, idx)
%%%%%%
% Input:
% |-A : given target matrix
% |-UList : set of U_1, U_2, ...U_r
% |-m_list : [m_1, m_2, ...m_r], n_list : [n_1, n_2, ...n_r]
% |-idx : gradient at idx-th matrix U dh/dU_idx 
% Output:
% |-hGrad : A\star(U_1\ot...X_idx..\otU_r)
global M N r

A_tmp = A;
 
if idx>1 % start from U_2
    %%%%%% reverse U_1 \ot U_2 \ot U_3 --> U_2 \ot U_3 \ot U_1 
        
    mBackSize = prod(m_list(1:idx-1));
    nBackSize = prod(n_list(1:idx-1));
        
    mForSize = prod(m_list(idx:end));
    nForSize = prod(n_list(idx:end));
        
    [S1, S2] = reverse_kron( [mForSize, nForSize], [mBackSize, nBackSize]);
    %%%% <A, U \ot V \ot W> = <A, S1 \ot V \ot W \ot U S2'> 
    %%%% S1 from \ot V \ot W <- [mForSize, nForSize] 
    %%%%
    %%%% S2 from U <- [mBackSize, nBackSize] 
    %%%% 
    %%%% <A, U \ot V \ot W> = <A, S1 \ot V \ot W \ot U S2'>  
    %%%%                    = <S1'*A*S2,V\ot(W\ot U)>
    
    A_tmp = S1'*A*S2;
    
end

%%%%%%%%% <A, U_1\ot....\ot U_n> = <\tilda A(U_1\ot...\U^_i\ot U_n), U_i>
%%%%%%%%% separate A to cell form
rowDist = ones(1, m_list(idx)).*(M/m_list(idx));
colDist = ones(1, n_list(idx)).*(N/n_list(idx));
    

Asep = mat2cell(A_tmp, rowDist, colDist);

%%%%%%%%% compute A(U_1\ot...\U^_i\ot U_n)
backUni = 1;

if idx+1<=r
    for j=idx+1:r
        backUni = kron(backUni, UList{j});
    end
end

if idx-1>=1
    for j=1:idx-1
        backUni = kron(backUni, UList{j});
    end
end

%%%% <S1'*A*S2,V\ot(W\ot U)> = <\inner(block[S1'*A*S2], (W\ot U)), V>

A_sep_conj = cellfun(@(a) conj(a), Asep, 'UniformOutput', 0);

dh = cellfun(@(a) sum(a.*backUni, 'all'), A_sep_conj, 'UniformOutput', 1);
    

hGrad_re = 2*real(dh);
hGrad_im = 2*imag(dh);
hGrad = dh;

