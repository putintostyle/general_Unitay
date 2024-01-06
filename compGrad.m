function [hGrad_re, hGrad_im, hGrad] = compGrad(A, UList, m_list, n_list, idx)
global M N r

A_tmp = A;
 
if idx>1
    %%%%%% reverse U_1 \ot U_2 \ot U_3 --> U_2 \ot U_3 \ot U_1 
        
    mBackSize = prod(m_list(1:idx-1));
    nBackSize = prod(n_list(1:idx-1));
        
    mForSize = prod(m_list(idx:end));
    nForSize = prod(n_list(idx:end));
        
    [S1, S2] = reverse_kron([mBackSize, nBackSize], [mForSize, nForSize]);
    A_tmp = S1'*A*S2;
    
end

%%%%%%%%% <A, U_1\ot....\ot U_n> = <A(U_1\ot...\U^_i\ot U_n), U_i>
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
    
A_sep_conj = cellfun(@(a) conj(a), Asep, 'UniformOutput', 0);
dh = cellfun(@(a) sum(a.*backUni, 'all'), A_sep_conj, 'UniformOutput', 1);
    
    % dh_conj = cellfun(@(a) sum(dot(a,conj(backUni))), conj(Asep), 'UniformOutput', 1);
    
% Wrintinger
hGrad_re = 2*real(dh);
hGrad_im = 2*imag(dh);
hGrad = dh;

