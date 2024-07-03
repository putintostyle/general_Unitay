global A M N r m_list n_list 
r = 2; % number of unitary matrices
 % given iteration number
J = sqrt(-1);
%% Generating size of unitary matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  we give U_i = m_i*n_i where m_i>=n_i
%%%  First determine the row number of U_i 
%%%  then the column number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_list = [3,  5]; % m_list = [m_1, m_2, ..., m_r] 3\leq m_i\leq 4
n_list = m_list;
% n_list = [];
% for i=1:r
%     n_list = [n_list, randi([3,m_list(i)])];
% end
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
    matrx = rand(m_list(i), n_list(i))+J*rand(m_list(i), n_list(i));
    [A0, S, V] = svd(matrx, "econ");
    
    
    A = kron(A, A0);
    %%%% append A_i to generating array
    AList{end+1} = A0; 

    %%%%%%%%%%%%%%% initial guess %%%%%%%%%%%%%%%%%%%
    %%%% size(A0) = (m_i, n_i) and U = U_1\ot U_2 \ot ...\ot U_r
    matrx = rand(m_list(i), n_list(i))+J*rand(m_list(i), n_list(i));
    [U0, S, V] = svd(matrx, "econ");
    
    U = kron(U, U0);
    %%%% append U_i to generating array
    UList{end+1} = U0;
end
% UList = AList;
%% Gradient Descent
y_init = [];
for i = 1:r
    Uinit = UList{i};
    UiR = real(Uinit);
    UiI = imag(Uinit);
    y_init = [y_init; UiR(:); UiI(:)];
end


T = 1000;
dt = 1e-2;

options = odeset('Abstol',1e-12,'Reltol',1e-12);
% options = odeset(options,'Stats','on');
% options = odeset(options,'Events',@myEventFcn);
dist_array = [];
t_tag = [];
tic
[t,y]= ode15s(@compGrad, [0:T],y_init,options);
toc
%% Extract the solution from y
UHist = {};
ResAll = [];
for time = 1:length(t)
    yn = y(time, :)';

    U_temp = {};
    cum_idx = 1;
    Un = 1;
    for i=1:r
        matrix_size = m_list(i)*n_list(i);
        s_idx = cum_idx;
        e_idx = s_idx+matrix_size-1;
        UR = reshape(yn(s_idx:e_idx), [m_list(i), n_list(i)]);

        s_idx = e_idx+1;
        e_idx = s_idx+matrix_size-1;
        UI = reshape(yn(s_idx:e_idx), [m_list(i), n_list(i)]);

        U_temp{end+1} = UR+sqrt(-1)*UI;
        
        Un = kron(Un, UR+sqrt(-1)*UI);
        cum_idx = e_idx+1;
    end

    UHist{end+1} = U_temp;
    ResAll = [ResAll, 1/2*norm(A-Un, 'fro')^2];
end
%% PLOT       
% first figure: error norm        
figure(1)
loglog(t, ResAll','-')
title('Evolution of residuals','Interpreter','latex','FontSize',15)
xlabel('Update steps','Interpreter','latex','FontSize',12)
ylabel('Residuals, $\frac{1}{2}||A-U_1\otimes\cdots\otimes U_5||_F^2$','Interpreter','latex','FontSize',12)
grid on
% exportgraphics(gcf,'example_1_closest_cpx.eps','Resolution',300);

%%
function gradFun = compGrad(t, y)
    %%%%%%
    % Input:
    % |-A : given target matrix
    % |-y : set of [U_1(:); U_2(:), ...U_r(:)]
    % |-m_list : [m_1, m_2, ...m_r], n_list : [n_1, n_2, ...n_r]
    % |-idx : gradient at idx-th matrix U dh/dU_idx 
    % Output:
    % |-hGrad : A\star(U_1\ot...X_idx..\otU_r)
    global A M N r m_list n_list
    UList = {};
    cum_idx = 1;
    for i=1:r
        
        matrix_size = m_list(i)*n_list(i);
        s_idx = cum_idx;
        e_idx = s_idx+matrix_size-1;
        UR = reshape(y(s_idx:e_idx), m_list(i), n_list(i));

        s_idx = e_idx+1;
        e_idx = s_idx+matrix_size-1;
        UI = reshape(y(s_idx:e_idx), m_list(i), n_list(i));

        UList{end+1} = UR+sqrt(-1)*UI;
        cum_idx = e_idx+1;
    end
    gradFun = [];
    A_tmp = A;
    for idx = 1:r
        U = UList{idx};
   
        if idx>1 % start from U_2
            %%%%%% reverse U_1 \ot U_2 \ot U_3 --> U_2 \ot U_3 \ot U_1 
                
            mBackSize = prod(m_list(1:idx-1));
            nBackSize = prod(n_list(1:idx-1));
                
            mForSize = prod(m_list(idx:end));
            nForSize = prod(n_list(idx:end));
                
            [S1, S2] = reverse_kron( [mForSize, nForSize], [mBackSize, nBackSize]);
            %%%% <A, U \ot V \ot W> = <A, S1 V \ot W \ot U S2'> 
            %%%% S1 from \ot V \ot W <- [mForSize, nForSize] 
            %%%%
            %%%% S2 from U <- [mBackSize, nBackSize] 
            %%%% 
            %%%% <A, U \ot V \ot W> = <A, S1 V \ot W \ot U S2'>  
            %%%%                    = <S1'*A*S2, V\ot(W\ot U)> 
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
        
        %%%% <S1'*A*S2,V\ot(W\ot U)>_R 
        %%%% Let block[S1'*A*S2] = \tildeA; (W\ot U) = WU
        %%%%  = <inner((\tildeA)_r, (WU)_r), V_r> - <inner((\tildeA)_r, (WU)_i), V_i>
        %%%%   + <inner((\tildeA)_i, (WU)_r), V_i> + <inner((\tildeA)_i, (WU)_i), V_r>
        
        %%%%  = <inner((\tildeA)_r, (WU)_r)+inner((\tildeA)_i, (WU)_i), V_r>
        %%%%   + <inner((\tildeA)_i, (WU)_r)-inner((\tildeA)_r, (WU)_i), V_i>
        
        %%%% d/dV_r = inner((\tildeA)_r, (WU)_r)+inner((\tildeA)_i, (WU)_i) 
        %%%%        = real(inner(cong(\tildeA), WU)
        
        %%%% d/dV_i = inner((\tildeA)_i, (WU)_r)-inner((\tildeA)_r, (WU)_i) 
        %%%%        = -imag(inner(cong(\tildeA), WU)
        
        A_sep_conj = cellfun(@(a) conj(a), Asep, 'UniformOutput', 0);
        
        dh_re = cellfun(@(a) real(sum(a.*backUni, 'all')), A_sep_conj, 'UniformOutput', 1);
        dh_im = -1*cellfun(@(a) imag(sum(a.*backUni, 'all')), A_sep_conj, 'UniformOutput', 1);
        
        hGrad_re = dh_re;
        hGrad_im = dh_im;
        hGrad = dh_re+sqrt(-1).*dh_im;
        
        hGrad = U*(((U'*hGrad)-((U'*hGrad))')/2);
        

        hGrad_re = real(hGrad);
        hGrad_im = imag(hGrad);
        gradFun = [gradFun; hGrad_re(:); hGrad_im(:)];
    end
end
