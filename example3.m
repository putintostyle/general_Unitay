function example3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In this example, we are working on a Toffoli gate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear, close All, clc
%global M N r


itnumb = 10; % iteration number
I = sqrt(-1);
TRnum = 1; % Repeat trail number
%% Generating size of unitary matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  we give U_i = m_i*n_i where m_i>=n_i
%%%  First determine the row number of U_i 
%%%  then the column number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% %%% [4\times 2]
% r = 2; % number of unitary matrices
% m_list = [4, 2]; % m_list = [m_1, m_2, ..., m_r] 3\leq m_i\leq 4
% n_list = [4, 2];
%%% [2\times 2 \times 2]
r = 3; % number of unitary matrices
m_list = [2, 2, 2]; % m_list = [m_1, m_2, ..., m_r] 3\leq m_i\leq 4
n_list = [2, 2, 2];

%%%%%%%%%%%%%%%%%%%%%%%%
M = prod(m_list);
N = prod(n_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define the initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% matrix A: Controlled controlled X gate (CCNOT or Toffoli gate) 
%%%  Create a controlled controlled X gate that acts on two control qubits 
%%%  with indices 1 and 2 and a target qubit with index 3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A = ccxGate(1,2,3); % Available for Matlab R2023a
A = [ 
     1     0     0     0     0     0     0     0;
     0     1     0     0     0     0     0     0;
     0     0     1     0     0     0     0     0;
     0     0     0     1     0     0     0     0;
     0     0     0     0     1     0     0     0;
     0     0     0     0     0     1     0     0;
     0     0     0     0     0     0     0     1;
     0     0     0     0     0     0     1     0];
A0 = A;

%% Set experience, 100 trails
trailResHistory = []; % store ultimate res with A;
trailResOrig = []; % store ultimate res with A0;
for trails = 1:TRnum
    UList = {}; % generate initial guesses U_1\ot U_2...U_r
    U = 1; % initialize U
    for j = 1:r
        %%%%%%%%%%%%%%% initial guess %%%%%%%%%%%%%%%%%%%
        %%%% size(A0) = (m_i, n_i) and U = U_1\ot U_2 \ot ...\ot U_r
        [U0, S, V] = svds(rand(m_list(j)), n_list(j));
        U = kron(U, U0);
        %%%% append U_i to generating array
        UList{end+1} = U0;
    end
    
    %% Polar decomposition 
    
    %%% Store U_i^p at each iteration p
    UItes = {UList};
    iteP = 1;
    %%% Give U = U_1\ot...U_r
    UNext = U;
    % trailRes = [1/2*norm(A-UNext, 'fro')^2];
    while iteP < itnumb
    
        UTrail =  UItes{end};
        
        for idx = 1:r
            [hGrad_re, hGrad_im, hGrad] = compGrad0(A, UTrail, m_list, n_list, idx, M, N, r); % compute the gradient dh/dU_i
            [U_polor, P_polor] = poldec0(hGrad); % do polar decomp
            
            UTrail{idx} = U_polor; % renew U_i^p -> U_i^{p+1}
            UNext = 1;
            for j=1:r
                UNext = kron(UNext, UTrail{j}); % compute new U = U_1^{p+1}\ot... U_i^{p+1}\ot U_{i+1}^{p}\ot U_r^p
            end
            
    
        end
    
        UItes{end+1} = UTrail;
        % trailRes = [trailRes, 1/2*norm(A-UNext, 'fro')^2];
        iteP = iteP+1;        
    end

   trailResHistory = [trailResHistory; norm(A-UNext, 'fro')^2] ;
   trailResOrig = [trailResOrig;norm(A0-UNext, 'fro')^2] ;        
  
end

%% PLOT       
% first figure: error norm    
% color = [205, 133, 63]./255; % define color in RGB space
%%%%%%% PLOT SCATTER %%%%%%%
figure(11)
plot(1:TRnum, trailResHistory, 'o');
xlabel('number of test','Interpreter','latex','FontSize',18)
ylabel('perturbation residuals','Interpreter','latex','FontSize',18)
 

figure(12)
[counts,bins] = hist(trailResHistory, 10);
% hist(trailResHistory, 10)
bar(bins(1:end), counts);
title('Collection of Pertubation Residuals ($\|\widehat{A} - \widehat{U}^{(1)}\otimes\cdots\otimes \widehat{U}^{(5)}\|_{F}^{2}$)','Interpreter','latex','FontSize',20)
xlabel('perturbation residuals','Interpreter','latex','FontSize',18)
ylabel('frequency of occurrence','Interpreter','latex','FontSize',18)
xlim([min(bins)-1e-15, max(bins(end))+1e-15]);
 

figure(21)
plot(1:TRnum, trailResOrig, 'o');
xlabel('number of test','Interpreter','latex','FontSize',18)
ylabel('residuals','Interpreter','latex','FontSize',18)
 


figure(22)
[counts,bins] = hist(trailResOrig, 10);
%hist(trailResOrig, 10)
bar(bins(1:end), counts);
title(' Collection of Residuals ($\|A - \widehat{U}^{(1)}\otimes\cdots\otimes \widehat{U}^{(5)}\|_{F}^{2}$)','Interpreter','latex','FontSize',20)
xlabel('residuals','Interpreter','latex','FontSize',18)
ylabel('frequency of occurrence','Interpreter','latex','FontSize',18)
saveas(gcf,'example_2bar_Ori1e_1')
saveas(gcf,'example_2bar_Ori1e_1','epsc')


% % %%%%%%% PLOT HISTOGRAM %%%%%%%
% % figure(2)
% % [counts,bins] = hist(trailResHistory, 100); % get counts and bin locations
% % % Computes the histogram of the data in the vector trailResHistory using 100 bins. 
% % % counts: the number of elements in each bin
% % % bins: contain the bin locations.
% % 
% % h = barh(bins,counts, 'facecolor', color);
% % %barh: Creates a horizontal bar chart using the bin locations (bins) as the y-coordinates 
% % % and the counts (counts) as the corresponding bar lengths. 
% % % The 'facecolor', color argument sets the color of the bars.
% % 
% % set(get(h,'Parent'),'xdir','r') % reverse histogram in x direction
% % % Reverses the direction of the x-axis in the histogram, making it go from right to left. 
% % % This is achieved by accessing the parent of the bar chart (get(h,'Parent')), 
% % % which is the axes object, and then setting the x-direction ('xdir') to 'r' (reverse).
% % xlabel('Trails')
% % ylabel('Residual')
% % axhist = gca;
% % 
% % 
% % 
% % %%%%%%% PLOT HISTOGRAM %%%%%%%
% % figure(3)
% % hold on
% % scatter(1:TRnum, trailResHistory, '.');
% % 
% % [counts,bins] = hist(trailResHistory, 100);
% % h = barh(bins,counts, 'facecolor', color);
% % set(get(h,'Parent'),'xdir','r');
% % xticks = axscatter.XTick;
% % %%%% set xticks %%%
% % ax1 = gca;
% % ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
% % set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% % set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
% % set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
% % % Set the x-tick and y-tick  labels for the second axes
% % XTick = linspace(min(axhist.XTick), max(axhist.XTick), length(axscatter.XTick));
% % set(ax2, 'XTickLabel', XTick(end:-1:1),'YTickLabel',axhist.YTick, 'XColor', color, 'YColor', color);
% % set(ax1, 'XTickLabel', axscatter.XTick(end:-1:1));
% % ax1.XLabel.String = "Trails" ;
% % ax1.YLabel.String = "Residual" ;
% % ax2.YLabel.String = "Residual" ;
% % ax2.XLabel.String = "Frequency" ;
% % exportgraphics(gcf,'example_2_closest.eps','Resolution',300);
% % 
% % %  saveas(gcf,'example_2_closest')
% % %  saveas(gcf,'example_2_closest','epsc')



    function [hGrad_re, hGrad_im, hGrad] = compGrad0(A, UList, m_list, n_list, idx, M, N, r)
        %%%%%%
        % Input:
        % |-A : given target matrix
        % |-UList : set of U_1, U_2, ...U_r
        % |-m_list : [m_1, m_2, ...m_r], n_list : [n_1, n_2, ...n_r]
        % |-idx : gradient at idx-th matrix U dh/dU_idx
        % Output:
        % |-hGrad : A\star(U_1\ot...X_idx..\otU_r)
        % global M N r

        A_tmp = A;

        if idx>1 % start from U_2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% The following code subroutine is to prepare perfect 
            %%%% shuffle entries S1 and S2 for the following action:
            %%%% For example r = 3 and idx = 2.
            %%%% <A, U \ot (V \ot W)> = <A, S1{(V \ot W)\ot U}S2'>
            %%%%                      = <S1'*A*S2,(V\ot W)\ot U>
            %%%% S1 for V \ot W: [mForSize, nForSize]
            %%%%
            %%%% S2 for U: [mBackSize, nBackSize]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mBackSize = prod(m_list(1:idx-1));
            nBackSize = prod(n_list(1:idx-1));

            mForSize = prod(m_list(idx:end));
            nForSize = prod(n_list(idx:end));

            [S1, S2] = reverse_kron0( [mForSize, nForSize], [mBackSize, nBackSize]);


            A_tmp = S1.'*A*S2; 
            
        end
        
        %%%%%%%%% <A, U_1\ot....\ot U_n> = <\tilda A(U_1\ot...\U^_i\ot U_n), U_i>
        
        % Fold A into a cell matrix with each cell of size
        % (M/m_list(idx))\times (N/n_list(idx)).
        rowDist = ones(1, m_list(idx)).*(M/m_list(idx));
        colDist = ones(1, n_list(idx)).*(N/n_list(idx));
        Asep = mat2cell(A_tmp, rowDist, colDist);

        
        %%% compute A(U_1\ot...\U^_i\ot U_n)
        backUni = 1;
        
        % compute U_{idx+1} \ot...\ot U_r
        if idx+1<=r
            for j=idx+1:r
                backUni = kron(backUni, UList{j});
            end
        end
        
        % compute U_{1} \ot...\ot U_{idx-1}\ot U_{idx+1} \ot...\ot U_r
        if idx-1>=1
            for j=1:idx-1
                backUni = kron(backUni, UList{j});
            end
        end
        
        %%%% <S1'*A*S2,V\ot(W\ot U)> = <\inner(block[S1'*A*S2], (W\ot U)), V>

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % A_sep_conj saves the complex conjugate of each element in the cell array
        % of Asep.
        % 'UniformOutput', 0: This specifies that the output of cellfun
        %  should be a cell array even if all elements are scalar. Without this,
        % MATLAB would attempt to concatenate the results into an array if
        % they are all scalars, which might cause issues when dealing with
        % different data types.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_sep_conj = cellfun(@(a) conj(a), Asep, 'UniformOutput', 0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The following is to calculate the action of
        % \overline{A}\circledast_\ell (U^{(1)}\otimes \ldots \otimes
        % \hat{X}^{(\ell)} \otimes \ldots \otimes U^{(r)}).
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sum(A_sep_conj.*backUni, 'all') takes a single argument "A_sep_conj",
        % multiplies each element of "A_sep_conj" element-wise.
        % 'UniformOutput', 1: This specifies that the output of cellfun
        % should be an array (matrix) since each result is a scalar. Without this,
        % MATLAB would attempt to concatenate the scalar results
        % into an array.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dh_re = cellfun(@(a) real(sum(a.*backUni, 'all')), A_sep_conj, 'UniformOutput', 1);
        dh_im = -1*cellfun(@(a) imag(sum(a.*backUni, 'all')), A_sep_conj, 'UniformOutput', 1);

        hGrad_re = dh_re;
        hGrad_im = dh_im;
        hGrad = dh_re+sqrt(-1).*dh_im;
        
    end




    function [S_1, S_2] = reverse_kron0(ASize, BSize)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This function is to find permutation matrices S_1 and S_2 such 
        % that  B\otimes A = S_1 (A\otimes B) S_2^\top, 
        % where A: m1xn1; B: m2xn2; S_1: m1xm2, S_2: n1xn2.
        % Reference: https://en.wikipedia.org/wiki/Kronecker_product
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m1 = ASize(1);
        n1 = ASize(2);

        m2 = BSize(1);
        n2 = BSize(2);


        I_m = eye(m1*m2);
        I_n = eye(n1*n2);
        %% This can be improved
        S_1 = [];
        for i=1:m2
            S_1 = [S_1;I_m(i:m2:m1*m2, :)];
        end

        S_2 = [];
        for i=1:n2
            S_2 = [S_2;I_n(i:n2:n1*n2, :)];
        end

        

    end






    function [U, H] = poldec0(A)
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


save('0220Toffoli.mat')

end