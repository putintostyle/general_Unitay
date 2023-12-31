
%%Plb = 5; Pub = 5;


%%for p = Plb:Pub
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Initial Setup
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     m1 = 20; n1 = 15;
% %     m2 = 20; n2 = 15;
    m1 = 5; n1 = 3;
    m2 = 4; n2 = 2;
    
    %m1=3; n1=2; m2=5; n2=4;
    
    
    m1m2 = m1*m2;
    n1n2 = n1*n2;
    
    
    itnumb = 100; % Number of required iterations 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Generate the initial matrix A.
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% % %     %%% [Case 1]
% % %     %%%
% % %     [Q1ini,tmp,tmp1] = svds(randn(m1),n1);
% % %     [Q2ini,tmp,tmp2] = svds(randn(m2),n2);
% % % 
% % %     A = kron(Q1ini,Q2ini);
% % %     normA = norm(A,'fro');

    
    %%% [Case 2]
    %%%
    A = randn(m1m2,n1n2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save R(A) \in \mathbb{R}^{m1n1\times m2n2}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rowDist = m2*ones(1,m1);
    colDist = n2*ones(1,n1);
    %%%
    %%% Divide A = [Aij] with Aij\in \mathbb{R}^{m2\times n2}
    %%%
    C = mat2cell(A, rowDist, colDist); 
    %%%
    %%% Change C into a matrix of m2-by-n2*(m1*n1) 
    %%%
    C = C(:)'; 
    %%%
    %%% reshape(cell2mat(C),m2*n2,[]): 
    %%% m2*n2 is to vec(Aij) and there are m1*n1 vectors of size m2*n2
    %%%
    RA = reshape(cell2mat(C),m2*n2,[])'; 

    P = 1; % The number of repeated trials.
    ResAll = [];
    
    for trial = 1:P
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% Generate an initial matrix Q2.
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Q2,S,V] = svds(randn(m2),n2);
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Initialize parameters for saving information.
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gQ2All = [];
        gQ1Q2All = []; Q1all = []; Q2all = []; 
        TimeRec = []; IteRec = [];
        
        iteP = 0;
        g2 = 1e5;
        tic
%        while g2 > 1e-8*(1+normA)% norm(Q02-Q12,'fro') > 1e-8*(norm(Q02,'fro'))
        while iteP < itnumb    
            
            ind = 2;
            [RAQ] = AStarmnQ(RA,Q2,m2,m1,n2,n1,ind);
            [Q1, P1] = poldec(RAQ);
            
            g1 = 1/2*norm(A-kron(Q1,Q2),'fro')^2;
            gQ1Q2All = [gQ1Q2All,g1];
            
            
            ind = 1;
            [RAQ] = AStarmnQ(RA,Q1,m2,m1,n2,n1,ind);
            [Q2, P2] = poldec(RAQ);
            
            g2 = 1/2*norm(A-kron(Q1,Q2),'fro')^2;
            
            gQ2All = [gQ2All,g2];
            gQ1Q2All = [gQ1Q2All,g2];
            
            
            iteP = iteP+1;
            Q1all = [Q1all, Q1(:)];
            Q2all = [Q2all, Q2(:)];
            
            
        end
        
        
        TimeFP = toc; 
        
        
        ResAll = [ResAll;gQ1Q2All];
        
        figure(1)
        semilogy(gQ1Q2All)
        
        
        figure(2)
        semilogy(ResAll','-')
        
        TimeRec = [TimeRec, TimeFP]
%        IteRec = [IteRec,iteP]
        
    
        
        
        
    end
    
%%end




%% Intermediate postprocessing
%
% 1) Plot the monotone decreasing gQ1Q2 values

figure(1), semilogx(gQ2All,'-');
title('Evolution of residuals','Interpreter','latex','FontSize',20)
xlabel('number of iterations','Interpreter','latex','FontSize',18)
ylabel('residuals','Interpreter','latex','FontSize',18)


figure(12), semilogx(gQ1Q2All,'-');
title('Evolution of residuals','Interpreter','latex','FontSize',20)
xlabel('number of iterations','Interpreter','latex','FontSize',18)
ylabel('residuals','Interpreter','latex','FontSize',18)

figure(2), semilogx(Q1all','-');
title('Evolution of $Q^{(1)}$','Interpreter','latex','FontSize',20)
xlabel('number of iterations','Interpreter','latex','FontSize',18)
ylabel('entries of $Q^{(1)}$','Interpreter','latex','FontSize',18)





figure(3), semilogx(Q2all','-');
title('Evolution of $Q^{(2)}$','Interpreter','latex','FontSize',20)
xlabel('number of iterations','Interpreter','latex','FontSize',18)
ylabel('entries of $Q^{(2)}$','Interpreter','latex','FontSize',18)





savefile = '0120E1PM';
%savefile = [savefile,num2str(p)];
save(savefile)