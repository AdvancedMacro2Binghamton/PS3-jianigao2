% PROGRAM NAME: ps4huggett.m
clear ;
close all;

% PARAMETERS
beta = .9932; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix


% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;%upper bound of grid points
num_a = 600;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.98;
q_max = 1;

% ITERATE OVER ASSET PRICES
aggsav = 1 ;
while abs(aggsav) >= 0.01 ;
    q_guess = (q_min + q_max) / 2;

    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.0001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
      
        vf=bsxfun(@plus,ret,permute(beta*PI*v_guess,[3,2,1]));
        
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn,pol_indx]=max(vf,[],2);
        v_tol=[max(abs(vfn(:,:,1)' - v_guess(1,:))) ; max(abs(vfn(:,:,2)' - v_guess(2,:)))];
        v_guess=[vfn(:,:,1)';vfn(:,:,2)'];
    end;
    
    % KEEP DECSISION RULE
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_fn=a(pol_indx); 
    
    % SET UP INITITAL DISTRIBUTION
    Mu=ones(2,num_a)/(2*num_a);
    
    % ITERATE OVER DISTRIBUTIONS
     m_tol=1;
     while m_tol>0.0001
        [emp_ind, a_ind, mass] = find(Mu > 0); % find non-zero indices
    
        MuNew = zeros(size(Mu));
   
        
       for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); % which a prime does the policy fn prescribe?
        MuNew(:, apr_ind) = MuNew(:, apr_ind) +[PI(emp_ind(ii),1)*Mu(emp_ind(ii),a_ind(ii));PI(emp_ind(ii),2)*Mu(emp_ind(ii),a_ind(ii))];
        % which mass of households goes to which exogenous state?
       end
       m_tol=max(max(abs(MuNew-Mu)));
       Mu=MuNew;
     end 
    
    aggsav=Mu(1,:)*a'+Mu(2,:)*a';
    if aggsav>0;
         q_min=q_guess;
    else q_max=q_guess;
    end
    q_min;
    q_max;
end

figure;
plot(a,vfn(:,:,1),a,vfn(:,:,2)),legend('Employed','Unemployed');

figure;
plot(a, pol_fn(1,:),a,pol_fn(2,:)),legend('Employed','Unemployed');


%%%%%for lorenz curve and gini index
p=[Mu(2,:);Mu(1,:)];
income=[repmat(y_s(2),1,num_a);repmat(y_s(1),1,num_a)];
g=gini(p,income,true);
title(['income gini index=',num2str(g)]);
wealth=[bsxfun(@plus,a,y_s(2));bsxfun(@plus,a,y_s(1))];
wealth(wealth<0)=0;
wg=gini(p,wealth,true);
title(['wealth gini index=',num2str(wg)]);