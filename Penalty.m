function [v,xbest,sol_MC,constr_viol] = Penalty(InputType, InputData, e1, max_time, tolerance)
tic;


%Parameters for computing eigenvector
%Set maximum number of outer iterations for computing eigenvalue
maxit = 300;
%tolerance = 0.1;
opts.isreal = 1;
opts.issym = 1;
opts.tol = tolerance;
opts.maxit = maxit;
%opts.p = 3;

%rng(10);

%%%%%%%%%
%Read input data or create random graph depending on the input
%%%%%%%%%
if find(InputType == 'R')
    disp('Creating random graph');
    %%%%%%%%%
    %Create a random graph
    %%%%%%%%%
    %Initial density of graph
    V = InputData;
    Edges = V*1.5;
    %Create random graph
    idx = randi(V,Edges,2);
    %Create edges to make the graph connected
    s = floor(rand(1,V-1).*(1:V-1))+1;
    t = 2:V;
    idx1 = [s',t'];
    %idx2 - matrix of all the edges of random connected graph
    idx2 = [idx;idx1];
    clear s
    clear t
    clear idx1
    clear idx

    %Create a symmetric matrix
    idx2 = [idx2; [idx2(:,2),idx2(:,1)]];
    %Delete repeated edges
    idx2 = unique(idx2,'rows');
    %Delete self-loops
    idx2(idx2(:,1)==idx2(:,2),:) = [];

    %Find degree of every node
    [a,b] = hist(idx2(:,1),unique(idx2(:,1)));
    val = [-1*ones(length(idx2),1);a'];
    idx2 = [idx2;[b,b]];
    clear a
    clear b

    %Create a sparse matrx
    C = sparse(idx2(:,1),idx2(:,2), val/4,V,V);
    
elseif find(InputType == 'G')
    disp('Creating graph from GSet');
    filename = load(InputData);
    C = filename.Problem.A;
    V = length(C);
    C = (spdiags(C*ones(V,1),0,V,V)-C)/4;
       
end

%%%%%%%%%
%Input parameters
%%%%%%%%%
alpha = V;
nsamples = 1;
mu = 1;
c = 4;
M = (mu*log(2*V))/e1;
beta = c*trace(C);
Cf = beta*M*alpha^2;
eta = 10^-4;
lambdamaxC = eigs(C,1,'LR');

%set initial point and create samples
%Create matrix of random samples
%Assume initial feasible point to be identity
disp('Creating initial samples');
x = 1:V;
for i = 1:nsamples
    z = randn(V,1);
    if i == 1
        Z = z;
    else
        Z = [Z z];
    end
end

%Create linear mapping of objective function [<C,X>, A(X)]
disp('Creating linear mapping');
t = 1;
gamma = 2/(t+2);
v = ones(V+1,1);
v(1) = trace(C);

%First iteration
%Compute the gradient of obj function
max_viol = max(abs(M.*(v(2:V+1)-1)));
penalty = cat(1,exp(M.*(v(2:V+1)-1)-max_viol),exp(M.*(1-v(2:V+1))-max_viol));
matrix_entry = (beta/sum(penalty))*(penalty(1:V)-penalty(V+1:2*V));
try
    opts.tol = (1/4)*((eta*gamma*Cf)/(alpha*lambdamaxC));
    if opts.tol >= 1
	    opts.tol = tolerance;
    end
    [u,l] = eigs(-sparse(x,x,matrix_entry)+C,1,'LR',opts);
catch
    [l,u] = eigifp(sparse(x,x,matrix_entry)-C,1);%,opt);
    l = -l;
end

%Compute update direction
h = zeros(V+1,1);
if l >= 0
    h(1) = alpha*u'*C*u;
    h(2:V+1) = alpha*u.^2;
end

%Initial constraint violation
constr_viol = abs((1-gamma)*v(2:V+1) + gamma*h(2:V+1) - ones(V,1));
%vecMaxPenalty = [];
%vecDualityGap = [];

%Run the loop until error is less than \epsilon
while ( ((h-v)'*[1;-matrix_entry] > e1*trace(C) ) && toc <= max_time)
    
    %Check the violation every few iterations
    if mod(t,1000) == 0
        disp(t);
        disp((h-v)'*[1;-matrix_entry]);
        disp(max_viol/M);
    end
    
    %vecDualityGap = [vecDualityGap; (h-v)'*[1;-matrix_entry]];
    
        
    %Generate samples from gradient and update the samples
    for i = 1:nsamples
        x1 = normrnd(0,1);
        w = u*x1;
        if l >= 0
            Z(:,i) = sqrt(1-gamma)*Z(:,i) + sqrt(gamma)*w;
        else
            Z(:,i) = sqrt(1-gamma)*Z(:,i);
        end
    end
    %Conditional gradient step
    v = (1-gamma)*v + gamma*h;
    t = t+1;
    gamma = 2/(t+2);
    

    %Compute constraint violation and change the value of M
    %This is done when the parameter value is too large
    Mfactor = 1.1;
    
    %Compute gradient at the new point
    max_viol = max(abs(M.*(v(2:V+1)-1)));
    penalty = cat(1,exp(M.*(v(2:V+1)-1)-max_viol),exp(M.*(1-v(2:V+1))-max_viol));
    matrix_entry = (beta/sum(penalty))*(penalty(1:V)-penalty(V+1:2*V));
    while find(matrix_entry == Inf) | find(matrix_entry == -Inf)
        M = M/Mfactor;
        max_viol = max(abs(M.*(v(2:V+1)-1)));
        penalty = cat(1,exp(M.*(v(2:V+1)-1)-max_viol),exp(M.*(1-v(2:V+1))-max_viol));
        matrix_entry = (beta/sum(penalty))*(penalty(1:V)-penalty(V+1:2*V));
    end
    while find(isnan(matrix_entry))
        M = M/Mfactor;
        max_viol = max(abs(M.*(v(2:V+1)-1)));
        penalty = cat(1,exp(M.*(v(2:V+1)-1)-max_viol),exp(M.*(1-v(2:V+1))-max_viol));
        matrix_entry = (beta/sum(penalty))*(penalty(1:V)-penalty(V+1:2*V));
    end
   
    %Compute eigenvector using eigs
    try
        opts.tol = (1/4)*((eta*gamma*Cf)/(alpha*lambdamaxC));
        if opts.tol >= 1
            opts.tol = 0.1;
        end
        %vecMaxPenalty = [vecMaxPenalty; max_viol];
            
        [u,l,f] = eigs(-sparse(x,x,matrix_entry)+C,1,'LR',opts);
    catch
        [l,u] =  eigifp(sparse(x,x,matrix_entry)-C,1);%,opt);
        l = -l;
    end
    %Compute update direction
    h = zeros(V+1,1);
    if l >= 0
        h(1) = alpha*u'*C*u;
        h(2:V+1) = alpha*u.^2;
    end
    
end

%Compute feasible solution with best value
m = -10000;
for i = 1:nsamples
    x = sign(Z(:,i));
    tmp = x'*C*x;
    if tmp > m
        m = tmp;
        xbest = x;
    end
end

toc
time = toc;

%Display output
sa = v(1) - beta*(1/sum(M))*log(sum(M.*(penalty(1:V)+penalty(V+1:2*V))));
sp = v(1);
s = 2*trace(C);
constr_viol = sum(abs(v(2:length(v)) -ones(V,1)))/V;
max_constr_viol = norm(v(2:length(v)) -ones(V,1),inf);

sol_MC = xbest'*C*xbest;

disp('Number of iterations:');
disp(t);
disp('Objective function value of augmented primal');
disp(sa);

disp('Average violation of constraints:');
disp(constr_viol);

disp('Maximum violation of constraints:');
disp(max_constr_viol);

disp('Solution of max cut problem:');
disp(xbest'*C*xbest);

%Write output to file
if InputType == 'R'
    InputData = int2str(InputData);
end

% fileID = fopen('Output/LogPenaltyOutput.txt','a');
% fprintf(fileID,'\n');
% fprintf(fileID,'%s&',datestr(now));
% fprintf(fileID,'%s&',InputType);
% fprintf(fileID,'%s&', InputData);
% %outd = [e1,tolerance];%, InputType, 
% %outt = '%.2f&%.2f&';
% %fprintf(fileID,outt,outd);
% outd = [sol_MC, t, time, max_constr_viol, constr_viol, sp];
% outt = '%d&%d&%.2f&%.3f&%.3f&%.3f';
% fprintf(fileID,outt,outd);
% fprintf(fileID, '\n');
% fclose(fileID);

% p = figure;
% %y = linspace(1,length(vecMaxPenalty),length(vecMaxPenalty));
% y = 1:numel(vecMaxPenalty);
% loglog(y(1:100:end),vecMaxPenalty(1:100:end)/M,'-o')
% title('Change in feasibility with iterations');
% ylabel('Max constraint violation');
% xlabel('Iteration number t');
% if find(InputType == 'R')
%    fname = strcat('Output/ConstrViol',int2str(V),'.png');
% else
%    fname = strcat('Output/ConstrViol',filename(end-5:end-4),'.png');
% end
% saveas(p,fname);
% 
% p = figure;
% y = 1:numel(vecDualityGap);
% loglog(y(1:100:end),vecDualityGap(1:100:end),'-o')
% title('Change in duality gap with iterations');
% ylabel('Duality gap');
% xlabel('Iteration number t');
% if find(InputType == 'R')
%    fname = strcat('Output/DualityGap',int2str(V),'.png');
% else
%    fname = strcat('Output/DualityGap',filename(end-5:end-4),'.png');
% end
% saveas(p,fname);

end
