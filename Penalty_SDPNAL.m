function [fn_val,sol_MC] = Penalty_SDPNAL(InputType, InputData, e1)
tic;


%rng(10);

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

n = V;

model = ccp_model('PenaltyMC1');
X = var_sdp(n, n);
model.add_variable(X);
model.maximize(inprod(C,X));
model.add_affine_constraint(map_diag(X) == ones(n,1));
model.setparameter ('tol', e1);
model.solve;

sol_SDP = get_value(X);
fn_val = trace(C*sol_SDP);

nsamples = 10;
R = mvnrnd(zeros(n,1),sol_SDP+1e-1*eye(n),nsamples);
m = -10000;
for i = 1:nsamples
    x = sign(R(i,:));
    tmp = x*C*x';
    if tmp > m
        m = tmp;
        xbest = x;
    end
end

sol_MC = xbest*C*xbest';

toc
time = toc;
disp(sol_MC);

%Write output to file
if InputType == 'R'
    InputData = int2str(InputData);
end
fileID = fopen('Output/Output_SDPNAL.txt','a');
fprintf(fileID,'\n');
fprintf(fileID,'%s&',datestr(now));
fprintf(fileID,'%s&',InputType);
fprintf(fileID,'%s&', InputData);
%outd = [e1,tolerance];%, InputType, 
%outt = '%.2f&%.2f&';
%fprintf(fileID,outt,outd);
outd = [sol_MC, time fn_val];
outt = '%d&%.2f&%.3f';
fprintf(fileID,outt,outd);
fprintf(fileID, '\n');
fclose(fileID);

end
