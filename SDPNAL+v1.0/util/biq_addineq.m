%%***********************************************************************
%% biq: generate SDP data for binary quadratic programming
%%      with additional valid inequality constraints
%% X = [Y,y; y', 1];
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

   function [blk,At,A,C,b,Bt,Bstruct,l,u] = biq_addineq(Q,options)

   if (nargin == 1); options = 0; end;
   Bt{1} = [];
   l = []; u = [];
   Bstruct.Bt = Bt;
   Bstruct.Bmap = @(X) []; 
   if isempty(Q)
      blk = []; At = []; C = []; b = []; 
      A =[]; Bt =[]; d = []; Bstruct =[];
      return;
   end
   if norm(Q-Q','fro') 
      Q = 0.5*(Q + Q'); 
   end
   
   n = size(Q,1);   
   blk{1,1} = 's'; blk{1,2} = n+1;   
%%
%% 
%%
   C{1} = [Q,sparse(n,1); sparse(1,n+1)];  
%%
%% diag(Y) = y; 
%% 
   n1 = n+1;   
   Acell = cell(1,n1); 
   for k = 1:n
      Acell{k} = spconvert([k,k,1;k,n1,-0.5; n1,k,-0.5;n1,n1,0]); 
   end
   Acell{n1} = spconvert([n1,n1,1]); 
   At = svec(blk,Acell,1); 
   b = [zeros(n,1); 1]; 
   A.At = At;
   A.Amap = @(X) diag(X{1}) - [X{1}(1:end-1,end);0];
   A.ATmap = @(y) ATmap(y);
%%
%% add the constraint <E,Y> >= (e'*y)^2
%%
   if (options == 1) %% Not useful. 
      e = ones(n,1); 
      tmp{1,1} = [sparse(n,n), 0.5*e; 0.5*e', 0];  
      tmp{1,2} = [ones(n,n), zeros(n,1); zeros(1,n), 0];  
      tmp{1,3} = sparse(n+1,n+1); 
      Atmp = svec(blk(1,:),tmp,1); 
      At{1} = [At{1}, Atmp{1}]; 
      b = [b; 0;0;1];
      C{2,1} = zeros(3,1); 
      blk{2,1} = 'q'; blk{2,2} = 3; 
      At{2,1} = [sparse(n1,3); 0,0,-1; -1,-1,0; 1,-1,0]'; 
   end
%%
%% add the constraint Y_{ij} <= y_i,Y_{ij} <= y_j, Y_{ij} >= y_i + y_j - 1 to form standard SDP
%%
   if (options == 2) 
      idx1 = [1:n];
      idx2 = [1:n];      
      len1 = length(idx1);
      len2 = length(idx2);
      m2 = len1*len2;
      Bcell = cell(1,m2); 
      BBcell = cell(1,m2);
      num = 0;
      for i = idx1
         for j = idx2
            if j ~= i
               Bcell{num+1} = spconvert([i,j,-0.5;j,i,-0.5;i,n1,0.5;n1,i,0.5;n1,n1,0]);  %% -Y_{ij} + y_i >= 0               
               BBcell{num+1} = spconvert([i,j,0.5;j,i,0.5;i,n1,-0.5;n1,i,-0.5;j,n1,-0.5;n1,j,-0.5;n1,n1,0]); %% Y_{ij} - y_i - y_j >= - 1               
               num = num+1;
            end
         end
      end
      m2 = num;
      Bcell = Bcell(1:m2);
      BBcell = BBcell(1:m2);
      Btmp1 = svec(blk,Bcell,1);  
      Btmp2 = svec(blk,BBcell,1); 
      Bt{1} = [Bt{1}, Btmp1{1}, Btmp2{1}]; 
      b = [b; zeros(m2,1); ones(m2,1)];
      C{2,1} = zeros(m2*2,1); 
      blk{2,1} = 'l'; blk{2,2} = m2*2; 
      Bt{2,1} = [sparse(m2*2,n1) speye(m2*2)]; 
   end
%%******************************************************
%% add the constraint -Y_{ij} + y_i >= 0, -Y_{ij} + y_j >= 0, Y_{ij} - y_i - y_j >= - 1
%%---- and Y_{ij}+Y{ik}+Y{jk}-y_i-y_j-y_k >= -1
%%
   if (options == 3) 
      idx1 = [1:1:n];
      idx2 = [1:1:n];
      idx3 = [];
      len1 = length(idx1);
      len2 = length(idx2);
      len3 = length(idx3);
      m2 = len1*len2;
      m3 = m2*len3;
      Bcell = cell(1,m2); 
      B1cell = cell(1,m2);
      BBcell = cell(1,m2);
      B3cell = cell(1,m3);
      num = 0;
      num1 = 0;
      for i = idx1
         for j = idx2
            if  i < j
               Bcell{num+1} = spconvert([i,j,-0.5;j,i,-0.5;i,n1,0.5;n1,i,0.5;n1,n1,0]);  %% -Y_{ij} + y_i >= 0
               B1cell{num+1} = spconvert([i,j,-0.5;j,i,-0.5;j,n1,0.5;n1,j,0.5;n1,n1,0]); %% -Y_{ij} + y_j >= 0
               BBcell{num+1} = spconvert([i,j,0.5;j,i,0.5;i,n1,-0.5;n1,i,-0.5;j,n1,-0.5;n1,j,-0.5;n1,n1,0]); %% Y_{ij} - y_i - y_j >= - 1
               num = num+1;
            end
            for k = idx3
               if j ~= i && k ~= i && k ~= j
                  B3cell{num1+1} = spconvert...
                     ([i,j,1;i,k,1;j,k,1;i,n1,-0.5;n1,i,-0.5;j,n1,-0.5;n1,j,-0.5;k,n1,-0.5;n1,k,-0.5;n1,n1,0]);
                  num1 = num1+1;
               end
            end
         end
      end
      m2 = num;
      Bcell = Bcell(1:m2);
      B1cell = B1cell(1:m2);
      BBcell = BBcell(1:m2);
      B3cell = B3cell(1:num1);
      Btmp1 = svec(blk,Bcell,1);  
      B1tmp1 = svec(blk,B1cell,1);
      Btmp2 = svec(blk,BBcell,1); 
      Btmp3 = svec(blk,B3cell,1);
      Bt{1} = [Btmp1{1}, B1tmp1{1}, Btmp2{1}]; 
      l = [zeros(m2*2,1); -ones(m2,1);];
      u = [ones(m2*2,1); zeros(m2,1);];
      Bstruct.Bt = Bt;
      Bstruct.Bmap = @(X) (mexsvec(blk,X{1})'*Bt{1,1})';
      Bstruct.BTmap = @(y) {mexsmat(blk,Bt{1,1}*y)};    
   end   
%%******************************************************
%% add the constraint -Y_{ij} + y_i >= 0, -Y_{ij} + y_j >= 0, Y_{ij} - y_i - y_j >= - 1
%%---- and Y_{ij}+Y{ik}+Y{jk}-y_i-y_j-y_k >= -1 to stardard form
%%
   if (options == 4) 
      idx1 = [1:1:n];
      idx2 = [1:1:n];
      idx3 = [];
      len1 = length(idx1);
      len2 = length(idx2);
      len3 = length(idx3);
      m2 = len1*len2;
      m3 = m2*len3;
      Bcell = cell(1,m2); 
      B1cell = cell(1,m2);
      BBcell = cell(1,m2);
      B3cell = cell(1,m3);
      num = 0;
      num1 = 0;
      for i = idx1
         for j = idx2
            if  i < j
               Bcell{num+1} = spconvert([i,j,-0.5;j,i,-0.5;i,n1,0.5;n1,i,0.5;n1,n1,0]);  %% -Y_{ij} + y_i >= 0
               B1cell{num+1} = spconvert([i,j,-0.5;j,i,-0.5;j,n1,0.5;n1,j,0.5;n1,n1,0]); %% -Y_{ij} + y_j >= 0
               BBcell{num+1} = spconvert([i,j,0.5;j,i,0.5;i,n1,-0.5;n1,i,-0.5;j,n1,-0.5;n1,j,-0.5;n1,n1,0]); %% Y_{ij} - y_i - y_j >= - 1
               num = num+1;
            end
            for k = idx3
               if j ~= i && k ~= i && k ~= j
                  B3cell{num1+1} = spconvert...
                     ([i,j,1;i,k,1;j,k,1;i,n1,-0.5;n1,i,-0.5;j,n1,-0.5;n1,j,-0.5;k,n1,-0.5;n1,k,-0.5;n1,n1,0]);
                  num1 = num1+1;
               end
            end
         end
      end
      m2 = num;
      Bcell = Bcell(1:m2);
      B1cell = B1cell(1:m2);
      BBcell = BBcell(1:m2);
      B3cell = B3cell(1:num1);
      Btmp1 = svec(blk,Bcell,1);  
      B1tmp1 = svec(blk,B1cell,1);
      Btmp2 = svec(blk,BBcell,1); 
      Btmp3 = svec(blk,B3cell,1);
      At{1} = [At{1}, Btmp1{1},B1tmp1{1}, Btmp2{1}]; 
      b = [b; zeros(m2*2,1); -ones(m2,1)];
      C{2,1} = zeros(m2*3,1); 
      blk{2,1} = 'l'; blk{2,2} = m2*3; 
      At{2,1} = [sparse(m2*3,n1), -speye(m2*3)];      
   end
%% add the constraint -Y_{ij} + y_i >= 0, -Y_{ij} + y_j >= 0, Y_{ij} - y_i - y_j >= - 1
%%---- and Y_{ij}+Y{ik}+Y{jk}-y_i-y_j-y_k >= -1
%%
   if (options == 5) 
      idx1 = [1:1:n];
      idx2 = [1:1:n];
      idx3 = [];
      len1 = length(idx1);
      len2 = length(idx2);
      len3 = length(idx3);
      m2 = len1*len2;
      m3 = m2*len3;
      Bcell = cell(1,m2); 
      B1cell = cell(1,m2);
      BBcell = cell(1,m2);
      B3cell = cell(1,m3);
      num = 0;
      num1 = 0;
      for i = idx1
         for j = idx2
            if  i < j
               Bcell{num+1} = spconvert([i,j,-0.5;j,i,-0.5;i,n1,0.5;n1,i,0.5;n1,n1,0]);  %% -Y_{ij} + y_i >= 0
               B1cell{num+1} = spconvert([i,j,-0.5;j,i,-0.5;j,n1,0.5;n1,j,0.5;n1,n1,0]); %% -Y_{ij} + y_j >= 0
               BBcell{num+1} = spconvert([i,j,0.5;j,i,0.5;i,n1,-0.5;n1,i,-0.5;j,n1,-0.5;n1,j,-0.5;n1,n1,0]); %% Y_{ij} - y_i - y_j >= - 1
               num = num+1;
            end
            for k = idx3
               if j ~= i && k ~= i && k ~= j
                  B3cell{num1+1} = spconvert...
                     ([i,j,1;i,k,1;j,k,1;i,n1,-0.5;n1,i,-0.5;j,n1,-0.5;n1,j,-0.5;k,n1,-0.5;n1,k,-0.5;n1,n1,0]);
                  num1 = num1+1;
               end
            end
         end
      end
      m2 = num;
      Bcell = Bcell(1:m2);
      B1cell = B1cell(1:m2);
      BBcell = BBcell(1:m2);
      B3cell = B3cell(1:num1);
      Btmp1 = svec(blk,Bcell,1);  
      B1tmp1 = svec(blk,B1cell,1);
      Btmp2 = svec(blk,BBcell,1); 
      Btmp3 = svec(blk,B3cell,1);
      Bt{1} = [At{1}, Btmp1{1}, B1tmp1{1}, Btmp2{1}]; 
      l = [b;zeros(m2*2,1); -ones(m2,1);];
      u = [b;ones(m2*2,1); zeros(m2,1);];
      At = {[]}; 
      b = [];
      Bstruct.Bt = Bt;
      Bstruct.Bmap = @(X) (mexsvec(blk,X{1})'*Bt{1,1})';
      Bstruct.BTmap = @(y) {mexsmat(blk,Bt{1,1}*y)};    
   end
%%*********************************************************************   
%%*********************************************************************   
      
   function Aty = ATmap(y)
        m = length(y);
        Aty = sparse([1:m],[1:m],y,m,m);
        Aty(m,1:m-1) = -0.5*y(1:m-1)';
        Aty(1:m-1,m) = -0.5*y(1:m-1);
        %Aty = (Aty + Aty')*0.5;
        Aty = {Aty};
%%*********************************************************************   
   


