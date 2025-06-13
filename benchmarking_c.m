addpath 'hierclust2nmf_v2'
addpath 'NNSVD-LRC_v2'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testcase B: tall log normal iid matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(51193);                 % some birthday
recompute_results = true;  % use precomputed values or compute from scratch



if recompute_results

% first dim : SPA, NNSVDLRC, NNDSVD, rand, QDR
% second dim: matrix sizes
% third dim : sampled data points
% fourth dim: init relative error, relative error, ANLS iters, ANLS convergence, init time, ANLS time

load datamatrices.mat;

m = 4; n = 4;
N = size(originalmatrices,1);

maxiter = 100;
tol_pow = -5;
tol = 10^tol_pow;
    


ALL_RESULTS = zeros(5,N,8);

compl = 25;

for i = 1:N

    if mod(i,round(N/4)) == 0, fprintf("%d %%... ", compl); compl =  compl + 25; end

    U = reshape(originalmatrices(i,:),m,n);
    Ubest = reshape(solutions(i,:),m,n);
    
    tic; % QDR
    [L0,R0] = qdrinit(U);
    ALL_RESULTS(5,i,5) = toc;
    tic; 
    [L,R,ALL_RESULTS(5,i,3),ALL_RESULTS(5,i,4)] = ANLS(U,L0,R0,maxiter,tol);
    ALL_RESULTS(5,i,6) = toc;
    ALL_RESULTS(5,i,1) = norm(L0*R0'-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(5,i,2) = norm(L*R'-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(5,i,7) = norm(L0*R0'-U,'fro')/norm(U,'fro');
    ALL_RESULTS(5,i,8) = norm(L*R'-U,'fro')/norm(U,'fro');

    tic; % SPA
    [L0,R0] = rank2nmf(U);
    ALL_RESULTS(1,i,5) = toc;
    tic; 
    [L,R,ALL_RESULTS(1,i,3),ALL_RESULTS(1,i,4)] = ANLS(U,L0,R0',maxiter,tol);
    ALL_RESULTS(1,i,6) = toc;
    ALL_RESULTS(1,i,1) = norm(L0*R0-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(1,i,2) = norm(L*R'-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(1,i,7) = norm(L0*R0-U,'fro')/norm(U,'fro');
    ALL_RESULTS(1,i,8) = norm(L*R'-U,'fro')/norm(U,'fro');

    tic; % NNSVDLRC
    [L0,R0] = NNSVDLRC(U,2);
    ALL_RESULTS(2,i,5) = toc;
    tic; 
    [L,R,ALL_RESULTS(2,i,3),ALL_RESULTS(2,i,4)] = ANLS(U,L0,R0',maxiter,tol);
    ALL_RESULTS(2,i,6) = toc;
    ALL_RESULTS(2,i,1) = norm(L0*R0-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(2,i,2) = norm(L*R'-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(2,i,7) = norm(L0*R0-U,'fro')/norm(U,'fro');
    ALL_RESULTS(2,i,8) = norm(L*R'-U,'fro')/norm(U,'fro');
    

    tic; % NNDSVD
    [L0,R0] = NNDSVD(U,2,0);
    ALL_RESULTS(3,i,5) = toc;
    tic; 
    [L,R,ALL_RESULTS(3,i,3),ALL_RESULTS(3,i,4)] = ANLS(U,L0,R0',maxiter,tol);
    ALL_RESULTS(3,i,6) = toc;
    ALL_RESULTS(3,i,1) = norm(L0*R0-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(3,i,2) = norm(L*R'-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(3,i,7) = norm(L0*R0-U,'fro')/norm(U,'fro');
    ALL_RESULTS(3,i,8) = norm(L*R'-U,'fro')/norm(U,'fro');
    

    tic; % random initialization
    L0 = rand(m,2); R0 = rand(n,2);
    ALL_RESULTS(4,i,5) = toc;
    tic; 
    [L,R,ALL_RESULTS(4,i,3),ALL_RESULTS(4,i,4)] = ANLS(U,L0,R0,maxiter,tol);
    ALL_RESULTS(4,i,6) = toc;
    ALL_RESULTS(4,i,1) = norm(L0*R0'-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(4,i,2) = norm(L*R'-Ubest,'fro')/norm(Ubest,'fro');
    ALL_RESULTS(4,i,7) = norm(L0*R0'-U,'fro')/norm(U,'fro');
    ALL_RESULTS(4,i,8) = norm(L*R'-U,'fro')/norm(U,'fro');

end


str_exp = sprintf("4x4_tol_1e%d_iters%d_case_c",tol_pow,maxiter);

writematrix(ALL_RESULTS,"RESULTS/" + str_exp + ".txt");

else 

str_exp = sprintf("4x4_tol_1e%d_iters%d_case_c",tol_pow,maxiter);

ALL_RESULTS = reshape(readmatrix("RESULTS/" + str_exp + ".txt"), 5, N, 8);

end



my_syms = ["--";":";"-.";"-*";"-"];
my_colors = ["#CC0000";"#FFCC00";"#663399";"#009900";"#1656DD"];



matrix_dims = matrix_dims(:,1);
my_colors = ["#CC0000";"#FFCC00";"#663399";"#009900";"#1656DD"];

% TIME PERFORMANCE

% plotting....
f = figure(1); clf;
f.Position = [10 10 900 1200];


res_count = 50;

ylabels = ["SPA", "NNSVDLRC", "NNDSVD","rand","QDR"];
marker_char = ['*','-'];

for j = 1:5
    [~,I] = sort(ALL_RESULTS(j,:,1),"descend");
    subplot(5,1,j); hold on;
    for k = 1:5
        plot(1:res_count, ALL_RESULTS(k,I(1:res_count),1),marker_char((j==k)+1),'LineWidth',3);
    end
    ylabel(ylabels(j));
    xlim([1 res_count]);
    colororder(my_colors);
end

legend('SPA','NNSVDLRC','NNDSVD','rand','QDR','Location','eastoutside');
xlabel("index");

sgtitle("Worst-case scaled distance between N^* and init");
fontsize(18,"points");
exportgraphics(f,"RESULTS/" + str_exp + "_accinit.png");
savefig("RESULTS/" + str_exp + "_accinit");


% plotting....
f = figure(2); clf;
f.Position = [10 10 900 1200];


res_count = 50;

ylabels = ["SPA", "NNSVDLRC", "NNDSVD","rand","QDR"];
marker_char = ['*','-'];

for j = 1:5
    [~,I] = sort(ALL_RESULTS(j,:,2),"descend");
    subplot(5,1,j); hold on;
    for k = 1:5
        plot(1:res_count, ALL_RESULTS(k,I(1:res_count),2),marker_char((j==k)+1),'LineWidth',3);
    end
    ylabel(ylabels(j));
    xlim([1 res_count]);
    colororder(my_colors);
end

legend('SPA','NNSVDLRC','NNDSVD','rand','QDR','Location','eastoutside');
xlabel("index");

sgtitle("Worst-case scaled distance to N^*");
fontsize(18,"points");
exportgraphics(f,"RESULTS/" + str_exp + "_acc.png");
savefig("RESULTS/" + str_exp + "_acc");



titles = ["SPA","NNSVDLRC","NNDSVD","rand","QDR"];
bounds = [inf, 1e-1, 1e-5, -inf];
for j = 1:5
fprintf("%s:\n",titles(j));
    for i = 1:3
        kk = sum((ALL_RESULTS(j,:,2) > bounds(i+1)) .* (ALL_RESULTS(j,:,2) <= bounds(i))); 
        fprintf("in [%.2f,%.2f) : %d kpl\n", bounds(i),bounds(i+1), kk); 
    end
end
    
