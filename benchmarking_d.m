addpath 'hierclust2nmf_v2'
addpath 'NNSVD-LRC_v2'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testcase "D": square lowrank matrices with noise
% Sorry: this is the test case B in the paper :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(51193);                 % some birthday
recompute_results = false;  % use precomputed values or compute from scratch

% square matrices
matrix_dims = [10;20;30;40;50;60;70;80;90;100];

N = 1000;           % number of test cases
maxiter = 1000;      % max ANLS iterations
tol_pow = -3;    
tol = 10^tol_pow;   % convergence criterion


if recompute_results

% first dim : SPA, NNSVDLRC, NNDSVD, rand, QDR
% second dim: matrix sizes
% third dim : sampled data points
% fourth dim: init relative error, relative error, ANLS iters, ANLS convergence, init time, ANLS time


ALL_RESULTS = zeros(5,size(matrix_dims,1),N,6);


for k = 1:size(matrix_dims,1)
    

    m = matrix_dims(k); n = m;

    fprintf("Computing %d x %d matrices... ", m, n);
    compl = 10;

    for i = 1:N

        if mod(i,N/10) == 0, fprintf("%d %%... ", compl); compl =  compl + 10; end

        U = generate_data(m,n,true,"lowrank");
 
        
        tic; % QDR
        [L0,R0] = qdrinit(U);
        ALL_RESULTS(5,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(5,k,i,3),ALL_RESULTS(5,k,i,4)] = ANLS(U,L0,R0,maxiter,tol);
        ALL_RESULTS(5,k,i,6) = toc;
        ALL_RESULTS(5,k,i,1) = norm(L0*R0'-U,'fro')/norm(U,'fro');
        ALL_RESULTS(5,k,i,2) = norm(L*R'-U,'fro')/norm(U,'fro');

        tic; % SPA
        [L0,R0] = rank2nmf(U);
        ALL_RESULTS(1,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(1,k,i,3),ALL_RESULTS(1,k,i,4)] = ANLS(U,L0,R0',maxiter,tol);
        ALL_RESULTS(1,k,i,6) = toc;
        ALL_RESULTS(1,k,i,1) = norm(L0*R0-U,'fro')/norm(U,'fro');
        ALL_RESULTS(1,k,i,2) = norm(L*R'-U,'fro')/norm(U,'fro');
        

        tic; % NNSVDLRC
        [L0,R0] = NNSVDLRC(U,2);
        ALL_RESULTS(2,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(2,k,i,3),ALL_RESULTS(2,k,i,4)] = ANLS(U,L0,R0',maxiter,tol);
        ALL_RESULTS(2,k,i,6) = toc;
        ALL_RESULTS(2,k,i,1) = norm(L0*R0-U,'fro')/norm(U,'fro');
        ALL_RESULTS(2,k,i,2) = norm(L*R'-U,'fro')/norm(U,'fro');
        

        tic; % NNDSVD
        [L0,R0] = NNDSVD(U,2,0);
        ALL_RESULTS(3,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(3,k,i,3),ALL_RESULTS(3,k,i,4)] = ANLS(U,L0,R0',maxiter,tol);
        ALL_RESULTS(3,k,i,6) = toc;
        ALL_RESULTS(3,k,i,1) = norm(L0*R0-U,'fro')/norm(U,'fro');
        ALL_RESULTS(3,k,i,2) = norm(L*R'-U,'fro')/norm(U,'fro');
        

        tic; % random initialization
        L0 = rand(m,2); R0 = rand(n,2);
        ALL_RESULTS(4,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(4,k,i,3),ALL_RESULTS(4,k,i,4)] = ANLS(U,L0,R0,maxiter,tol);
        ALL_RESULTS(4,k,i,6) = toc;
        ALL_RESULTS(4,k,i,1) = norm(L0*R0'-U,'fro')/norm(U,'fro');
        ALL_RESULTS(4,k,i,2) = norm(L*R'-U,'fro')/norm(U,'fro');
        


    end

    fprintf("\n");

end


str_exp = sprintf("nxn_tol_1e%d_iters%d_case_d",tol_pow,maxiter);

writematrix(ALL_RESULTS,"RESULTS/" + str_exp + ".txt");

else 

str_exp = sprintf("nxn_tol_1e%d_iters%d_case_d",tol_pow,maxiter);

ALL_RESULTS = reshape(readmatrix("RESULTS/" + str_exp + ".txt"), 5, size(matrix_dims,1), N, 6);

end



my_syms = ["--";":";"-.";"-*";"-"];
my_colors = ["#CC0000";"#FFCC00";"#663399";"#009900";"#1656DD"];

dists_scaled_init = ALL_RESULTS(:,:,:,1) ./ min(ALL_RESULTS(:,:,:,1),[],1);
dists_scaled = ALL_RESULTS(:,:,:,2) ./ min(ALL_RESULTS(:,:,:,2),[],1);




% plotting....
f = figure(1); clf;
f.Position = [10 10 1800 600];
subplot(1,2,1); hold on;
for k = 1:5
    plot(matrix_dims, reshape(median(ALL_RESULTS(k,:,:,5) + ALL_RESULTS(k,:,:,6),3),size(matrix_dims)), 'LineWidth',3);
end
linestyleorder(my_syms);
colororder(my_colors);

legend('SPA','NNSVDLRC','NNDSVD','rand','QDR','Location','northwest');
xlabel("n"); ylabel("time (s)");

xlim([min(matrix_dims) max(matrix_dims)]);

fontsize(18,"points");
title("Median time performance");


subplot(1,2,2); hold on;
for k = [1,2,3,5]
    plot(matrix_dims, max(dists_scaled(k,:,:),[],3), 'LineStyle', my_syms(k), 'color', my_colors(k), 'LineWidth',3);
end

legend('SPA','NNSVDLRC','NNDSVD','QDR','Location','northeast');
xlabel("n");
xlim([min(matrix_dims) max(matrix_dims)]);
fontsize(18,"points");
title("Worst ratio between distances");

exportgraphics(f,"RESULTS/" + str_exp + ".png");
savefig("RESULTS/" + str_exp);


fprintf("\n\nSummary of the results\n\n");


for i = 1:size(matrix_dims,1)
    fprintf("n = %d & SPA & NNSVDLRC & NNDSVD & rand & QDR \\\\\n", matrix_dims(i));
    fprintf("mean time & %s", join(string(round(mean(ALL_RESULTS(:,i,:,5)+ALL_RESULTS(:,i,:,6),3),3))," & "));
    fprintf("\\\\\n");
    fprintf("max time & %s", join(string(round(max(ALL_RESULTS(:,i,:,5)+ALL_RESULTS(:,i,:,6),[],3),3))," & "));
    fprintf("\\\\\n");
    fprintf("mean acc & %s", join(string(round(mean(dists_scaled(:,i,:),3),3))," & "));
    fprintf("\\\\\n");
    fprintf("max acc & %s", join(string(round(max(dists_scaled(:,i,:),[],3),3))," & "));
    fprintf("\\\\\n");
    fprintf("mean acc init & %s", join(string(round(mean(dists_scaled_init(:,i,:),3),3))," & "));
    fprintf("\\\\\n");
    fprintf("max acc init & %s", join(string(round(max(dists_scaled_init(:,i,:),[],3),3))," & "));
    fprintf("\\\\\n");
    fprintf("\n\n\n");
end