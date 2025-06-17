addpath 'hierclust2nmf_v2'
addpath 'NNSVD-LRC_v2'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testcase A: square log normal iid matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(51193);                 % some birthday
recompute_results = false;  % use precomputed values or compute from scratch

% square matrices
matrix_dims = [100;200;300;400;500;600;700;800;900;1000;1100;1200];

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
norms_U = zeros(N,k);

for k = 1:size(matrix_dims,1)
    

    n = matrix_dims(k);

    fprintf("Computing %d x %d matrices... ", n, n);
    compl = 25;

    for i = 1:N

        if mod(i,N/4) == 0, fprintf("%d %%... ", compl); compl =  compl + 25; end

        U = generate_data(n,n,true);
        norms_U(i,k) = norm(U,'fro');
 
        
        tic; % QDR
        [L0,R0] = qdrinit(U);
        ALL_RESULTS(5,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(5,k,i,3),ALL_RESULTS(5,k,i,4)] = ANLS(U,L0,R0,maxiter,tol);
        ALL_RESULTS(5,k,i,6) = toc;
        ALL_RESULTS(5,k,i,1) = norm(L0*R0'-U,'fro');
        ALL_RESULTS(5,k,i,2) = norm(L*R'-U,'fro');


        tic; % SPA
        [L0,R0] = rank2nmf(U);
        ALL_RESULTS(1,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(1,k,i,3),ALL_RESULTS(1,k,i,4)] = ANLS(U,L0,R0',maxiter,tol);
        ALL_RESULTS(1,k,i,6) = toc;
        ALL_RESULTS(1,k,i,1) = norm(L0*R0-U,'fro');
        ALL_RESULTS(1,k,i,2) = norm(L*R'-U,'fro');


        tic; % NNSVDLRC
        [L0,R0] = NNSVDLRC(U,2);
        ALL_RESULTS(2,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(2,k,i,3),ALL_RESULTS(2,k,i,4)] = ANLS(U,L0,R0',maxiter,tol);
        ALL_RESULTS(2,k,i,6) = toc;
        ALL_RESULTS(2,k,i,1) = norm(L0*R0-U,'fro');
        ALL_RESULTS(2,k,i,2) = norm(L*R'-U,'fro');


        tic; % NNDSVD
        [L0,R0] = NNDSVD(U,2,0);
        ALL_RESULTS(3,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(3,k,i,3),ALL_RESULTS(3,k,i,4)] = ANLS(U,L0,R0',maxiter,tol);
        ALL_RESULTS(3,k,i,6) = toc;
        ALL_RESULTS(3,k,i,1) = norm(L0*R0-U,'fro');
        ALL_RESULTS(3,k,i,2) = norm(L*R'-U,'fro');


        tic; % random initialization
        L0 = rand(n,2); R0 = rand(n,2);
        ALL_RESULTS(4,k,i,5) = toc;
        tic; 
        [L,R,ALL_RESULTS(4,k,i,3),ALL_RESULTS(4,k,i,4)] = ANLS(U,L0,R0,maxiter,tol);
        ALL_RESULTS(4,k,i,6) = toc;
        ALL_RESULTS(4,k,i,1) = norm(L0*R0'-U,'fro');
        ALL_RESULTS(4,k,i,2) = norm(L*R'-U,'fro');


    end

    fprintf("\n");

end


str_exp = sprintf("nxn_tol_1e%d_iters%d_case_a",tol_pow,maxiter);

writematrix(ALL_RESULTS,"RESULTS/" + str_exp + ".txt");
writematrix(norms_U,"RESULTS/"+str_exp+ "_norms.txt");

else 

str_exp = sprintf("nxn_tol_1e%d_iters%d_case_a",tol_pow,maxiter);

ALL_RESULTS = reshape(readmatrix("RESULTS/" + str_exp + ".txt"), 5, size(matrix_dims,1), N, 6);
norms_U = readmatrix("RESULTS/"+str_exp+ "_norms.txt");

end



my_syms = ["--";":";"-.";"-*";"-"];
my_colors = ["#CC0000";"#FFCC00";"#663399";"#009900";"#1656DD"];





% TIME PERFORMANCE


f = figure(1); clf;
f.Position = [10 10 900 800];

subplot(2,1,1);
hold on;
for i  = 1:5
    plot(matrix_dims, reshape(median(ALL_RESULTS(i,:,:,5)+ALL_RESULTS(i,:,:,6),3),size(matrix_dims)), 'LineWidth',3);
end
linestyleorder(my_syms);
colororder(my_colors);

legend('SPA','NNSVDLRC','NNDSVD','rand','QDR','Location','northwest');


ylabel("time (s)");
xlim([min(matrix_dims) max(matrix_dims)]); ylim([-0.01, 0.23]);
fontsize(20,"points");
title("Median time performance for nxn matrices");



subplot(2,1,2); hold on;

for i  = 1:5
    plot(matrix_dims, reshape(median(ALL_RESULTS(i,:,:,5)+ALL_RESULTS(i,:,:,6),3),size(matrix_dims)), 'LineWidth',3);
end
linestyleorder(my_syms);
colororder(my_colors);


xlabel("n"); ylabel("time (s)");
xlim([min(matrix_dims) max(matrix_dims)]);
ylim([0 0.03]);
fontsize(20,"points");
title("Zoom in");




exportgraphics(f,"RESULTS/" + str_exp + "_time.png");
savefig("RESULTS/" + str_exp + "_time");





% INITIALIZATION ACCURACY

% scale distances

dists_scaled_init = ALL_RESULTS(:,:,:,1) ./ min(ALL_RESULTS(:,:,:,1),[],1);
dists_scaled = ALL_RESULTS(:,:,:,2) ./ min(ALL_RESULTS(:,:,:,2),[],1);




f = figure(2); clf;
f.Position = [10 10 1200 600];


subplot(2,2,1); hold on;
for i = [1,3,5]
    plot(matrix_dims, mean(dists_scaled_init(i,:,:),3), 'LineStyle', my_syms(i), 'color', my_colors(i), 'LineWidth',3);
end

xlim([min(matrix_dims) max(matrix_dims)]);
fontsize(20,"points");
title("Mean relative errors of initializations");

subplot(2,2,2); hold on;
for i  = [1,5]
    plot(matrix_dims, mean(dists_scaled_init(i,:,:),3), 'LineStyle', my_syms(i), 'color', my_colors(i), 'LineWidth',3);
end


xlim([min(matrix_dims) max(matrix_dims)]);
fontsize(20,"points");
title("Zoom in");

subplot(2,2,3); hold on;
for i = [1,3,5]
    plot(matrix_dims, max(dists_scaled_init(i,:,:),[],3), 'LineStyle', my_syms(i), 'color', my_colors(i), 'LineWidth',3);
end


xlabel("n");
xlim([min(matrix_dims) max(matrix_dims)]);
fontsize(20,"points");
title("Max relative errors of initializations");

subplot(2,2,4); hold on;
for i  = [1,5]
    plot(matrix_dims, max(dists_scaled_init(i,:,:),[],3), 'LineStyle', my_syms(i), 'color', my_colors(i), 'LineWidth',3);
end

xlabel("n");
xlim([min(matrix_dims) max(matrix_dims)]);
ylim([0.999,1.005]);
fontsize(20,"points");
title("Zoom in");

exportgraphics(f,"RESULTS/" + str_exp + "_acc_init.png");
savefig("RESULTS/" + str_exp + "_acc_init");




% ACCURACY AFTER ANLS

f = figure(3); clf;
f.Position = [1000 10 1000 600];

subplot(2,2,1); hold on;
for i  = [1,3,5]
    scatter(matrix_dims, mean(dists_scaled(i,:,:),3), 'filled', 'MarkerFaceColor', my_colors(i));
end
xlim([min(matrix_dims) max(matrix_dims)]); 
fontsize(20,"points");
title("Median relative errors");


subplot(2,2,2);  hold on;
for i  = [1,5]
    scatter(matrix_dims, mean(dists_scaled(i,:,:),3), 'filled','MarkerFaceColor',  my_colors(i));
end
xlim([min(matrix_dims) max(matrix_dims)]); 
fontsize(20,"points");
title("Zoom in");

subplot(2,2,3); hold on;
for i  = [1,3,5]
    scatter(matrix_dims, max(dists_scaled(i,:,:),[],3), 'filled', 'MarkerFaceColor', my_colors(i));
end
xlabel("n");
xlim([min(matrix_dims) max(matrix_dims)]); 
fontsize(20,"points");
title("Maximum relative errors");


subplot(2,2,4);  hold on;
for i  = [1,5]
    scatter(matrix_dims, max(dists_scaled(i,:,:),[],3), 'filled','MarkerFaceColor', my_colors(i));
end

xlabel("n");
xlim([min(matrix_dims) max(matrix_dims)]); 
fontsize(20,"points");
title("Zoom in");


exportgraphics(f,"RESULTS/" + str_exp + "_acc.png");
savefig("RESULTS/" + str_exp + "_acc");





% FREQUENCY OF BEING FASTEST


f = figure(4); clf;
f.Position = [10 10 900 900];

I = zeros(size(matrix_dims,1),N);
for i = 1:size(matrix_dims,1)
    [~,II] = min(ALL_RESULTS(:,i,:,5)+ALL_RESULTS(:,i,:,6),[],1);
    I(i,:) = reshape(II,1,N);
end

hold on;
for i  = 1:5
    plot(matrix_dims, sum(I==i,2)/(N/100), 'LineWidth',3);
end
linestyleorder(my_syms);
colororder(my_colors);

legend('SPA','NNSVDLRC','NNDSVD','rand','QDR','Location','east');

xlabel("n"); ylabel("%");
xlim([min(matrix_dims) max(matrix_dims)]);
fontsize(20,"points");
title("Frequency of being the fastest algorithm");

exportgraphics(f,"RESULTS/" + str_exp + "_fastest.png");
savefig("RESULTS/" + str_exp + "_fastest");




% ITERATIONS

f = figure(5); clf;
f.Position = [10, 10, 900, 900];
subplot(3,1,1); hold on; xlim([min(matrix_dims) max(matrix_dims)]); title("Median number of iterations");
subplot(3,1,2); hold on; xlim([min(matrix_dims) max(matrix_dims)]); ylim([0 10]); title("Zoom in");
subplot(3,1,3); hold on; xlim([min(matrix_dims) max(matrix_dims)]); title("% of cases with >" + string(maxiter) + " iters");


for i = 1:5
    subplot(3,1,1); % median
    plot(matrix_dims,median(ALL_RESULTS(i,:,:,3),3), 'LineWidth',3);

    subplot(3,1,2); % zoom
    plot(matrix_dims,median(ALL_RESULTS(i,:,:,3),3),'LineWidth',3);

    subplot(3,1,3); % not converge
    plot(matrix_dims,sum(ALL_RESULTS(i,:,:,3)==maxiter+1,3)./(N/100),'LineWidth',3);
end
linestyleorder(my_syms);
colororder(my_colors);

subplot(3,1,2); 
legend('SPA','NNSVDLRC','NNDSVD','rand','QDR','Location','eastoutside');


fontsize(20,"points");
exportgraphics(f,"RESULTS/" + str_exp + "_iters.png");
savefig("RESULTS/" + str_exp + "_iters");



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