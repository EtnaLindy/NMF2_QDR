addpath 'hierclust2nmf_v2'
addpath 'NNSVD-LRC_v2'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testcase C: 4x4 integer matrices with known optima
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(51193);                 % some birthday
recompute_results = false;  % use precomputed values or compute from scratch

load datamatrices.mat;

N = size(originalmatrices,1);

m = 4; n = 4;
maxiter = 1000;
tol_pow = -5;
tol = 10^tol_pow;


if recompute_results

% first dim : SPA, NNSVDLRC, NNDSVD, rand, QDR
% second dim : sampled data points
% fourth dim: init relative error to optimum, relative error to optimum, 
%               ANLS iters, ANLS convergence, 
%               init time, ANLS time,
%               init relative error to U, relative error to U

ALL_RESULTS = zeros(5,N,8);

compl = 10;

for i = 1:N

    if mod(i,round(N/10)) == 0, fprintf("%d %%... ", compl); compl =  compl + 10; end

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
my_colors = ["#CC0000";"#FFCC00";"#663399";"#1656DD"];



% TIME PERFORMANCE

% plotting....
f = figure(1); clf;
f.Position = [10 10 900 900];


res_count = 50;

ylabels = ["SPA", "NNSVDLRC", "NNDSVD","rand","QDR"];
marker_char = ['*','-'];
inds_exclude_rand = [1,2,3,5];

for j = inds_exclude_rand
    [~,I] = sort(ALL_RESULTS(j,:,1),"descend");
    subplot(4,1,min(j,4)); hold on;
    for k = inds_exclude_rand
        plot(1:res_count, ALL_RESULTS(k,I(1:res_count),1),marker_char((j==k)+1),'LineWidth',3);
    end
    ylabel(ylabels(j));
    xlim([1 res_count]);
    colororder(my_colors);
end

legend('SPA','NNSVDLRC','NNDSVD','QDR','Location','eastoutside');
xlabel("index");



sgtitle("Worst-case scaled distance between N^* and init");
fontsize(20,"points");
exportgraphics(f,"RESULTS/" + str_exp + "_accinit.png");
savefig("RESULTS/" + str_exp + "_accinit");


% plotting....
f = figure(2); clf;
f.Position = [10 10 900 900];


res_count = 50;

ylabels = ["SPA", "NNSVDLRC", "NNDSVD","rand","QDR"];
marker_char = ['*','-'];

for j = inds_exclude_rand
    [~,I] = sort(ALL_RESULTS(j,:,2),"descend");
    subplot(4,1,min(j,4)); hold on;
    for k = inds_exclude_rand
        plot(1:res_count, ALL_RESULTS(k,I(1:res_count),2),marker_char((j==k)+1),'LineWidth',3);
    end
    ylabel(ylabels(j));
    xlim([1 res_count]);
    colororder(my_colors);
end

legend('SPA','NNSVDLRC','NNDSVD','QDR','Location','eastoutside');
xlabel("index");

sgtitle("Worst-case scaled distance to N^*");
fontsize(20,"points");
exportgraphics(f,"RESULTS/" + str_exp + "_acc.png");
savefig("RESULTS/" + str_exp + "_acc");



titles = ["SPA","NNSVDLRC","NNDSVD","rand","QDR"];
bound_str = ["> 0.1", "[0.1,1e-5[", "<= 1e-5"];
bounds = [inf, 1e-1, 1e-5, -inf];

data = zeros(5,3);

for j = 1:5
fprintf("%s:\n",titles(j));
    for i = 1:3
        data(j,i) = sum((ALL_RESULTS(j,:,2) > bounds(i+1)) .* (ALL_RESULTS(j,:,2) <= bounds(i))); 
        fprintf("%s : \t%d kpl\n", bound_str(i), data(j,i)); 
    end
    fprintf("\n");
end




fprintf("Number of cases where certain method gives the closest initial point:\n");

[~,I] = min(ALL_RESULTS(:,:,1),[],1);

fprintf("SPA:\t\t%d\n", sum(I==1));
fprintf("NNSVDLRC:\t%d\n", sum(I==2));
fprintf("NNDSVD:\t\t%d\n", sum(I==3));
fprintf("rand:\t\t%d\n", sum(I==4));
fprintf("QDR:\t\t%d\n", sum(I==5));