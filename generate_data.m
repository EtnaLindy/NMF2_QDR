% Generate mxn datamatrix
% onlynegative : bool           - regenerate the matrix until rank-2 approximation is not nonnegative (default false)
% datatype : {  "lowrank",      - rank-2 nonnegative matrix with noise
%               "integer",      - integer matrix with entries summing up to 1000
%               "simplex",      - stochastic matrix
%               default         - lognormal independent entries
%             }


function [U] = generate_data(m,n,onlynegative,datatype,noise_std)
    
if nargin < 3, onlynegative=false; end
if nargin < 4, datatype = "generic"; end
if nargin < 5, noise_std = 0.2; end

while true % keep generating until rank-2 approximation is not nonnegative

if datatype == "lowrank"
    
    p = round(0.5*m);
    q = round(0.5*n);

    Psi = zeros(p,1);
    Psi = [Psi ; rand(m-p,1)*pi/2]; Psi = Psi(randperm(m));
    Phi = pi/2*ones(q,1);
    Phi = [rand(n-q,1)*pi/2 ; Phi]; Phi = Phi(randperm(n));

    Du = rand(m,1); % weights
    Dv = rand(n,1);

    L0 = round([cos(Psi),sin(Psi)],3);
    R0 = round([cos(Phi),sin(Phi)],3);

    noise = abs(randn(m,n))*noise_std;

    U = Du.*L0*(Dv.*R0)' + noise;
    

elseif datatype == "integer"

    maxval = 1000;
    x0 = [0 ; sort(randi(maxval,n*m-1,1)) ; maxval];
    U0 = zeros(m*n,1);
    for j = 1:m*n
        U0(j) = x0(j+1)-x0(j);
    end
    U = reshape(U0(randperm(n*m)),m,n);

elseif datatype == "simplex"

    x0 = [0 ; sort(rand(n*m-1,1)) ; 1];
    U0 = zeros(m*n,1);
    for j = 1:m*n
        U0(j) = x0(j+1)-x0(j);
    end
    U = reshape(U0,m,n);

else

    U = exp(sqrt(log(max(m,n)))*randn(m,n));

end

if ~(onlynegative && is_nonneg2(U))
    return;
end


end

end

