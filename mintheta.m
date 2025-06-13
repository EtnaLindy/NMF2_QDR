% We attempt to find the minimizer for the function

% E(theta) = sum_i du_i^2 sin^2(max(0,theta - psi - pi/2)) 
%          + sum_j dv_j^2 sin^2(max(0,phi - theta))

% such that theta in [0,pi/2]

%--------------------------------------------------


function [theta1,fmin1] = mintheta(psi,phi,du,dv)

% preprocessing; filter out unimportant angles and sort the rest

du0 = du(psi <= 0);
psi0 = psi(psi <= 0);

dv0 = dv(phi >= 0);
phi0 = phi(phi >= 0);

% bounds for theta
a1 = min(psi0)+pi/2; b1 = max(phi0);

du0 = du0(psi0+pi/2 <= b1);
psi0 = psi0(psi0+pi/2 <= b1);  m = size(psi0,1);
dv0 = dv0(phi0 >= a1);
phi0 = phi0(phi0 >= a1);

% sort all the angles, store indices so that we know which is which
thetas = [psi0 + pi/2 ; phi0];
[thetas,I] = sort(thetas);

% function definitions
f = @(theta) sum(du(theta >= psi+pi/2).^2 .* sin(theta-psi(theta >= psi+pi/2)-pi/2).^2) + sum(dv(phi>=theta).^2 .* sin(phi(phi>=theta) - theta).^2);
a_func = @(theta) sum(du(theta >= psi+pi/2).^2 .* (1 - 2*sin(psi(theta >= psi+pi/2)).^2)) + sum(dv(theta <= phi).^2 .* (1 - 2*cos(phi(theta <= phi)).^2));
b_func = @(theta) sum(du(theta >= psi+pi/2).^2  .* sin(2*psi(theta >= psi+pi/2))) - sum(dv(theta <= phi).^2 .* sin(2*phi(theta <= phi)));


if a1 >= b1 % already nonnegative, nothing needs to be done
    theta1 = b1; % choose the smallest possible
    fmin1 = 0; % error equal to zero
    return;
end

% initial values for alpha and beta
alpha = sum(du(thetas(1) >= psi+pi/2).^2 .* (1 - 2*sin(psi(thetas(1) >= psi+pi/2)).^2)) + sum(dv(thetas(1) <= phi).^2 .* (1 - 2*cos(phi(thetas(1) <= phi)).^2));
beta = sum(du(thetas(1) >= psi+pi/2).^2  .* sin(2*psi(thetas(1) >= psi+pi/2))) - sum(dv(thetas(1) <= phi).^2 .* sin(2*phi(thetas(1) <= phi)));

for i = 2:size(thetas,1)

    % compute the root of the derivative for [theta_(i-1) , theta_i]
    if alpha == 0
        sol = pi/4;
    else
        sol = atan(beta/alpha)/2;
        if sol < 0, sol = sol + pi/2; end
    end 

    if (sol <= thetas(i) && sol >= thetas(i-1)) % the unique minimum found

        theta1 = sol; fmin1 = f(theta1); 
        return;

    elseif (sol < a1) % the minimum can be found at the end point of the tight interval (numerical error)

        theta1 = a1; fmin1 = f(theta1); 
        return;

    elseif (sol > b1) % the minimum can be found at the end point of the tight interval (numerical error)

        theta1 = b1; fmin1 = f(theta1); 
        return;

    else % compute new alpha and beta

        if I(i) <= m % point is in psi
            alpha = alpha + du0(I(i))^2*(1-2*sin(psi0(I(i)))^2);
            beta = beta + du0(I(i))^2*sin(2*psi0(I(i)));
        else % point is in phi
            alpha = alpha - dv0(I(i)-m)^2*(1-2*cos(phi0(I(i)-m))^2);
            beta = beta + dv0(I(i)-m)^2*sin(2*phi0(I(i)-m));
        end

    end
end

theta1 = -1; fmin1 = -1;
fprintf("We should never end up here. Are you sure that the first columns are nonnegative? \n\n")

end

