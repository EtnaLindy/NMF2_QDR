
load datamatrices;

N = size(originalmatrices,1);

for i = 1:N
    U = reshape(originalmatrices(i,:),4,4);
    [L,S,R]=svds(U,2);
    U0 = L*S*R';
    R=R*sqrt(S);
    L=L*sqrt(S);
    if R(1,1)<0, R=-R;L=-L; end
    Psi0=atan2(L(:,2),L(:,1));
    Phi0=atan2(R(:,2),R(:,1));
    Du0 = sqrt(sum(L.^2,2));
    Dv0 = sqrt(sum(R.^2,2));
    figure(1); clf;
    Thetaplot(Psi0,Phi0,Du0,Dv0);
    
    alpha1 = mintheta(Psi0,Phi0,Du0,Dv0);
    alpha2 = mintheta(Phi0,Psi0,Dv0,Du0);
    Psi = min(alpha2,max(alpha1-pi/2,Psi0));
    Phi = min(alpha1,max(alpha2-pi/2,Phi0));
    Du = Du0.*cos(Psi-Psi0);
    Dv = Dv0.*cos(Phi-Phi0);
    figure(2); clf;
    Thetaplot(Psi,Phi,Du,Dv);

    disp(cos(Psi));
    
    Ubest = reshape(solutions(i,:),4,4);
    
    s=sqrt(cos(alpha2-alpha1)); 
    psi_new = min(max(Psi,alpha1-pi/2),alpha2);
    phi_new = min(max(Phi,alpha2-pi/2),alpha1);
    du = Du.*cos(Psi-psi_new);
    dv = Dv.*cos(Phi-phi_new);
    L=du.*[cos(alpha1-psi_new),sin(alpha2-psi_new)]/s;
    R=dv.*[cos(alpha2-phi_new),sin(alpha1-phi_new)]/s;
    D=diag(diag(L'*L)./diag(R'*R))^(1/4); %balance
    L=L/D;R=R*D;
    
    figure(3); clf;
    PlotAngles(L,R);

    pause;

end



function PlotAngles(L,R)

    Psi = atan2(L(:,2),L(:,1));
    Phi = atan2(R(:,2),R(:,1));
    Du = sqrt(sum(L.^2,2));
    Dv = sqrt(sum(R.^2,2));

    scale=max(max(Du),max(Dv));
    % plot a circle of radius scale in 100 dots
    theta=linspace(0,2*pi+0.1,100);
    plot(scale*exp(1i*theta),'k');
    hold on;
    scatter(Du.*cos(Psi),Du.*sin(Psi),'db','filled');
    scatter(Dv.*cos(Phi),Dv.*sin(Phi),'or','filled');
    plot(scale*[-1,1],[0,0],'-k');
    plot([0,0],scale*[-1,1],'-k');
    pbaspect([1 1 1]); xlim([-scale scale]); ylim([-scale scale]);
    xlabel('First component');
    ylabel('Second component');
    legend('Scaled circle','U','V');
    legend('Location','southwest');
    title('Rows of U and V with angular coordinates');
    hold off;

end

function [scale] = Thetaplot(Psi,Phi,Du,Dv)
% 
%  function [scale] = Thetaplot(Psi,Phi,Du,Dv)
%  
%  plots the points (given in in polar coordinates) 
%    Du.*exp(i.Psi) and Dv.*exp(i.Phi) 
%  on a disc of radius scale=max(Max(Du),max(Dv))
%
%  Paul Van Dooren, 22 November 2024
%
scale=max(max(Du),max(Dv));
% plot a circle of radius scale in 100 dots
theta=linspace(0,2*pi+0.1,100);
plot(scale*exp(1i*theta),'k');
hold on;
scatter(Du.*cos(Psi),Du.*sin(Psi),'db','filled');
scatter(Dv.*cos(Phi),Dv.*sin(Phi),'or','filled');
plot(scale*[-1,1],[0,0],'-k');
plot([0,0],scale*[-1,1],'-k');
plot([0,scale*cos(max(Psi)-pi/2)],[0,scale*sin(max(Psi)-pi/2)],':r','LineWidth',2);
plot([0,scale*cos(max(Psi))],[0,scale*sin(max(Psi))],':b','LineWidth',2);
plot([0,scale*cos(max(Phi)-pi/2)],[0,scale*sin(max(Phi)-pi/2)],':b','LineWidth',2);
plot([0,scale*cos(max(Phi))],[0,scale*sin(max(Phi))],':r','LineWidth',2);
viol = abs(Psi-Phi') > pi/2;
inds = reshape(1:16,4,4);
I = inds(viol);
if size(I,1) > 0
psi_inds = mod(I-1,4)+1;
phi_inds = round((I+1)/4);
scatter(Du(psi_inds).*cos(Psi(psi_inds)),Du(psi_inds).*sin(Psi(psi_inds)),100,'db');
scatter(Dv(phi_inds).*cos(Phi(phi_inds)),Dv(phi_inds).*sin(Phi(phi_inds)),100,'or');
end
pbaspect([1 1 1]); xlim([-scale scale]); ylim([-scale scale]);
xlabel('First component');
ylabel('Second component');
legend('Scaled circle','U','V');
legend('Location','southwest');
title('Rows of U and V with angular coordinates');
hold off;

end