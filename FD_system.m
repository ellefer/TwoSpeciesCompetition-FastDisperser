%% parametrs 
restart = 0;
N = 200;
x = [-1:2/N:1];
dx = 2/N;
k=1;
d1 = 0.2; %p-laplacian u 
d2 = 0; %linear u
d3 = 0; %p-laplacian v
d4 =.3; %linear v
q1 = 0;
q2 = 0;%drift 
q = 0;
pu = 7/4; % p for u 
pv = 7/4; % p for v
eps = .0001;
%r = ones(N+1,1);
%r = [ones(floor((N+1)/6),1);.01*ones(ceil(5*(N+1)/6),1)];
r = 1.25*(-x+1)';
%% ICs 

if restart == 1


    v = [ones(floor((1*N+1)/6),1);zeros(ceil(5*(N+1)/6),1)];%x+1;%
    u =  [zeros(floor(5*(N+1)/6),1);.1*ones(ceil((N+1)/6),1)];%(ones(1,N+1));%exp(-10*(x+.5).^2)';

    %u = (x+1.1);
    %v = .25*(x+1.1); 
    t0 = 0;
else
    t0 = t(end);
end
tmax =50;

fast = 0;


%% UXX
uxx = zeros(N+1,N+1);
for i = 2:N
    for j = 1:N+1
        if i == j+1
            uxx(i,j) = 1/(dx^2);
        elseif i==j
            uxx(i,j) = -(2/dx^2);
        elseif i == j-1
            uxx(i,j) = 1/(dx^2);
        end
    end
end
uxx(1,1) = -(1+dx*q/(d1+d2))/(dx^2);
uxx(1,2) = 1/(dx^2);
uxx(N+1,N+1) = -1/(dx^2);
uxx(N+1,N) = 1/(dx^2);

vxx = uxx;
vxx(1,1) =-(1+dx*q/(d3 + d4))/(dx^2);

%% UX

ux = zeros(N+1,N+1);
for i = 2:N
    for j = 1:N+1
        if i == j+1
            ux(i,j) = -1/dx;
            ux(i,j+1) = 1/dx;
        end
    end
end
ux(1,1) = q/(d1+d2);
vx = ux;
vx(1,1) = q/(d3+d4);
ux(N+1,N+1) = 0;

%% Time Step



[t,uv] = ode45(@(t,y) RHS(y,N,ux,vx,uxx,vxx,d1,d2,d3,d4,q,pu, pv,eps,r), [t0 t0+tmax], [u,v]);

plotgap = floor(length(t)/100);
for i = 1:100
    figure(1)
    plot(x,r,'Color',[0.93,0.69,0.13],'Linewidth',2)
    hold on
    plot(x,uv((i-1)*plotgap+1,1:N+1),'--','Color',[0.00,0.45,0.74], 'Linewidth',2)
    plot(x,uv((i-1)*plotgap+1,N+2:end), 'Color',[0.85,0.33,0.10],'Linewidth',2)
    hold off
    legend("Resources","fast (p = "+ num2str(pu)+ ")" + " d1 = " + num2str(d1),"slow, (p = " + num2str(pv) + ")" + " d3 = " + num2str(d3) )
    %ylim([0 2])
    xlabel('x')
    ylabel(['u(',num2str(t((i-1)*plotgap+1)),',x)'])
    drawnow
    filename = "q"+ replace(num2str(q),".","-") + "d1" + replace(num2str(d1),".","-") ...
       + "d2" + replace(num2str(d2),".","-")   + "d3" + replace(num2str(d3),".","-") ...
       +  "d4" + replace(num2str(d4),".","-")   + "pu" + replace(num2str(pu),".","-") ...
       +  "pv" + replace(num2str(pv),".","-") + "t"+ ...
        num2str(floor(t((i-1)*plotgap+1)))+".fig";
    % if i == 1
    %     savefig("classic_nodrift/"+filename)
    % elseif i == 7
    %     savefig("classic_nodrift/"+filename)
    % elseif i == 100
    %     savefig("classic_nodrift/"+filename)
    % end
    %filename = "fast"+ num2str(fast)+ "q"+ replace(num2str(q),".","-")+ "t"+ ...
    %num2str(floor(t((i-1)*plotgap+1)))+"p1-1.fig";
    % if i == 1
    %     savefig("hom/"+filename)
    % elseif i == 25
    %     savefig("hom/"+filename)
    % elseif i == 100
    %     savefig("hom/"+filename)
    % end
    pause(.05)

end
u = uv(end,1:N+1);
v = uv(end,N+2:end);



function [res] = RHS(uv,N,ux,vx,uxx,vxx,d1,d2,d3,d4, q,pu,pv,eps,r)
    u = uv(1:N+1);
    v = uv(N+2:end);
    dx = 2/N;

    secu = uxx*u;
    firu = ux*u; % might be able to take first entry to the (1/p-1) power to make it p-Laplacian
    firu(1) =  firu(1)^(1/(pu-1));
    secu(1) = (-1/dx^2)*(u(1) + dx*firu(1)) + u(2)/dx^2;

    denomu = firu.^2 + eps;

    diffu = d1.*secu.*(denomu).^((pu-4)/2).*((pu-1).*firu.^2 + eps) + d2.*(uxx*u);
    ut = diffu - q*(firu) +u.*(r-u-v);
    
    secv = vxx*v;
    firv = vx*v;
    firv(1) =  firv(1)^(1/(pv-1));
    secv(1) = (-1/dx^2)*(v(1) + dx*firv(1)) + v(2)/dx^2;

    denomv = firv.^2 + eps;
    

    diffv = d3.*secv.*(denomv).^((pv-4)/2).*((pv-1).*firv.^2 + eps) + d4.*(vxx*v);
    vt = diffv - q*(firv) + v.*(r-u-v);

    res= [ut;vt];
end
