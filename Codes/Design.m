
%This code solves the optimization problem considered in the paper via a
%line search on the scalar mu. The codes requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. Any SDP solver can be used.    

clear all;

Lambda =[1 0; 0 sqrt(2)];

H=[0 1.1;
   1 0]; 

B=eye(2);

n=max(size(H));

m=min(size(B));


prompt = 'Insert the resolution for the line search ';

n_points=input(prompt);

if(n_points>1)

prompt = 'Insert the lower bound for the line search ';

mu_min = input(prompt);

prompt = 'Insert the uppper bound for the line search ';

mu_max = input(prompt);

mu_v=linspace(mu_min,mu_max,n_points);

else
    
prompt = 'Insert the value of mu ';

mu_min = input(prompt); 

mu_v=linspace(mu_min,mu_min,1);

end

l=max(size(mu_v));
gamma_v=nan(l,1);

for i=1:l

mu=mu_v(i);


P=diag(sdpvar(n,1));

Q=sdpvar(n,n,'full');

theta=sdpvar(1,1,'full');

c=sdpvar(1,1,'full'); 

Y=sdpvar(m,n,'full');

M=[Q+Q'+Lambda*P, -Q'*H-Y'*B', -Y'*B';
   (-Q'*H-Y'*B')', -exp(-mu)*P*Lambda, zeros(n,n);
   (-Y'*B')', zeros(n,n),-theta*eye(n)];

Optimization=[P eye(n); eye(n) c*eye(n)];

problem=[M<=-1e-8*eye(max(size(M))), theta>=0, P>=1e-9*eye(n), c>=0, Optimization>=0];
    
options=sdpsettings('solver','mosek','verbose',0);

solution=solvesdp(problem, theta+c,options);

if(solution.problem==0)

P=double(P);

Q=double(Q);
K=(inv(Q'*B))*double(Y);

theta=double(theta);

gamma_v(i)=sqrt(theta/min(eig(P)))*exp(mu/2);


end

%%
if(l>1)
    
[gamma_best,index_best]=min(gamma_v);

mu=mu_v(index_best);

P=diag(sdpvar(n,1));

Q=sdpvar(n,n,'full');

theta=sdpvar(1,1,'full');

c=sdpvar(1,1,'full'); 

Y=sdpvar(m,n,'full');

M=[Q+Q'+Lambda*P, -Q'*H-Y'*B', -Y'*B';
   (-Q'*H-Y'*B')', -exp(-mu)*P*Lambda, zeros(n,n);
   (-Y'*B')', zeros(n,n),-theta*eye(n)];

Optimization=[P eye(n); eye(n) c*eye(n)];

problem=[M<=-1e-8*eye(max(size(M))), theta>=0, P>=1e-9*eye(n), c>=0, Optimization>=0];
    
options=sdpsettings('solver','mosek','verbose',0);

solution=solvesdp(problem, theta+c,options);

if(solution.problem==0)

P=double(P);

Q=double(Q);
K=(inv(Q'*B))*double(Y);

theta=double(theta);

gamma_v(i)=sqrt(theta/min(eig(P)))*exp(mu/2);
    
end
end

global H B K Lambda

end

