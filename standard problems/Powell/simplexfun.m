function [x,z,flg,sgma]=simplexfun(A,A1,b,c,m,n,n1,cb,xx)

% A,b are the matric in A*x=b

% c is the matrix in max z=c*x

% A1 is the matric in simplex table

% m is the numbers of row in A and n is the con number in A

% n1 is the nubers of artificial variables,and artificial variables are default as the last % n1 variables in x.

% cb is the worth coefficient matrix for basic variables

% xx is the index matrix for basic variables

% B1 is the invers matrix for the basic matrix in simplex table.The initial

% matrix is default as the last m con in the matrix A.

x=zeros(n,1);

z=0;

B1=A1(:,n-m+1:n);

sgma1=c-(cb*B1)*A;

[masg,kk]=max(sgma1);

k=kk(1);

flg=0;

ll=0;

while (masg>0)&&(ll<20)

ll=ll+1;

thita=1000+zeros(m,1);

for i=1:m

if A1(i,k)>0

thita(i)=A1(i,k)\b(i);

end

end

[r8,c8]=find(thita>999);

if sum(c8)

[mith,rr]=min(thita);

r=rr(1);

aa=A1(r,k);

for i=1:m

if i==r

b(r)=b(r)/aa;

for j=1:n

A1(r,j)=A1(r,j)/aa ;

end

end

end

for i=1:m

if i~=r

cc=A1(i,k)

b(i)=b(i)-b(r)*cc;

for j=1:n

A1(i,j)=A1(i,j)-A1(r,j)*cc;

end

end

end

cb(r)=c(k);

xx(r)=k;

B1=A1(:,n-m+1:n);

sgma1=c-(cb*B1)*A;

[masg,kk]=max(sgma1);

k=kk(1);

thita=100+zeros(m,1);

else

flg=3;

masg=-1;

x='unbound solution';

z='inf';

end

end

if flg~=3

if n1==0

sgma1=c-(cb*B1)*A

[rc,ccc]=find(sgma1

if sum(rc)==n-m

flg=1;

else

flg=2;

end

x=zeros(n,1);

for i=1:m

x(xx(i))=b(i);

end

z=c*x;

else

