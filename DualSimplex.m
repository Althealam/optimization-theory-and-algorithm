function [xm,zm,Table,flag]=DualSimplex(A,b,c,ind_B)
%% 初始化 %%
[m,n]=size(A);
ind_N=setdiff(1:n,ind_B);
flag=4;
zm=0;xm=[];x=[];
Table=[length(b)+1,length(c)+1]; %定义单纯形表，行增加检验数，列增加基变量索引和右端系数

%% 迭代求解 %%
while flag==4
    %1.给定初始对偶可行的基本解
    x(ind_B)=b;x(ind_N)=0;
    Table(1:length(b),length(c)+2)=b; %单纯形表最左边为基变量下标，最右边为系数
    Table(1:length(b),1)=ind_B; %单纯形表最左边为基变量索引
    Table(1:length(b),2:length(c)+1)=A; %单纯形表中间的部分为系数矩阵
    Table(length(b)+1,length(c)+2)=c*x'; 
        
      %2.计算b=B-1*b，满足最优性准则时则停止计算，否则继续迭代
    if all(b>=0)
        xm=x;
        zm=c*xm';
        disp('有最优解x:');disp(xm);
        disp('最优值z:');disp(zm);
        flag=0;
        return;
    end
    
    %3.未满足最优性准则时则继续迭代
    %3.1确定出基变量
    [~,out1]=min(b); %r1为出基变量在ind_B中的索引
    out=ind_B(out1); %r为出基变量的下标
    
    %3.2确定出基变量后判断是否可以找到入基变量
    if all(b(out1)<0) & all(A(out1,:)>=0)
        disp('无可行解');
        flag=1;
        return
    end
    
    %3.3确定入基变量
    %（1）计算检验数
    cB=c(ind_B);
    sigma=zeros(1,n);
    sigma(ind_N)=cB*A(:,ind_N)-c(ind_N);
    sigma(ind_B)=0;
    Table(length(b)+1,1:(length(c)+1))=[0,sigma];
    disp(Table);
    matrixwrite(Table);
    
    %（2）确定最小检验数的下标
    theta=zeros(1,n);
    theta(ind_N)=sigma(ind_N)./A(out1,ind_N);
    theta(theta==-Inf)=Inf;
    theta(ind_B)=Inf;
    [~,in1]=min(theta(ind_N));
    in=in1; %k为入基变量的下标，r为出基变量的下标
    
    %（3）换基
    ind_B(out1)=in;
    ind_N=setdiff(1:n,ind_B);
    
    %3.3旋转运算
    b(out1)=b(out1)/A(out1,in); %出基变量的右端系数
    A(out1,:)=A(out1,:)/A(out1,in); %主行的元素
    for  i=1:m
            if i~=out1
                b(i)=b(i)-A(i,in)*b(out1);
                A(i,:)=A(i,:)-A(i,in)*A(out1,:);
            end
        end
end
