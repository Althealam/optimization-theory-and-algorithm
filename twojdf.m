function [xm,zm,Table,flag]=twojdf(A,b,c,ind_B)
%使用条件：有单位子块；右边元素非负；等式约束

%% phase1 %%
%1.初始化%
[m,n]=size(A);
ind_N=setdiff(1:n,ind_B);
flag=4;
zm=0;xm=[];x=[];
Table=zeros(length(b)+1,length(c)+2);
ind_B0=ind_B; %定义一个ind_B0用于phase2中去掉人工变量
c0=[];c0(ind_B)=1;c0(ind_N)=0;
%设置c0为phase1中目标函数的价格系数

%2.求解目标函数（人工变量求和）的最优解并保留其进行行变换后的矩阵
while flag==4
        %1.定义初始解
        x(ind_B)=b;x(ind_N)=0;
        Table(1:length(b),length(c0)+2)=b; %单纯形表最左边为基变量下标，最右边为系数
        Table(1:length(b),1)=ind_B; %单纯形表最左边为基变量索引
        Table(1:length(b),2:length(c0)+1)=A; %单纯形表中间的部分为系数矩阵

        %2.求检验数
        cB=c0(ind_B);
        sigma=zeros(1,n);
        sigma(ind_N)=cB*A(:,ind_N)-c0(ind_N);
        sigma(ind_B)=0;
        Table(length(b)+1,1:(length(c0)+1))=[0,sigma]; %Table的最后一行为检验数
        Table(length(b)+1,length(c0)+2)=cB*x(ind_B)'; %Table的右下角的数字为目标函数值
        disp(Table);
        matrixwrite(Table);

        %3.确定入基变量in
        [~,in1]=max(sigma(ind_N)); 
        in=ind_N(in1); 

        %4.确定出基变量out
        theta=b./A(:,in);
        theta(theta<0)=Inf; %将小于0的数字设置为无穷大
        [~,out1]=min(theta); 
        out=ind_B(out1);
        %注意：k1与r1为sigma与theta对应的分量索引
        %     k1=k为主列，r为主行

        %5.换基
        %入基变量为k=x3，出基变量为r=x5
        %入基变量索引为k1=2，出基变量索引为r1=2
        ind_B(out1)=in; %将出基变量out换为入基变量in
        ind_N=setdiff(1:n,ind_B);
        
        %6.判断解的情况
        %（1）判断是否为最优解
        if all(sigma<=0)
            xm=x;
            zm=c*xm';
            flag=0;
            disp('有有限最优解x:');disp(xm);
            disp('最优值z:');disp(zm);
            return;
        end

        %（2）判断无界解
        if all(A(:,in)<=0) && all(sigma<=0)
            x=[];
            flag=1; %当对应的检验数最大的一列的分量非正数时，问题不存在有限最优解
            disp('有无界解');
            return;
        end

        %（3）判断无可行解
        if all(sigma<=0) && x(ind_B0)~=0
            flag=2;disp('无可行解');
            return;
        end
%在满足标准型的线性规划中不会出现这种情况

        %7.旋转运算，计算新的系数矩阵等数值
%         %(1)先变换非基变量的矩阵
%         A(:,ind_N)=A(:,ind_B)\A(:,ind_N);
%         %(2)新的右端变量系数
%         b=A(:,ind_B)\b;
%         %(3)变换基变量的矩阵
%         A(:,ind_B)=eye(m,m);
     b(out1)=b(out1)/A(out1,in); %出基变量的右端系数
        A(out1,:)=A(out1,:)/A(out1,in); %主行的元素
        %主元原来的值为 A(out1,in1)
        h=[];
        for  i=1:m
            if i~=out1
                b(i)=b(i)-A(i,in)*b(out1);
                A(i,:)=A(i,:)-A(i,in)*A(out1,:);
            end
        end
     end
%注意：两阶段法中的phase1就相当于改变c后利用原始单纯形法

%% phase2 %%
%1.求出去掉人工变量后的系数矩阵与价格系数
A(:,ind_B0)=[];
c(ind_B0)=[];

%2.重新初始化
[m,n]=size(A);
ind_N=setdiff(1:n,ind_B);

%3.再次利用原始单纯形法
[xm,zm,flag]=dcxf(A,b,c,ind_B);
end




