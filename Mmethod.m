function [xm,zm,Table,flag]=Mmethod(A,b,c,ind_B)
%% 初始化 %%
[m,n]=size(A);
ind_B0=ind_B;ind_N=setdiff(1:n,ind_B);
flag=4;
zm=0;xm=[];x=[];M=1e4;
Table=zeros(length(b)+1,length(c)+2);

    %% 迭代求解 %%
    while flag==4
        %1.定义初始解与目标函数值
        x(ind_B)=b;x(ind_N)=0;
        Table(1:length(b),length(c)+2)=b; %单纯形表最左边为基变量下标，最右边为系数
        Table(1:length(b),1)=ind_B; %单纯形表最左边为基变量索引
        Table(1:length(b),2:length(c)+1)=A; %单纯形表中间的部分为系数矩阵
        
        %2.求检验数
        sigma=zeros(1,n);
        cB=c(ind_B);
        sigma(ind_N)=cB*A(:,ind_N)-c(ind_N);
        sigma(ind_B)=0;
        Table(length(b)+1,1:(length(c)+1))=[0,sigma]; %Table的最后一行为检验数
        Table(length(b)+1,length(c)+2)=c(ind_B)*x(ind_B)'; %Table的右下角的数字为目标函数值
        disp(Table);
        matrixwrite(Table);

        %3.确定入基变量in
        for i=1:length(c)
            if i~=ind_B0 & i~=ind_B %保证不让人工变量进基
            [~,in1]=max(sigma(i)); 
            end
        end
        in=ind_N(in1); 
        %in1为入基变量索引，in为入基变量下标

        %4.确定出基变量out
        theta=b./A(:,in);
        theta(theta<0)=Inf; %将小于0的数字设置为无穷大
        [~,out1]=min(theta); 
        out=ind_B(out1);
        %out1为出基变量索引，out为出基变量下标（out>out1）

        %5.换基
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
        if all(A(:,in)<=0)
            x=[];
            flag=1;
            disp('有无界解');
            return;
        end
        
        %（3）判断无可行解
        if all(sigma<=0) && x(ind_B0)~=0
            flag=2;disp('无可行解');
            return;
        end

        %7.旋转运算，计算新的系数矩阵等数值
        %(1)先变换非基变量的矩阵
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
end
