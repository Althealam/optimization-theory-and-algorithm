%% 大M法 %%
%(1)
clear;clc;
M=10000;
c=[1,1,-3,0,0,M,M];
A=[1,-2,1,1,0,0,0;2,1,-4,0,-1,1,0;1,0,-2,0,0,0,1];
b=[11;3;1];
ind_B=[4,6,7];
[xm,zm,Table,flag]= Mmethod(A,b,c,ind_B);

%% 单纯形法 %%
clear;clc;
c=[1,-2,1,0,0,0]; %最小化问题的价格系数
A=[1,1,-2,1,0,0;2,-1,4,0,1,0;-1,2,-4,0,0,1];
b=[10;8;4];
ind_B=[4,5,6];
[xm,zm,Table,flag]=dcxf(A,b,c,ind_B);

%% 两阶段法 %%
clear;clc;
c=[-2,1,0,0,0,0,0,0]; %phase2 （去掉人工变量）
A=[1,1,-1,0,0,1,0,0;
    1,-1,0,-1,0,0,1,0;
    1,0,0,0,1,0,0,1];
b=[2;1;3];
ind_B=[6,7,8]; %phase1中的基变量索引
[xm,zm,flag]=twojdf(A,b,c,ind_B);

%% 对偶单纯形法 %%
clear;clc;
c=[12,8,16,12,0,0];
b=[-2;-3];
A=[-2,-1,-4,0,1,0;-2,-2,0,-4,0,1];
ind_B=[5,6];
[xm,zm,flag]=DualSimplex(A,b,c,ind_B);