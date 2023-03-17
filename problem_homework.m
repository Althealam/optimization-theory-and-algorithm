%% 第一题 %%
%（1）对偶单纯形法
clear;clc;
c=[3,2,4,8,0,0];
A=[-2,5,3,-5,1,0;-1,-2,-5,-6,0,1];
b=[3;-8];
ind_B=[5,6];
[xm,zm,Table,flag]=DualSimplex(A,b,c,ind_B);

%（2）大M法
clear;clc;
M=10000;
c=[3,2,4,8,0,0,M];
A=[-2,5,3,-5,1,0,1;1,2,5,6,0,-1,1];
b=[3;8];
ind_B=[5,7];
[xm,zm,Table,flag]=SimplexMin(A,b,c,ind_B,1);

%（3）两阶段法
clear;clc
c=[3,2,4,8,0,0,0,0];
A=[-2,5,3,-5,1,0,1,0;1,2,5,6,0,-1,0,1];
b=[3;8];
ind_B=[7,8];
[xm,zm,Table,flag]=SimplexMin(A,b,c,ind_B,2);

%% 第二题 %%
%（1）大M法
clear;clc;
M=10000;
c=[-3,1,3,-1,M,M,M];
A=[1,2,-1,1,1,0,0;1,-1,2,-1,0,1,0;2,-2,3,3,0,0,1];
b=[0;6;9];
ind_B=[5,6,7];
[xm,zm,Table,flag]=SimplexMin(A,b,c,ind_B,1);

%（2）两阶段法
clear;clc;
c=[-3,1,3,-1,0,0,0];
A=[1,2,-1,1,1,0,0;1,-1,2,-1,0,1,0;2,-2,3,3,0,0,1];
ind_B=[5,6,7];
b=[0;6;9];
[xm,zm,Table,flag]=SimplexMin(A,b,c,ind_B,2);

%% 第三题 %%
%（1）对偶单纯形法
clear;clc;
c=[2,-3,0,0,0];
A=[-2,1,1,1,0;-1,1,-1,0,1];
b=[-3;-2];
ind_B=[4,5];
[xm,zm,Table,flag]=DualSimplex(A,b,c,ind_B);

%（2）大M法
clear;clc;
M=10000;
c=[2,-3,0,0,0,M,M];
A=[2,-1,-1,-1,0,1,0;1,-1,1,0,-1,0,1];
b=[3;2];
ind_B=[6,7];
[xm,zm,Table,flag]=SimplexMin(A,b,c,ind_B,1);

%（3）两阶段法
clear;clc;
c=[2,-3,0,0,0,1,1];
A=[2,-1,-1,-1,0,1,0;1,-1,1,0,-1,0,1];
b=[3;2];
ind_B=[6,7];
[xm,zm,Table,flag]=SimplexMin(A,b,c,ind_B,2);