function [xm,zm,Table,flag]=Simplex(A,b,c,ind_B,method)

if method==0 %原始单纯形法
    [xm,zm,Table,flag]=dcxf(A,b,c,ind_B);
end

if method==1 %大M法
    [xm,zm,Table,flag]=Mmethod(A,b,c,ind_B);
end

if method==2 %两阶段法
    [xm,zm,Table,flag]=twojdf(A,b,c,ind_B);
end
    
    