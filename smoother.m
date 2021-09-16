function new_Force=smoother(file_Force, sizeKernel)

 % z.B.: sizeKernel=5; anz schritte pro richtung
new_Force(1:sizeKernel,1)=file_Force(1:sizeKernel);
for a=(1+sizeKernel):1:(length(file_Force)-sizeKernel)
    index=(a-sizeKernel):1:(a+sizeKernel);
    new_Force(a,1)= mean(file_Force(index));
end
new_Force(end:(end+sizeKernel),1)=file_Force((end-sizeKernel):end);
end

