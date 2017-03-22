function output=optim_kernel(data,series,maxdev,bcx_0,bcx_end,bcy_0,bcy_end)

% optimization kernel for use with the FSSFgui and FSSFkernel

datasize=size(data);
limit=floor(datasize(1)^0.5);
if limit>15
    limit=15;
else
end
output.Max_Dev=[]; k=2;
output.optim_order=[];
output.coeff=[];
output.zfitted=[];
output.optim_order=[];
for k=1:2:35
    iter_output=feval(@kernel,data,k,series,bcx_0,bcx_end,bcy_0,bcy_end);
    if (abs(iter_output.Max_Dev))>(abs(maxdev))
    output.Max_Dev=[output.Max_Dev; iter_output.Max_Dev];
    output.coeff=iter_output.coeff;
    output.zfitted=iter_output.zfitted;
    output.optim_order=[output.optim_order; k];   
    else
    output.Max_Dev=[output.Max_Dev; iter_output.Max_Dev];
    output.coeff=iter_output.coeff;
    output.zfitted=iter_output.zfitted;
    output.optim_order=[output.optim_order; k];   
    break    
    end
end

        
        
        
        