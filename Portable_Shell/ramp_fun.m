function ramp = ramp_fun(f_start,df,f_end,t_start,t_end,t)
    
    f_vec = [f_start:df:f_end];
    dt = (t_end-t_start)/length(f_vec);
    
    if(t  < t_start)
        ramp = f_start;  
    elseif((t >= t_start) && (t < t_end))
        ramp = f_vec(floor((t-t_start)/dt)+1); % upwards sweep
    elseif(t >= t_end)
        ramp = f_end;
    end  
    
end