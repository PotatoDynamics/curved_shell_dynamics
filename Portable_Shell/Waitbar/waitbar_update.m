function [status] = waitbar_update(t, waitbar_instance)
global tend t_old
status = 0;
if floor(t) > floor(t_old)
    t_old = floor(t);
    waitbar(t/tend, waitbar_instance, ['Simulation at t = ' num2str(floor(t))]);
    tic;
elseif toc >= 1
    waitbar(t/tend, waitbar_instance, ['Simulation at t = ' num2str(t)]);
    tic;
end
if getappdata(waitbar_instance, 'canceling')
    tend = t;
    status = 1;
end
end