function status = solver_output(t, Y, flag, waitbar_instance, cleanup)
status = 0;
if strcmp(flag, 'init')
elseif strcmp(flag, 'done')
    delete(cleanup);
else
    [status] = waitbar_update(t(end), waitbar_instance);
end
end