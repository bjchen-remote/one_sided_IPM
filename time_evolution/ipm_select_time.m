function data = ipm_select_time(data,selector,mode)
%IPM_SELECT_TIME Select a saved native frame for ordinary MATLAB plotting.
if nargin<3, mode='nearest'; end
data=ipm.output.selectTime(data,selector,mode);
end
