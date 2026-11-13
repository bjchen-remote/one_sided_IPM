function weight=compactCutoff(x,start,stop)
%COMPACTCUTOFF Smooth cutoff equal to one before START and zero after STOP.
u=(x-start)/(stop-start);
weight=double(u<=0);
inside=u>0&u<1;
z=u(inside);
rise=exp(-1./z);
fall=exp(-1./(1-z));
weight(inside)=fall./(rise+fall);
end
