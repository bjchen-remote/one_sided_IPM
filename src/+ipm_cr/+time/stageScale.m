function s=stageScale(previous,z,physicalIncrement)
%STAGESCALE One scale representation for RHS evaluation and accepted commit.
% The RK physical-time increment is supplied before it is added to the large
% clock. Recovering it from z(3)-previous.physicalTime would lose small steps.
validateattributes(physicalIncrement,{'numeric'},{'scalar','real','finite','nonnegative'});
assert(numel(z)==6&&all(isfinite(z)),'ipm_cr:ExplicitScale','Use six finite RK scalars.');
s=previous;s.logC_l=z(1);s.logC_omega=z(2);s.X_shift=z(4);
s.canonicalTime=z(5);s.densityOrigin=z(6);
high=previous.physicalTime;low=previous.physicalTimeLow;
sumHigh=high+physicalIncrement;b=sumHigh-high;
error=(high-(sumHigh-b))+(physicalIncrement-b);
t=low+error;s.physicalTime=sumHigh+t;
s.physicalTimeLow=(sumHigh-s.physicalTime)+t;
end
