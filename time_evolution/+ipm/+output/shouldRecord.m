function record = shouldRecord(state,nextOutput)
%IPM.OUTPUT.SHOULDRECORD Test scheduled and terminal output clocks.

tau = state.scale.canonicalTime;
time = state.config.time;
record = tau+10*eps(tau) >= nextOutput || ...
    tau >= time.finalTime || ...
    state.scale.physicalTime >= time.physicalFinalTime;
end
