function active = isActive(state)
%IPM.EVOLVE.ISACTIVE True while both clocks and the step budget remain open.

time = state.config.time;
active = state.scale.canonicalTime < time.finalTime && ...
    state.scale.physicalTime < time.physicalFinalTime && ...
    state.step < time.maxSteps;
end
