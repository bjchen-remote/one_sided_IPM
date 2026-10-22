function window = wallWindowL2(wallSlope,ops,settings)
%IPM.EVOLVE.WALLWINDOWL2 Fixed compact wall-slope L2 functional.

window = ipm.evolve.wallWindowMoment(wallSlope,ops,settings,2);
end
