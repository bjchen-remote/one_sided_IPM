function plotResult(result)
%IPM.OUTPUT.PLOTRESULT Plot the terminal physical fields and trusted diagnostics.

result = ipm.output.validate(result);
gridData = result.grid;
physical = result.physical;
common = result.history.common;
figure('Name','IPM on the upper half-plane');
tiledlayout(2,2,'TileSpacing','compact');
positiveSide = gridData.physicalX >= 0;
if ~any(positiveSide)
    positiveSide = true(size(gridData.physicalX));
end
displayX = gridData.comovingPhysicalX(positiveSide);
[displayGridX,displayGridY] = meshgrid(displayX,gridData.physicalY);

axisHandle = nexttile;
surf(axisHandle,displayGridX,displayGridY, ...
    physical.omega(:,positiveSide),'EdgeColor','none');
view(axisHandle,2); axis(axisHandle,'image'); axis(axisHandle,'tight');
colorbar(axisHandle);
xlabel('physical x_1 relative to peak'); ylabel('physical x_2');
title(sprintf('physical rho_{x_1} on x_1>=0, tau=%.4g', ...
    result.state.canonicalTime));

axisHandle = nexttile;
speed = hypot(physical.u1,physical.u2);
surf(axisHandle,displayGridX,displayGridY, ...
    speed(:,positiveSide),'EdgeColor','none');
view(axisHandle,2); axis(axisHandle,'image'); axis(axisHandle,'tight');
colorbar(axisHandle);
xlabel('physical x_1 relative to peak'); ylabel('physical x_2');
title('physical |u| in peak frame');

nexttile;
semilogy(common.physicalTime,common.physicalGradInf,'LineWidth',1.5);
grid on; xlabel('physical t'); ylabel('||grad rho||_infinity');
title('gradient growth');

nexttile;
yyaxis left;
plot(common.canonicalTau,common.canonicalCX,'LineWidth',1.5);
hold on;
plot(common.canonicalTau,common.canonicalCY,'LineWidth',1.5);
plot(common.canonicalTau,common.canonicalCOmega,'LineWidth',1.5);
ylabel('canonical scaling rates');
yyaxis right;
plot(common.canonicalTau,common.gradientInnerFraction, ...
    '--','LineWidth',1);
hold on;
plot(common.canonicalTau,common.gradientNearFraction, ...
    '--','LineWidth',1);
plot(common.canonicalTau,common.gradientOuterFraction, ...
    '--','LineWidth',1);
ylim([0,1]); ylabel('gradient-energy fraction');
grid on; xlabel('canonical tau');
legend('c_x','c_y','c_omega','E inner','E near','E outer', ...
    'Location','best');
title('rates and scale transfer');
end
