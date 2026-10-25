function test_report(folder,stem,report)
%TEST_REPORT Atomically publish one named MAT/JSON report, without arithmetic.
% Assertions and scientific acceptance rules remain in each independent test.
assert(ischar(stem)&&~isempty(regexp(stem,'^[A-Za-z][A-Za-z0-9_]*$','once')), ...
    'ipm_cr:TestReportStem','Use one report basename.');
temporary=fullfile(folder,[stem,'.writing.mat']);save(temporary,'report','-v7.3');
move_checked(temporary,fullfile(folder,[stem,'.mat']));
ipm_cr.io.json(fullfile(folder,[stem,'.json']),report);
end
function move_checked(from,to)
[ok,message]=movefile(from,to,'f');assert(ok,'ipm_cr:TestReportWrite','%s',message);
end
