function tests = testIpr
tests = functiontests(localfunctions);
end

function testBothLocalHelpers(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
folders={fullfile(root,'20230316','Unbounded'), ...
         fullfile(root,'20230316','MultifractalPerturbation')};
vectors={[1;0;0],ones(4,1)/2,[1,1i]/sqrt(2)};
for k=1:numel(folders)
    oldPath=path;
    restore=onCleanup(@()path(oldPath));
    addpath(folders{k});
    verifyEqual(testCase,which('Ipr'),fullfile(folders{k},'Ipr.m'));
    for j=1:numel(vectors)
        psi=vectors{j};
        expected=sum((psi.*conj(psi)).^2);
        verifyEqual(testCase,Ipr(psi),expected,'AbsTol',1e-15);
    end
    clear restore;
end
end

