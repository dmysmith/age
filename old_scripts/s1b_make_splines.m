% Create smoothing splines for age analysis

% Load nda file
ndafile = '/space/syn50/1/data/ABCD/d9smith/age/nda5.0.txt';
nda = readtable(ndafile);

agevec = nda.interview_age;
knots = [10:2:20]; % Todo - where to set knots?
pp = spline(knots, eye(length(knots)));
bf  = ppval(pp, agevec);
% mean(bf,2, "omitnan")
bf_demean = bf - mean(bf,2, "omitnan");

bfstr = cell(length(knots), 1);
bfstr_demean = bfstr;
for i=1:length(knots)
    bfstr{i} = ['bf', num2str(i)];
    bfstr_demean{i} = ['bf_demean_', num2str(i)];
end

T = [nda, array2table(bf', 'VariableNames',bfstr), array2table(bf_demean', 'VariableNames',bfstr_demean)];

outfile = '/space/syn50/1/data/ABCD/d9smith/age/nda5.0_withbfs.txt';
writetable(T, outfile);

if 0 % code from Anders
    agevec = linspace(10,20,101);
    knots = [10:2:20];
    pp = spline(knots, eye(length(knots)));
    bf  = ppval(pp, agevec);

    figure; imagesc(bf);
    figure; plot(bf');
end