% Ground truth of model's parameters

Flag.inxRun.init = (inxRun == 1);

if (Flag.inxRun.init)
    para = struct();
end

Radius = setting.Radius(inxRun)

para.MU = [-1 1; 1  1; 1 -1; -1 -1]' * Radius; 

para.offset = [1 1] * setting.offset;

para.MU(1,:) = para.MU(1,:) + para.offset(1);
para.MU(2,:) = para.MU(2,:) + para.offset(2);

para.SIGMA = eye(2); % independent Gaussian clusters 