% DROPLET_INCLINE30_GRAVITY
% 单文件 MATLAB 程序：真实 30 度斜面 + 竖直向下重力 + Shan-Chen 型伪势 LBM。
% 在 MATLAB 中打开本文件，点击 Run。无需其他 .m 文件或工具箱。
% 参考 Mohammed Fouad Kamel Ibrahim / UNSW 的 D2Q9 MRT + PR EOS 算例。
% 保留其 PR/MRT、tune 和 Q 框架；斜壁、润湿、重力和后处理在此重新实现。
% 这是二维截面模型，不是三维液滴或 IBM。外框为封闭中性反弹边界。
% 坡面由左上向右下倾斜，重力在实验室坐标中为 (0,-g)。
% 默认几何虚拟密度润湿边界使用你标定的平衡角；不是动态角经验公式。
% 斜面边界已经不同于原平壁代码，仍须独立做静态角与网格收敛验证。
% 运行耗时只在屏幕/内存显示，不写入结果文件。

cfg.steps = 40000;
cfg.nx = 402;
cfg.ny = 302;
cfg.angleDeg = 30;
cfg.wallHeightLeft = 260;       % 真实斜线 y=260-x*tan(30 deg)。
cfg.dropContactX = 100;         % 初始圆帽的对称轴与坡面的交点 x。
cfg.radius = 32;
cfg.initialAngleDeg = 70;       % 仅初始化圆帽形状，后续接触角由流动决定。
cfg.interfaceWidth = 5;
cfg.Gw = -0.8;                 % geometric模式仅是标签；source_force模式才使用此力参数。
cfg.wettingMode = 'geometric'; % 'geometric'：角度控制；'source_force'：原壁面力对照。
cfg.equilibriumAngleDeg = 69.107172; % 默认几何边界用此角，Gw仅作标定标签。
cfg.ghostReferenceDistance = 2; % 虚拟固体密度从流体侧参考面插值。
cfg.gravity = 1e-5;            % 格子单位加速度；未做 SI 单位映射。
cfg.gravityStart = 10000;      % 前10000步无重力弛豫。
cfg.gravityRamp = 2000;        % 然后2000步余弦平滑加载。
cfg.relaxationMode = 'source'; % 'source': 末流体节点全域松弛；'local': 局部黏度。
cfg.wallScheme = 'interpolated'; % 'interpolated': Bouzidi; 'halfway': 阶梯反弹对照。
cfg.conserveBoundaryMass = true; % 插值边界的质量差补到静止分布，不改动量。
cfg.plot = true;
cfg.plotEvery = 200;
cfg.sampleEvery = 200;
cfg.saveEvery = 5000;
cfg.runtimeUpdateSeconds = 1;
cfg.fitMinDistance = 2.5;      % 相对真实斜面的法向距离，单位：格点。
cfg.fitMaxDistance = 9;
cfg.fitOrder = 2;
cfg.boundaryClearance = 8;     % 液滴接近外框时停止，避免把碰壁误认为坡面滑动。
cfg.maxAllowedSpeed = 0.25;    % 明显离开低马赫数范围时停止并提示。
cfg.maxRelativeMassDrift = 1e-3; % 质量漂移超过0.1%时停止，检查边界方案。
cfg.returnPopulations = false;
cfg.outputDir = fullfile(fileparts(mfilename('fullpath')),'incline30_results');
% cfg.outputDir = '';         % 只显示，不写任何结果。

result = run_incline_case(cfg);
disp(result.history);
fprintf('Stopped at step %d: %s\n',result.step,result.stopReason);
fprintf('Wall-clock elapsed: %s (%.2f s)\n',...
    elapsed_string(result.elapsedSeconds),result.elapsedSeconds);
if ~isempty(cfg.outputDir), fprintf('Results: %s\n',cfg.outputDir); end

function result = run_incline_case(cfg)
runClock = tic;
validate_config(cfg);
p = cfg;
p.rhoV = 0.34011; p.rhoL = 6.34793;
p.nuV = 0.27462; p.nuL = 0.065;
p.Tc = 0.0729; p.T = 0.85*p.Tc;
p.a = 2/49; p.b = 2/21; p.R = 1; p.omega = 0.344;
p.cs = 1/sqrt(3); p.cs2 = p.cs^2; p.G = -1;
p.tune = 0.107; p.KQ = 0.01;
p.betaPR = (1+(0.37464+1.54226*p.omega-0.26992*p.omega^2)*...
    (1-sqrt(p.T/p.Tc)))^2;
cx = [0 1 0 -1 0 1 -1 -1 1];
cy = [0 0 1 0 -1 1 1 -1 -1];
opposite = [1 4 5 2 3 8 9 6 7];
w = [4/9 repmat(1/9,1,4) repmat(1/36,1,4)];
wi = [0 repmat(1/3,1,4) repmat(1/12,1,4)];
M = [1 1 1 1 1 1 1 1 1; -4 -1 -1 -1 -1 2 2 2 2;
    4 -2 -2 -2 -2 1 1 1 1; 0 1 0 -1 0 1 -1 -1 1;
    0 -2 0 2 0 1 -1 -1 1; 0 0 1 0 -1 1 1 -1 -1;
    0 0 -2 0 2 1 1 -1 -1; 0 1 -1 1 -1 0 0 0 0;
    0 0 0 0 0 1 -1 1 -1];
Minv = bsxfun(@rdivide,M.',[9 36 36 6 12 6 12 4 4]);
[X,Y] = ndgrid(0:p.nx-1,0:p.ny-1);
c = cosd(p.angleDeg); sn = sind(p.angleDeg);
normalDistance = sn*X+c*(Y-p.wallHeightLeft);
alongCoordinate = c*X-sn*(Y-p.wallHeightLeft);
slopeSolid = normalDistance<=0;
outerSolid = false(p.nx,p.ny);
outerSolid([1 end],:) = true; outerSolid(:,[1 end]) = true;
solid = slopeSolid | outerSolid;
fluid = ~solid; ids = find(fluid); n = numel(ids);
[ix,iy] = ind2sub(size(fluid),ids);
neighbors = zeros(n,9); neighbors(:,1) = ids;
for k=2:9, neighbors(:,k) = ids+cx(k)+p.nx*cy(k); end
plan = streaming_plan(ids,neighbors,fluid,slopeSolid,outerSolid,...
    normalDistance,opposite,p.wallScheme);
wallPlan = virtual_wall_plan(neighbors,fluid,slopeSolid,outerSolid,...
    normalDistance,X,Y,p);

% Initial smooth circular cap, clipped by the actual geometric plane.
contact = [p.dropContactX; p.wallHeightLeft-p.dropContactX*tand(p.angleDeg)];
center = contact-p.radius*cosd(p.initialAngleDeg)*[sn;c];
distance = hypot(X-center(1),Y-center(2));
fraction = 0.5*(1-tanh(2*(distance-p.radius)/p.interfaceWidth));
rho = zeros(size(X)); rho(fluid) = p.rhoV+(p.rhoL-p.rhoV)*fraction(fluid);
[pressure,psi] = pr_eos(virtual_density(rho,wallPlan,p),p,0);
[ax,ay,axFF,ayFF] = all_forces(rho,psi,neighbors,fluid,slopeSolid,ids,p,0);
d = rho(ids); u = zeros(n,1); v = u;
% Match zero physical velocity with the half-force definition at initialization.
f = equilibrium(d,-0.5*ax./d,-0.5*ay./d,cx,cy,w);
historyRows = zeros(ceil(p.steps/p.sampleEvery)+2,18); nh = 0;
initialMass = sum(d); initialAlong = NaN;
h = []; titlePrefix = ''; lastTimerUpdate = -Inf;
stopReason = 'requested step count reached';
if ~isempty(p.outputDir) && ~exist(p.outputDir,'dir'), mkdir(p.outputDir); end
fprintf('Actual %.1f deg plane; gravity=(0,-%.4g), Gw=%g\n',...
    p.angleDeg,p.gravity,p.Gw);
fprintf('Grid %d x %d, %d fluid nodes, wall=%s, relaxation=%s\n',...
    p.nx,p.ny,n,p.wallScheme,p.relaxationMode);
fprintf('Wetting=%s, prescribed equilibrium angle=%.4f deg (not dynamic angle).\n',...
    p.wettingMode,p.equilibriumAngleDeg);
fprintf('No gravity until step %d, ramp over %d steps.\n',p.gravityStart,p.gravityRamp);

for step=0:p.steps
    if step>0
        [qxx,qyy,qxy] = q_components(psi,ids,neighbors,cx,cy,wi,p);
        nu = p.nuL*(min(max(d,p.rhoV),p.rhoL)-p.rhoV)/(p.rhoL-p.rhoV)...
            +p.nuV*(p.rhoL-min(max(d,p.rhoV),p.rhoL))/(p.rhoL-p.rhoV);
        shear = 1./(3*nu+0.5);
        if strcmpi(p.relaxationMode,'source'), shear(:) = shear(end); end
        rates = repmat([1 .8 .8 1 1.1 1 1.1 1 1],n,1);
        rates(:,8) = shear; rates(:,9) = shear;
        u2 = u.^2+v.^2;
        meq = [d d.*(-2+3*u2) d.*(1-3*u2) d.*u -d.*u d.*v -d.*v ...
            d.*(u.^2-v.^2) d.*u.*v];
        work = u.*ax+v.*ay;
        correction = 12*p.tune*(axFF.^2+ayFF.^2)./...
            (psi(ids).^2*((1/.8)-.5)*p.cs*p.cs);
        forcing = [zeros(n,1) 6*work+correction -6*work-correction ...
            ax -ax ay -ay 2*(u.*ax-v.*ay) u.*ay+v.*ax];
        qc = zeros(n,9);
        qc(:,2) = 1.5*rates(:,2).*(qxx+qyy);
        qc(:,3) = -1.5*rates(:,3).*(qxx+qyy);
        qc(:,8) = -shear.*(qxx-qyy); qc(:,9) = -shear.*qxy;
        moments = ordered_transform(f,M);
        moments = moments-rates.*(moments-meq)+(1-rates/2).*forcing+qc;
        post = ordered_transform(moments,Minv);
        f = stream_and_reflect(post,plan,opposite,p.conserveBoundaryMass);
        d = zeros(n,1); mx = d; my = d;
        for k=1:9
            d = d+f(:,k); mx = mx+cx(k)*f(:,k); my = my+cy(k)*f(:,k);
        end
        if any(~isfinite(d)) || any(d<=0)
            error('Incline:Density','Invalid density at step %d. Reduce gravity/check wall settings.',step);
        end
        rho(ids) = d;
        [pressure,psi] = pr_eos(virtual_density(rho,wallPlan,p),p,step);
        [ax,ay,axFF,ayFF] = all_forces(rho,psi,neighbors,fluid,slopeSolid,ids,p,step);
        u = (mx+0.5*ax)./d; v = (my+0.5*ay)./d;
        if any(~isfinite(u)) || any(~isfinite(v))
            error('Incline:Velocity','Non-finite velocity at step %d.',step);
        end
    end
    sampleNow = step==0 || mod(step,p.sampleEvery)==0 || step==p.steps;
    plotNow = p.plot && (step==0 || mod(step,p.plotEvery)==0 || step==p.steps);
    saveNow = ~isempty(p.outputDir) && step>0 && ...
        (mod(step,p.saveEvery)==0 || step==p.steps);
    dense = d>(p.rhoV+p.rhoL)/2;
    hitBoundary = any(dense & (ix<=p.boundaryClearance | ix>=p.nx-p.boundaryClearance ...
        | iy>=p.ny-p.boundaryClearance));
    tooFast = max(hypot(u,v))>p.maxAllowedSpeed;
    massDrift = (sum(d)-initialMass)/initialMass;
    massProblem = abs(massDrift)>p.maxRelativeMassDrift;
    stopNow = hitBoundary || tooFast || massProblem || ~any(dense);
    if hitBoundary
        stopReason = 'liquid approached the finite outer boundary';
    elseif tooFast
        stopReason = 'speed exceeded the configured low-Mach check';
    elseif massProblem
        stopReason = 'boundary mass drift exceeded the configured tolerance';
    elseif ~any(dense)
        stopReason = 'no resolved dense-liquid core remained';
    end
    if stopNow, sampleNow = true; plotNow = p.plot; saveNow = ~isempty(p.outputDir); end
    if sampleNow || plotNow || saveNow
        m = observe_drop(rho,u,v,ids,fluid,X,Y,alongCoordinate,normalDistance,p);
        if isnan(initialAlong), initialAlong = m.alongCenter; end
        m.displacement = m.alongCenter-initialAlong;
    end
    if sampleNow
        nh = nh+1;
        historyRows(nh,:) = [step gravity_fraction(step,p) sum(d) ...
            (sum(d)-initialMass)/initialMass m.liquidArea m.centerX m.centerY ...
            m.alongCenter m.displacement m.meanAlongVelocity m.sReceding m.sAdvancing ...
            m.thetaReceding m.thetaAdvancing max(hypot(u,v)) min(d) ...
            m.fitRmsReceding m.fitRmsAdvancing];
        fprintf('%6d: g/g0=%.2f, ds=%+.3f, U_s=%+.4g, thetaRec/Adv=%.2f/%.2f, mass drift=%+.3g%%\n',...
            step,gravity_fraction(step,p),m.displacement,m.meanAlongVelocity,...
            m.thetaReceding,m.thetaAdvancing,100*(sum(d)-initialMass)/initialMass);
    end
    if plotNow
        if isempty(h) || ~ishghandle(h.figure)
            h = create_plot(rho,solid,X,Y,p);
        end
        titlePrefix = sprintf('Actual %.0f deg slope | step %d | g/g_0=%.2f | Rec/Adv=%.1f/%.1f deg',...
            p.angleDeg,step,gravity_fraction(step,p),m.thetaReceding,m.thetaAdvancing);
        set(h.image,'CData',rho.');
        set(h.interface,'XData',m.interface(1,:),'YData',m.interface(2,:));
        set(h.center,'XData',m.centerX,'YData',m.centerY);
        set(h.receding,'XData',historyRows(1:nh,1),'YData',historyRows(1:nh,13));
        set(h.advancing,'XData',historyRows(1:nh,1),'YData',historyRows(1:nh,14));
        set(h.displacement,'XData',historyRows(1:nh,1),'YData',historyRows(1:nh,9));
        set(h.title,'String',{titlePrefix; ['Wall-clock elapsed: ' elapsed_string(toc(runClock))]});
        drawnow;
    end
    if saveNow
        snapshot = make_result(step,rho,u,v,ids,pressure,solid,normalDistance,p,m,...
            historyRows(1:nh,:),stopReason);
        save(fullfile(p.outputDir,sprintf('snapshot_t%06d.mat',step)),'snapshot');
    end
    elapsed = toc(runClock);
    if elapsed-lastTimerUpdate>=p.runtimeUpdateSeconds && ~isempty(h) && ishghandle(h.title)
        set(h.title,'String',{titlePrefix; ['Wall-clock elapsed: ' elapsed_string(elapsed)]});
        drawnow limitrate; lastTimerUpdate = elapsed;
    end
    if stopNow
        warning('Incline:Stopped','Step %d: %s.',step,stopReason);
        break;
    end
end
result = make_result(step,rho,u,v,ids,pressure,solid,normalDistance,p,m,...
    historyRows(1:nh,:),stopReason);
if p.returnPopulations, result.fFluid = f; result.fluidIndices = ids; end
result.elapsedSeconds = toc(runClock);
if ~isempty(p.outputDir)
    diskResult = rmfield(result,'elapsedSeconds');
    fileData = struct('result',diskResult);
    save(fullfile(p.outputDir,'final_result.mat'),'-struct','fileData');
    writetable(result.history,fullfile(p.outputDir,'history.csv'));
    if ~isempty(h) && ishghandle(h.figure)
        set(h.title,'String',titlePrefix); % No elapsed runtime in saved images.
        print(h.figure,fullfile(p.outputDir,'final_state.png'),'-dpng','-r150');
        set(h.title,'String',{titlePrefix; ['Wall-clock elapsed: ' elapsed_string(toc(runClock))]});
    end
end
end

function validate_config(p)
assert(p.nx>=40 && p.ny>=40 && p.steps>=0 && p.steps==floor(p.steps));
assert(p.angleDeg>=0 && p.angleDeg<75,'Use a slope angle between 0 and 75 degrees.');
assert(p.gravity>=0 && p.gravityRamp>0 && p.gravityStart>=0);
assert(p.radius>5 && p.interfaceWidth>0 && p.initialAngleDeg>10 && p.initialAngleDeg<170);
assert(p.sampleEvery>=1 && p.plotEvery>=1 && p.saveEvery>=1);
assert(p.fitMaxDistance>p.fitMinDistance && p.fitOrder>=1 && p.fitOrder<=3);
assert(any(strcmpi(p.relaxationMode,{'source','local'})));
assert(any(strcmpi(p.wallScheme,{'interpolated','halfway'})));
assert(any(strcmpi(p.wettingMode,{'geometric','source_force'})));
assert(p.equilibriumAngleDeg>10 && p.equilibriumAngleDeg<170);
assert(p.ghostReferenceDistance>=sqrt(2));
assert(p.wallHeightLeft>5 && p.wallHeightLeft<p.ny-5);
assert(p.wallHeightLeft-(p.nx-1)*tand(p.angleDeg)>2,...
    'The sloped wall must remain above the bottom outer boundary.');
end

function plan = streaming_plan(ids,neighbors,fluid,slopeSolid,outerSolid,dist,opp,scheme)
map = zeros(size(fluid)); map(ids) = 1:numel(ids);
plan = repmat(struct('interior',[],'sourceRows',[],'near',[],'far',[],...
    'qNear',[],'qFar',[],'backRows',[]),1,9);
for k=1:9
    source = neighbors(:,opp(k)); sourceFluid = fluid(source);
    q = 0.5*ones(numel(ids),1);
    slopeLink = ~sourceFluid & slopeSolid(source) & ~outerSolid(source);
    if strcmpi(scheme,'interpolated')
        q(slopeLink) = dist(ids(slopeLink))./(dist(ids(slopeLink))-dist(source(slopeLink)));
    end
    back = neighbors(:,k); backFluid = fluid(back);
    q(~sourceFluid & q<0.5 & ~backFluid) = 0.5;
    plan(k).interior = find(sourceFluid);
    plan(k).sourceRows = map(source(sourceFluid));
    plan(k).near = find(~sourceFluid & q<0.5);
    plan(k).far = find(~sourceFluid & q>=0.5);
    plan(k).qNear = q(plan(k).near);
    plan(k).qFar = q(plan(k).far);
    plan(k).backRows = map(back(plan(k).near));
end
end

function f = stream_and_reflect(post,plan,opp,conserveMass)
f = zeros(size(post));
massCorrection = zeros(size(post,1),1);
for k=1:9
    a = plan(k); outgoing = opp(k);
    f(a.interior,k) = post(a.sourceRows,k);
    r = a.near; q = a.qNear;
    f(r,k) = 2*q.*post(r,outgoing)+(1-2*q).*post(a.backRows,outgoing);
    r = a.far; q = a.qFar;
    f(r,k) = post(r,outgoing)./(2*q)+(2*q-1)./(2*q).*post(r,k);
    if conserveMass
        boundary = [a.near;a.far];
        massCorrection(boundary) = massCorrection(boundary)+...
            f(boundary,k)-post(boundary,outgoing);
    end
end
% Halfway reflection is globally mass conserving. Subtract the interpolation
% difference from c0=0 locally; this preserves both components of momentum.
% This is an explicit boundary variant requiring its own wetting validation.
if conserveMass, f(:,1) = f(:,1)-massCorrection; end
end

function fraction = gravity_fraction(step,p)
if p.gravity==0, fraction = 0; return; end
q = min(max((step-p.gravityStart)/p.gravityRamp,0),1);
fraction = 0.5*(1-cos(pi*q));
end

function [ax,ay,axFF,ayFF] = all_forces(rho,psi,neighbors,fluid,slopeSolid,ids,p,step)
z = psi(neighbors); mask = double(fluid(neighbors)); base = psi(ids);
if strcmpi(p.wettingMode,'geometric')
    % Full isotropic interaction includes the virtual wall density. No extra
    % Gw force is added: doing both would double-count the wetting condition.
    axFF = -p.G*base.*(4*(z(:,2)-z(:,4))+z(:,6)-z(:,7)+z(:,9)-z(:,8))/36;
    ayFF = -p.G*base.*(4*(z(:,3)-z(:,5))+z(:,6)-z(:,9)+z(:,7)-z(:,8))/36;
    ax = axFF;
    ay = ayFF-rho(ids)*p.gravity*gravity_fraction(step,p);
    return;
end
tx1 = (4*z(:,2)-4*z(:,4)).*mask(:,4).*mask(:,2);
tx2 = (z(:,6)-z(:,7)).*mask(:,7).*mask(:,6);
tx3 = (z(:,9)-z(:,8)).*mask(:,9).*mask(:,8);
[tx2,tx3] = source_stencil(tx1,tx2,tx3);
ty1 = (4*z(:,3)-4*z(:,5)).*mask(:,3).*mask(:,5);
ty2 = (z(:,6)-z(:,9)).*mask(:,9).*mask(:,6);
ty3 = (z(:,7)-z(:,8)).*mask(:,7).*mask(:,8);
[ty2,ty3] = source_stencil(ty1,ty2,ty3);
axFF = -p.G*base.*(tx1+tx2+tx3)/36;
ayFF = -p.G*base.*(ty1+ty2+ty3)/36;
wall = double(slopeSolid(neighbors));
wx = (wall(:,2)-wall(:,4))/3+(wall(:,6)-wall(:,7))/12+(wall(:,9)-wall(:,8))/12;
wy = (wall(:,3)-wall(:,5))/3+(wall(:,6)-wall(:,9))/12+(wall(:,7)-wall(:,8))/12;
ax = axFF-(p.Gw/3)*base.*wx;
ay = ayFF-(p.Gw/3)*base.*wy-rho(ids)*p.gravity*gravity_fraction(step,p);
end

function a = virtual_wall_plan(neighbors,fluid,slopeSolid,outerSolid,N,X,Y,p)
ghosts = unique(neighbors(~fluid(neighbors)));
a.ids = ghosts;
[i,j] = ind2sub(size(fluid),ghosts);
cx = [1 0 -1 0 1 -1 -1 1]; cy = [0 1 0 -1 1 1 -1 -1];
a.neighbors = zeros(numel(ghosts),8);
for k=1:8
    ii = min(max(i+cx(k),1),p.nx); jj = min(max(j+cy(k),1),p.ny);
    a.neighbors(:,k) = sub2ind(size(fluid),ii,jj);
end
a.fluidNeighbors = double(fluid(a.neighbors));
a.count = sum(a.fluidNeighbors,2);
a.slopeRows = find(slopeSolid(ghosts) & ~outerSolid(ghosts));
g = ghosts(a.slopeRows);
distance = p.ghostReferenceDistance-N(g);
tangentShift = distance*cotd(p.equilibriumAngleDeg);
refX = X(g)+distance*sind(p.angleDeg);
refY = Y(g)+distance*cosd(p.angleDeg);
a.plus = sample_plan(refX+tangentShift*cosd(p.angleDeg),...
    refY-tangentShift*sind(p.angleDeg),fluid);
a.minus = sample_plan(refX-tangentShift*cosd(p.angleDeg),...
    refY+tangentShift*sind(p.angleDeg),fluid);
end

function a = sample_plan(x,y,fluid)
[nx,ny] = size(fluid);
inside = x>=0 & x<=nx-1 & y>=0 & y<=ny-1;
i = min(max(floor(x)+1,1),nx-1); j = min(max(floor(y)+1,1),ny-1);
tx = min(max(x-(i-1),0),1); ty = min(max(y-(j-1),0),1);
a.ids = [sub2ind([nx ny],i,j) sub2ind([nx ny],i+1,j) ...
    sub2ind([nx ny],i,j+1) sub2ind([nx ny],i+1,j+1)];
a.weights = [(1-tx).*(1-ty) tx.*(1-ty) (1-tx).*ty tx.*ty];
a.valid = inside & all(fluid(a.ids),2);
end

function extended = virtual_density(rho,a,p)
extended = rho;
if strcmpi(p.wettingMode,'source_force'), return; end
% Neutral zero-gradient extension on the outer box, also the fallback near
% the two slope/box corners. All ghost nodes touch at least one fluid node.
values = sum(rho(a.neighbors).*a.fluidNeighbors,2)./max(a.count,1);
plus = sum(rho(a.plus.ids).*a.plus.weights,2);
minus = sum(rho(a.minus.ids).*a.minus.weights,2);
valid = a.plus.valid & a.minus.valid;
if p.equilibriumAngleDeg<=90, wetting = max(plus,minus);
else, wetting = min(plus,minus); end
rows = a.slopeRows(valid); values(rows) = wetting(valid);
extended(a.ids) = values;
end

function [qx,qy,qxy] = q_components(psi,ids,neighbors,cx,cy,wi,p)
base = psi(ids); xx = zeros(size(base)); yy = xx; xy = xx;
for k=2:9
    delta = psi(neighbors(:,k))-base;
    xx = xx+wi(k)*cx(k)^2*delta;
    yy = yy+wi(k)*cy(k)^2*delta;
    xy = xy+wi(k)*cx(k)*cy(k)*delta;
end
qx = (p.G*p.KQ/6)*base.*xx;
qy = (p.G*p.KQ/6)*base.*yy;
qxy = (p.G*p.KQ/6)*base.*xy;
end

function [pressure,psi] = pr_eos(rho,p,step)
pressure = rho*p.R*p.T./(1-p.b*rho)...
    -p.a*p.betaPR*rho.^2./(1+2*p.b*rho-p.b^2*rho.^2);
arg = 2*(pressure-rho*p.cs2)/(p.G*p.cs2);
if any(~isfinite(arg(:))) || any(arg(:)<0)
    error('Incline:EOS','Invalid real pseudopotential at step %d, min argument %.5g.',step,min(arg(:)));
end
psi = sqrt(arg);
end

function f = equilibrium(d,u,v,cx,cy,w)
f = zeros(numel(d),9); u2 = u.^2+v.^2;
for k=1:9
    cu = cx(k)*u+cy(k)*v;
    f(:,k) = w(k)*d.*(1+3*cu+4.5*cu.^2-1.5*u2);
end
end

function result = ordered_transform(values,matrix)
result = zeros(size(values));
for k=1:9, result = result+values(:,k)*matrix(:,k).'; end
end

function [b,c] = source_stencil(a,b,c)
zero = a==0; b(zero) = 0; c(zero) = 0;
copyC = ~zero & b==0 & c~=0; copyB = ~zero & c==0 & b~=0;
b(copyC) = c(copyC); c(copyB) = b(copyB);
end

function m = observe_drop(rho,u,v,ids,fluid,X,Y,S,N,p)
d = rho(ids); dense = d>(p.rhoV+p.rhoL)/2;
weight = max(d-p.rhoV,0).*dense; total = sum(weight);
m = struct('centerX',NaN,'centerY',NaN,'alongCenter',NaN,'liquidArea',sum(dense),...
    'meanAlongVelocity',NaN,'sReceding',NaN,'sAdvancing',NaN,...
    'thetaReceding',NaN,'thetaAdvancing',NaN,'fitRmsReceding',NaN,...
    'fitRmsAdvancing',NaN,'interface',zeros(2,0),'displacement',NaN);
if total<=0, return; end
m.centerX = sum(weight.*X(ids))/total; m.centerY = sum(weight.*Y(ids))/total;
m.alongCenter = sum(weight.*S(ids))/total;
m.meanAlongVelocity = sum(weight.*(u*cosd(p.angleDeg)-v*sind(p.angleDeg)))/total;
field = rho; field(~fluid) = NaN;
level = (p.rhoV+p.rhoL)/2;
contours = contourc(0:p.nx-1,0:p.ny-1,field.',[level level]);
pos = 1; best = -Inf; arc = zeros(2,0);
while pos<size(contours,2)
    count = round(contours(2,pos)); segment = contours(:,pos+(1:count));
    pos = pos+count+1;
    ns = sind(p.angleDeg)*segment(1,:)+cosd(p.angleDeg)*(segment(2,:)-p.wallHeightLeft);
    if count>=10 && min(ns)<p.fitMinDistance+2
        lengthArc = sum(hypot(diff(segment(1,:)),diff(segment(2,:))));
        if lengthArc>best, arc = segment; best = lengthArc; end
    end
end
if isempty(arc), return; end
m.interface = arc;
along = cosd(p.angleDeg)*arc(1,:)-sind(p.angleDeg)*(arc(2,:)-p.wallHeightLeft);
normal = sind(p.angleDeg)*arc(1,:)+cosd(p.angleDeg)*(arc(2,:)-p.wallHeightLeft);
[~,apex] = max(normal);
first = 1:apex; second = apex:numel(normal);
if mean(along(first))>mean(along(second)), tmp = first; first = second; second = tmp; end
[m.thetaReceding,m.sReceding,m.fitRmsReceding] = fit_contact(along(first),normal(first),1,p);
[m.thetaAdvancing,m.sAdvancing,m.fitRmsAdvancing] = fit_contact(along(second),normal(second),-1,p);
% N is supplied as the geometric distance field for consistency checks.
assert(all(N(ids)>0));
end

function [theta,sContact,rmsFit] = fit_contact(s,n,side,p)
theta = NaN; sContact = NaN; rmsFit = NaN;
keep = n>=p.fitMinDistance & n<=p.fitMaxDistance;
s = s(keep); n = n(keep);
if numel(n)<p.fitOrder+3 || max(n)-min(n)<2, return; end
origin = mean(n); scale = max(std(n),1);
z = (n-origin)/scale;
coeff = polyfit(z,s,p.fitOrder);
sContact = polyval(coeff,-origin/scale);
derivative = polyval(polyder(coeff),-origin/scale)/scale;
theta = acosd(max(-1,min(1,side*derivative/sqrt(1+derivative^2))));
rmsFit = sqrt(mean((polyval(coeff,z)-s).^2));
end

function h = create_plot(rho,solid,X,Y,p)
h.figure = figure('Name','Inclined droplet: actual 30 degree wall','Color','w',...
    'Position',[100 80 1050 800]);
h.axes = subplot(2,2,[1 2],'Parent',h.figure);
h.image = imagesc(h.axes,0:p.nx-1,0:p.ny-1,rho.');
set(h.image,'AlphaData',double(~solid).');
set(h.axes,'YDir','normal','Color',[.55 .55 .55]);
axis(h.axes,'equal'); axis(h.axes,'tight'); hold(h.axes,'on');
caxis(h.axes,[p.rhoV p.rhoL]); %#ok<CAXIS> R2016b compatibility.
colorbar(h.axes); colormap(h.figure,parula);
plot(h.axes,X(:,1),p.wallHeightLeft-X(:,1)*tand(p.angleDeg),'k-','LineWidth',1.5);
h.interface = plot(h.axes,NaN,NaN,'w-','LineWidth',1.4);
h.center = plot(h.axes,NaN,NaN,'r+','MarkerSize',8,'LineWidth',1.5);
h.title = title(h.axes,'');
xlabel(h.axes,'x (lattice units)'); ylabel(h.axes,'y (lattice units)');
ax = subplot(2,2,3,'Parent',h.figure); hold(ax,'on'); grid(ax,'on');
h.receding = plot(ax,NaN,NaN,'-','LineWidth',1.2);
h.advancing = plot(ax,NaN,NaN,'--','LineWidth',1.2);
xlabel(ax,'Lattice step'); ylabel(ax,'Local apparent angle (deg)');
legend(ax,'Receding / uphill','Advancing / downhill','Location','best');
ax = subplot(2,2,4,'Parent',h.figure); grid(ax,'on');
h.displacement = plot(ax,NaN,NaN,'-','LineWidth',1.2);
xlabel(ax,'Lattice step'); ylabel(ax,'Downhill centroid displacement');
assert(all(size(X)==size(Y)));
end

function result = make_result(step,rho,u,v,ids,pressure,solid,N,p,m,rows,reason)
U = zeros(size(rho)); V = U; U(ids) = u; V(ids) = v;
history = array2table(rows,'VariableNames',{'step','gravityFraction','fluidMass',...
    'relativeMassDrift','liquidArea','centerX','centerY','alongCenter',...
    'downhillDisplacement','meanDownhillVelocity','recedingContactS',...
    'advancingContactS','thetaRecedingDeg','thetaAdvancingDeg','maxSpeed','rhoMin',...
    'recedingFitRms','advancingFitRms'});
result = struct('step',step,'rho',rho,'u',U,'v',V,'pressure',pressure,...
    'solid',solid,'normalDistance',N,'parameters',p,'measurement',m,...
    'history',history,'stopReason',reason);
end

function text = elapsed_string(seconds)
t = floor(seconds);
text = sprintf('%02d:%02d:%02d',floor(t/3600),mod(floor(t/60),60),mod(t,60));
end
