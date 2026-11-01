
steps = 40000;
initialPhase = 'liquid';       % 'vapor': 原气泡；'liquid': 气体中的液滴。
Gw = -1;                     % 壁面相互作用强度，不是角度（度）。
relaxationMode = 'source';    % 'source': 与 C++ 相同；'local': 逐点黏度。
nx = 362;
ny = 202;
radius = 40;
centerY = 10;
doPlot = true;
plotEvery = 50;
sampleEvery = 1000;
saveEvery = 10000;
runtimeUpdateSeconds = 1;    % Update the on-screen wall-clock timer every second.
returnPopulations = false;   % 调试时可设 true，返回分布函数。
outputDir = fullfile(fileparts(mfilename('fullpath')), ...
    'matlab_output_single_file');
% outputDir = '';            % 取消注释此行可关闭所有磁盘输出。

% ------------------- 启动计算 -------------------
result = solve_static_wetting( ...
    'Steps',steps,'Nx',nx,'Ny',ny,'Radius',radius,'CenterY',centerY, ...
    'Gw',Gw,'InitialPhase',initialPhase,'RelaxationMode',relaxationMode, ...
    'OutputDir',outputDir,'Plot',doPlot,'PlotEvery',plotEvery, ...
    'SampleEvery',sampleEvery,'SaveEvery',saveEvery, ...
    'ReturnPopulations',returnPopulations,'RuntimeUpdateSeconds',runtimeUpdateSeconds);

disp(result.history);
fprintf('Final apparent liquid-side angle: %.3f deg\n',result.angle.thetaLiquid);
fprintf('Final apparent vapor-side angle:  %.3f deg\n',result.angle.thetaVapor);
fprintf('Circle-fit relative RMS residual: %.4g\n',result.angle.relativeFitError);
fprintf('Total wall-clock run time: %s (%.2f s)\n', ...
    format_runtime(result.elapsedSeconds),result.elapsedSeconds);
if ~isempty(outputDir), fprintf('Outputs: %s\n',outputDir); end

% ------------------- 以下为本文件的局部函数 -------------------

function result = solve_static_wetting(varargin)
%SOLVE_STATIC_WETTING Local solver; all helper functions are in this file.
%
% Original: "final f-s 20-5-2022.cpp", Mohammed Fouad Kamel Ibrahim,
% UNSW Canberra. This is an isothermal D2Q9 MRT pseudopotential model with
% Peng-Robinson EOS, tune correction and Q tensor; it is not a CH/IBM model.
% The default reproduces the SOURCE, including its wall update ordering,
% asymmetric Q-term wrap and global (last-node) shear relaxation rate.
% Arrays use rho(xIndex,yIndex); C++ (i,j) becomes MATLAB (i+1,j+1).
% No optional toolboxes are required. MATLAB R2016b or newer is required.

clockRun = tic; % Actual elapsed runtime, including initialization and plotting.
ip = inputParser;
ip.addParameter('Steps',40000,@(v) integer_at_least(v,0));
ip.addParameter('Nx',362,@(v) integer_at_least(v,12));
ip.addParameter('Ny',202,@(v) integer_at_least(v,12));
ip.addParameter('Gw',0,@(v) isnumeric(v) && isscalar(v) && isfinite(v));
ip.addParameter('Radius',40,@(v) isnumeric(v) && isscalar(v) && v>0);
ip.addParameter('CenterY',10,@(v) isnumeric(v) && isscalar(v) && isfinite(v));
ip.addParameter('InitialPhase','vapor',@(v) any(strcmpi(v,{'vapor','liquid'})));
ip.addParameter('RelaxationMode','source',@(v) any(strcmpi(v,{'source','local'})));
ip.addParameter('Plot',true,@(v) islogical(v) && isscalar(v));
ip.addParameter('PlotEvery',1000,@(v) integer_at_least(v,1));
ip.addParameter('SampleEvery',1000,@(v) integer_at_least(v,1));
ip.addParameter('SaveEvery',10000,@(v) integer_at_least(v,1));
ip.addParameter('Verbose',true,@(v) islogical(v) && isscalar(v));
ip.addParameter('OutputDir',fullfile(fileparts(mfilename('fullpath')),...
    'matlab_output'),@(v) ischar(v) || (isstring(v) && isscalar(v)));
ip.addParameter('WriteTecplot',true,@(v) islogical(v) && isscalar(v));
ip.addParameter('ReturnPopulations',false,@(v) islogical(v) && isscalar(v));
ip.addParameter('WallY',2.5,@(v) isnumeric(v) && isscalar(v) && isfinite(v));
ip.addParameter('FitExclude',3,@(v) isnumeric(v) && isscalar(v) && v>=0);
ip.addParameter('RuntimeUpdateSeconds',1,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>0);
ip.parse(varargin{:});
p = ip.Results;
p.InitialPhase = lower(char(p.InitialPhase));
p.RelaxationMode = lower(char(p.RelaxationMode));
p.OutputDir = char(p.OutputDir); % '' disables ALL disk output (useful for tests).
p.rhoV = 0.34011; p.rhoL = 6.34793;
p.nuV = 0.27462; p.nuL = 0.065;
p.Tc = 0.0729; p.T = 0.85*p.Tc;
p.R = 1; p.omega = 0.344; p.a = 2/49; p.b = 2/21;
p.G = -1; p.cs = 1/sqrt(3); p.cs2 = p.cs^2; p.tune = 0.107; p.KQ = 0.01;
p.betaPR = (1+(0.37464+1.54226*p.omega-0.26992*p.omega^2)*...
    (1-sqrt(p.T/p.Tc)))^2;
nx = p.Nx; ny = p.Ny;
ix = 2:nx-1; iy = 4:ny-2; % C++ i=1..IM-2, j=3..JM-3.
cx = [0 1 0 -1 0 1 -1 -1 1];
cy = [0 0 1 0 -1 1 1 -1 -1];
w = [4/9 repmat(1/9,1,4) repmat(1/36,1,4)];
wi = [0 repmat(1/3,1,4) repmat(1/12,1,4)];
M = [1 1 1 1 1 1 1 1 1;
    -4 -1 -1 -1 -1 2 2 2 2;
    4 -2 -2 -2 -2 1 1 1 1;
    0 1 0 -1 0 1 -1 -1 1;
    0 -2 0 2 0 1 -1 -1 1;
    0 0 1 0 -1 1 1 -1 -1;
    0 0 -2 0 2 1 1 -1 -1;
    0 1 -1 1 -1 0 0 0 0;
    0 0 0 0 0 1 -1 1 -1];
Minv = bsxfun(@rdivide,M.',[9 36 36 6 12 6 12 4 4]);
[X,Y] = ndgrid(0:nx-1,0:ny-1);
inside = hypot(X-(nx-2)/2,Y-p.CenterY)<p.Radius;
inside([1 end],:) = false; inside(:,[1 end]) = false;
if strcmp(p.InitialPhase,'vapor')
    rho = p.rhoL*ones(nx,ny); rho(inside) = p.rhoV;
else
    rho = p.rhoV*ones(nx,ny); rho(inside) = p.rhoL;
end
u = zeros(nx,ny); v = u; fx = u; fy = u; fx1 = u; fy1 = u;
f = bsxfun(@times,rho,reshape(w,1,1,9)); pop = f;
solid = zeros(nx,ny); solid(:,[1:3 ny-1:ny]) = 1;
qx = u; qy = u; qxy = u;
[pressure,psi] = eos(rho,p,0);
hist = zeros(ceil(p.Steps/p.SampleEvery)+1,9); nh = 0;
fig = []; runtimeText = []; runtimeTitlePrefix = '';
lastRuntimeUpdate = -Inf; runtimeLog = zeros(0,2);
if ~isempty(p.OutputDir) && ~exist(p.OutputDir,'dir'), mkdir(p.OutputDir); end
if p.Verbose
    fprintf('UNSW static wetting: %d x %d, %d steps, Gw=%g\n',nx,ny,p.Steps,p.Gw);
    fprintf('Initial inclusion: %s; relaxation: %s; gravity: absent\n',...
        p.InitialPhase,p.RelaxationMode);
end

for step = 0:p.Steps
    if step>0
        % EOS -> Qterm -> COLLID -> STREAM -> MACRO -> fluid/solid force.
        [~,psi] = eos(rho,p,step-1);
        [qx,qy,qxy] = q_tensor(psi,p,ix,iy,cx,cy,wi);
        d = rho(ix,iy); d = d(:);
        ux = u(ix,iy); ux = ux(:); vy = v(ix,iy); vy = vy(:);
        ax = fx(ix,iy); ax = ax(:); ay = fy(ix,iy); ay = ay(:);
        ax1 = fx1(ix,iy); ax1 = ax1(:); ay1 = fy1(ix,iy); ay1 = ay1(:);
        ps = psi(ix,iy); ps = ps(:);
        qxx = qx(ix,iy); qxx = qxx(:);
        qyy = qy(ix,iy); qyy = qyy(:);
        qxyc = qxy(ix,iy); qxyc = qxyc(:);
        nu = p.nuL*(min(max(d,p.rhoV),p.rhoL)-p.rhoV)/(p.rhoL-p.rhoV)...
            +p.nuV*(p.rhoL-min(max(d,p.rhoV),p.rhoL))/(p.rhoL-p.rhoV);
        s = repmat([1 0.8 0.8 1 1.1 1 1.1 1 1],numel(d),1);
        if strcmp(p.RelaxationMode,'source')
            % Original COLLID leaves one GLOBAL s after its first nested loop.
            s(:,8:9) = 1/(3*nu(end)+0.5);
        else
            s(:,8:9) = repmat(1./(3*nu+0.5),1,2);
        end
        speed2 = ux.^2+vy.^2;
        meq = [d d.*(-2+3*speed2) d.*(1-3*speed2) ...
            d.*ux -d.*ux d.*vy -d.*vy d.*(ux.^2-vy.^2) d.*ux.*vy];
        work = ux.*ax+vy.*ay;
        correction = 12*p.tune*(ax1.^2+ay1.^2)./...
            (ps.^2*((1/0.8)-0.5)*p.cs*p.cs);
        forceM = [zeros(size(d)) 6*work+correction -6*work-correction ...
            ax -ax ay -ay 2*(ux.*ax-vy.*ay) ux.*ay+vy.*ax];
        qc = zeros(numel(d),9);
        qc(:,2) = 1.5*s(:,2).*(qxx+qyy);
        qc(:,3) = -1.5*s(:,3).*(qxx+qyy);
        qc(:,8) = -s(:,8).*(qxx-qyy); qc(:,9) = -s(:,8).*qxyc;
        moments = ordered_transform(reshape(f(ix,iy,:),[],9),M);
        moments = moments-s.*(moments-meq)+(1-s/2).*forceM+qc;
        pop(ix,iy,:) = reshape(ordered_transform(moments,Minv),numel(ix),numel(iy),9);
        pop(1,:,:) = pop(nx-1,:,:); pop(nx,:,:) = pop(2,:,:);
        for k = 1:9
            f(ix,2:ny-1,k) = pop(ix-cx(k),(2:ny-1)-cy(k),k);
        end
        % Source updates SOLID pop after streaming, for use in the next step.
        opposite = [1 4 5 2 3 8 9 6 7]; walls = [1:3 ny-1:ny];
        pop(ix,walls,2:9) = f(ix,walls,opposite(2:9));
        rho = zeros(nx,ny); u = rho; v = rho;
        df = zeros(numel(ix),numel(iy)); momentumX = df; momentumY = df;
        for k = 1:9
            df = df+f(ix,iy,k);
            momentumX = momentumX+f(ix,iy,k)*cx(k);
            momentumY = momentumY+f(ix,iy,k)*cy(k);
        end
        df(df<=0) = 1e-30; % Source density floor, not a stability repair.
        rho(ix,iy) = df;
        u(ix,iy) = momentumX./df;
        v(ix,iy) = momentumY./df;
        rho(1,:) = rho(nx-1,:); rho(nx,:) = rho(2,:);
        u(1,:) = u(nx-1,:); u(nx,:) = u(2,:);
        v(1,:) = v(nx-1,:); v(nx,:) = v(2,:);
        [pressure,psi] = eos(rho,p,step);
        [fx,fy,fx1,fy1] = interaction_force(psi,solid,p,ix,iy);
        u(ix,iy) = u(ix,iy)+0.5*fx(ix,iy)./(df+1e-30);
        v(ix,iy) = v(ix,iy)+0.5*fy(ix,iy)./(df+1e-30);
        if any(~isfinite(u(:))) || any(~isfinite(v(:)))
            error('UNSW:Unstable','Non-finite velocity at step %d.',step);
        end
    end
    sampleNow = step==0 || mod(step,p.SampleEvery)==0 || step==p.Steps;
    plotNow = p.Plot && (step==0 || mod(step,p.PlotEvery)==0 || step==p.Steps);
    saveNow = ~isempty(p.OutputDir) && step>0 && ...
        (mod(step,p.SaveEvery)==0 || step==p.Steps);
    if sampleNow || plotNow || saveNow
        angle = measure_static_angle(rho,p.rhoV,p.rhoL,...
            p.InitialPhase,p.WallY,p.FitExclude);
    end
    if sampleNow
        d = rho(ix,iy); speed = hypot(u(ix,iy),v(ix,iy));
        nh = nh+1;
        hist(nh,:) = [step sum(d(:)) min(d(:)) max(d(:)) max(speed(:)) ...
            angle.thetaLiquid angle.thetaVapor angle.relativeFitError toc(clockRun)];
        if p.Verbose
            fprintf('%6d/%d: rho=[%.5g, %.5g], umax=%.3g, thetaL=%.3f deg, %.1fs\n',...
                step,p.Steps,min(d(:)),max(d(:)),max(speed(:)),...
                angle.thetaLiquid,toc(clockRun));
        end
    end
    if plotNow
        if isempty(fig) || ~ishghandle(fig)
            fig = figure('Name','UNSW static contact angle','Color','w');
        end
        set(groot, 'CurrentFigure', fig); % Update without raising the window.
        clf(fig);
        imagesc(1:nx-2,3:ny-3,rho(ix,iy).');
        axis xy equal tight; hold on; colorbar;
        caxis([p.rhoV p.rhoL]); colormap(parula); %#ok<CAXIS> R2016b compatibility.
        if ~isempty(angle.interface)
            plot(angle.interface(1,:),angle.interface(2,:),'w-','LineWidth',1.4);
        end
        plot([1 nx-2],[p.WallY p.WallY],'k--','LineWidth',1);
        xlabel('x (lattice units)'); ylabel('y (lattice units)');
        runtimeTitlePrefix = sprintf('t=%d, %s inclusion, liquid-side angle=%.2f deg',...
            step,p.InitialPhase,angle.thetaLiquid);
        elapsedNow = toc(clockRun);
        runtimeText = title({runtimeTitlePrefix; sprintf(...
            'Wall-clock elapsed: %s (%.1f s)',format_runtime(elapsedNow),elapsedNow)});
        drawnow;
    end
    if saveNow
        snapshot = struct('step',step,'rho',rho,'u',u,'v',v,...
            'psi',psi,'pressure',pressure,'angle',angle,'parameters',p);
        save(fullfile(p.OutputDir,sprintf('snapshot_t%06d.mat',step)),'snapshot');
        if p.WriteTecplot, write_tecplot(rho,psi,pressure,step,p.OutputDir); end
    end
    elapsedNow = toc(clockRun);
    if step==0 || elapsedNow-lastRuntimeUpdate>=p.RuntimeUpdateSeconds || step==p.Steps
        runtimeLog(end+1,:) = [step elapsedNow]; %#ok<AGROW> One row per clock update.
        lastRuntimeUpdate = elapsedNow;
        if p.Plot && ~isempty(runtimeText) && ishghandle(runtimeText)
            set(runtimeText,'String',{runtimeTitlePrefix; sprintf(...
                'Wall-clock elapsed: %s (%.1f s)',format_runtime(elapsedNow),elapsedNow)});
            drawnow limitrate;
        end
    end
end
history = array2table(hist(1:nh,:),'VariableNames',{'step','fluidMass',...
    'rhoMin','rhoMax','maxSpeed','thetaLiquidDeg','thetaVaporDeg',...
    'relativeCircleFitError','elapsedSeconds'});
result = struct('step',p.Steps,'rho',rho,'u',u,'v',v,'psi',psi,...
    'pressure',pressure,'fx',fx,'fy',fy,'fxFluid',fx1,'fyFluid',fy1,...
    'qx',qx,'qy',qy,'qxy',qxy,'angle',angle,'history',history,'parameters',p);
result.elapsedSeconds = toc(clockRun);
result.runtimeHistory = array2table(runtimeLog,...
    'VariableNames',{'step','elapsedSeconds'});
if p.ReturnPopulations, result.f = f; result.pop = pop; end
if ~isempty(p.OutputDir)
    % Keep runtime measurements in memory/on screen, not in output files.
    diskResult = rmfield(result,{'elapsedSeconds','runtimeHistory'});
    diskResult.history.elapsedSeconds = [];
    writetable(diskResult.history,fullfile(p.OutputDir,'history.csv'));
    fileData = struct('result',diskResult);
    save(fullfile(p.OutputDir,'final_result.mat'),'-struct','fileData');
    if p.Plot && ~isempty(fig) && ishghandle(fig)
        savedCaption = get(runtimeText,'String');
        set(runtimeText,'String',runtimeTitlePrefix);
        saveas(fig,fullfile(p.OutputDir,'final_density.png'));
        set(runtimeText,'String',savedCaption);
    end
end
end

function text = format_runtime(seconds)
wholeSeconds = floor(seconds);
text = sprintf('%02d:%02d:%02d',floor(wholeSeconds/3600),...
    mod(floor(wholeSeconds/60),60),mod(wholeSeconds,60));
end

function ok = integer_at_least(v,lowerBound)
ok = isnumeric(v) && isscalar(v) && isfinite(v) && v>=lowerBound && v==floor(v);
end

function transformed = ordered_transform(values,matrix)
% The source tests gradients for EXACT zero. Preserve its addition order;
% BLAS matrix products can change a zero into roundoff and switch a stencil.
transformed = zeros(size(values));
for k = 1:9
    transformed = transformed+values(:,k)*matrix(:,k).';
end
end

function [pressure,psi] = eos(rho,p,step)
pressure = rho*p.R*p.T./(1-p.b*rho)...
    -p.a*p.betaPR*rho.^2./(1+2*p.b*rho-p.b^2*rho.^2);
arg = 2*(pressure-rho*p.cs2)/(p.G*p.cs2);
if any(~isfinite(arg(:))) || any(arg(:)<0)
    error('UNSW:UnstableEOS',...
        'Invalid real pseudopotential at step %d. Min sqrt argument: %.6g.',...
        step,min(arg(:)));
end
psi = sqrt(arg); % Do not hide instability with abs(), real(), or clipping.
end

function [qx,qy,qxy] = q_tensor(psi,p,ix,iy,cx,cy,wi)
base = psi(ix,iy); xx = zeros(size(base)); yy = xx; xy = xx;
for k = 2:9
    xp = ix+cx(k); yp = iy+cy(k);
    % Exact translation of Qterm: xp<1 -> IM-1 in zero-based C++.
    % The leftmost node therefore samples the RIGHT GHOST, not left ghost.
    xp(xp<2) = p.Nx; xp(xp>p.Nx) = 2;
    yp(yp<2) = p.Ny; yp(yp>p.Ny) = 2;
    delta = psi(xp,yp)-base;
    xx = xx+wi(k)*cx(k)^2*delta;
    yy = yy+wi(k)*cy(k)^2*delta;
    xy = xy+wi(k)*cx(k)*cy(k)*delta;
end
qx = zeros(size(psi)); qy = qx; qxy = qx;
qx(ix,iy) = (p.G*p.KQ/6)*base.*xx;
qy(ix,iy) = (p.G*p.KQ/6)*base.*yy;
qxy(ix,iy) = (p.G*p.KQ/6)*base.*xy;
end

function [fx,fy,fx1,fy1] = interaction_force(psi,solid,p,ix,iy)
fluid = 1-solid; base = psi(ix,iy);
tx1 = (4*psi(ix+1,iy)-4*psi(ix-1,iy)).*fluid(ix-1,iy).*fluid(ix+1,iy);
tx2 = (psi(ix+1,iy+1)-psi(ix-1,iy+1)).*fluid(ix-1,iy+1).*fluid(ix+1,iy+1);
tx3 = (psi(ix+1,iy-1)-psi(ix-1,iy-1)).*fluid(ix+1,iy-1).*fluid(ix-1,iy-1);
[tx2,tx3] = source_stencil(tx1,tx2,tx3);
ty1 = (4*psi(ix,iy+1)-4*psi(ix,iy-1)).*fluid(ix,iy+1).*fluid(ix,iy-1);
ty2 = (psi(ix+1,iy+1)-psi(ix+1,iy-1)).*fluid(ix+1,iy-1).*fluid(ix+1,iy+1);
ty3 = (psi(ix-1,iy+1)-psi(ix-1,iy-1)).*fluid(ix-1,iy+1).*fluid(ix-1,iy-1);
[ty2,ty3] = source_stencil(ty1,ty2,ty3);
fx1 = zeros(size(psi)); fy1 = fx1; fx = fx1; fy = fx1;
fx1(ix,iy) = -p.G*base.*(tx1+tx2+tx3)/36;
fy1(ix,iy) = -p.G*base.*(ty1+ty2+ty3)/36;
wallx = (solid(ix+1,iy)-solid(ix-1,iy))/3 ...
    +(solid(ix+1,iy+1)-solid(ix-1,iy+1))/12 ...
    +(solid(ix+1,iy-1)-solid(ix-1,iy-1))/12;
wally = (solid(ix,iy+1)-solid(ix,iy-1))/3 ...
    +(solid(ix+1,iy+1)-solid(ix+1,iy-1))/12 ...
    +(solid(ix-1,iy+1)-solid(ix-1,iy-1))/12;
fx(ix,iy) = fx1(ix,iy)-(p.Gw/3)*base.*wallx;
fy(ix,iy) = fy1(ix,iy)-(p.Gw/3)*base.*wally;
end

function [b,c] = source_stencil(a,b,c)
zero = a==0; b(zero) = 0; c(zero) = 0;
copyC = ~zero & b==0 & c~=0;
copyB = ~zero & c==0 & b~=0;
b(copyC) = c(copyC); c(copyB) = b(copyB);
end

function write_tecplot(rho,psi,pressure,step,folder)
[nx,ny] = size(rho); [x,y] = ndgrid(1:nx-2,1:ny-2);
d = rho(2:nx-1,2:ny-1); ps = psi(2:nx-1,2:ny-1);
pr = pressure(2:nx-1,2:ny-1);
filename = fullfile(folder,sprintf('Tecplot_t%d.dat',step));
[fid,message] = fopen(filename,'w');
if fid<0, error('UNSW:Output','Cannot write %s: %s',filename,message); end
closer = onCleanup(@() fclose(fid));
fprintf(fid,'title = "Fluid"\nvariables = "X","Y","D","psi","P_EOS_2","P_EOS"\n');
fprintf(fid,'zone T="ZONE 1", I=%d, J=%d, K=1, f=point\n',nx-2,ny-2);
% P_EOS_2 is always zero in the original source; keep its output column.
fprintf(fid,'%d %d %.16g %.16g %.16g %.16g\n',...
    [x(:) y(:) d(:) ps(:) zeros(numel(d),1) pr(:)].');
end

function a = measure_static_angle(rho,rhoV,rhoL,phase,wallY,exclude)
%UNSW_MEASURE_CONTACT_ANGLE Apparent static angle from a fitted density contour.
% rho has dimensions [Nx,Ny]. Coordinates are the original zero-based lattice
% coordinates. wallY=2.5 is an EXPLICIT postprocessing convention, not a
% calibrated hydrodynamic wall position. Returns both liquid and vapor angles.
% The fit assumes one approximately circular cap and is not a local dynamic
% contact-angle model. NaN means a usable wall-intersecting cap was not found.
if nargin<5, wallY = 2.5; end
if nargin<6, exclude = 3; end
a = struct('thetaLiquid',NaN,'thetaVapor',NaN,'thetaInclusion',NaN,...
    'center',[NaN NaN],'radius',NaN,'relativeFitError',NaN,...
    'contactX',[NaN NaN],'interface',zeros(2,0),'wallY',wallY,...
    'densityLevel',(rhoV+rhoL)/2,'status','No usable attached cap');
[nx,ny] = size(rho);
c = contourc(1:nx-2,3:ny-3,rho(2:nx-1,4:ny-2).',...
    [a.densityLevel a.densityLevel]);
pos = 1; bestScore = -Inf; arc = zeros(2,0);
while pos<size(c,2)
    n = round(c(2,pos)); segment = c(:,pos+(1:n)); pos = pos+n+1;
    if n>=10 && min(segment(2,:))<=wallY+exclude+2
        score = sum(hypot(diff(segment(1,:)),diff(segment(2,:))));
        if score>bestScore, arc = segment; bestScore = score; end
    end
end
if isempty(arc), return; end
a.interface = arc;
points = arc(:,arc(2,:)>=wallY+exclude & arc(2,:)<ny-6);
if size(points,2)<10, a.status = 'Too few points away from the wall'; return; end
origin = mean(points,2); z = bsxfun(@minus,points,origin);
design = [2*z(1,:).' 2*z(2,:).' ones(size(z,2),1)];
if rank(design)<3, a.status = 'Degenerate circle fit'; return; end
fit = design\sum(z.^2,1).';
radius2 = fit(3)+fit(1)^2+fit(2)^2;
if radius2<=0, a.status = 'Invalid fitted radius'; return; end
a.center = (origin+fit(1:2)).'; a.radius = sqrt(radius2);
distances = hypot(points(1,:)-a.center(1),points(2,:)-a.center(2));
a.relativeFitError = sqrt(mean((distances-a.radius).^2))/a.radius;
ratio = (wallY-a.center(2))/a.radius;
if abs(ratio)>=1, a.status = 'Fitted circle does not cross the chosen wall'; return; end
a.thetaInclusion = acosd(ratio);
if strcmpi(phase,'vapor')
    a.thetaVapor = a.thetaInclusion; a.thetaLiquid = 180-a.thetaVapor;
else
    a.thetaLiquid = a.thetaInclusion; a.thetaVapor = 180-a.thetaLiquid;
end
halfWidth = sqrt(a.radius^2-(wallY-a.center(2))^2);
a.contactX = a.center(1)+[-halfWidth halfWidth];
a.status = 'Circle-fit apparent angle (not a local dynamic angle)';
end

