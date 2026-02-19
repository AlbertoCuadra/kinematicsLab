classdef KinematicsLab < matlab.apps.AppBase
    % Streamline / Pathline / Streakline explorer for 2D velocity fields
    %
    % Educational tool to illustrate core concepts in kinematics and fluid
    % mechanics for steady and unsteady flows. The app provides an interactive
    % environment to compare:
    %
    %   (i)   the instantaneous velocity field (quiver),
    %   (ii)  streamlines at the current time,
    %   (iii) a pathline of a single Lagrangian tracer particle, and
    %   (iv)  a streakline formed by continuously releasing particles from a
    %         fixed location.
    %
    % The goal is to highlight the differences between streamlines, pathlines,
    % and streaklines, and how/when they coincide (e.g., in steady flows).
    %
    % Author: Alberto Cuadra Lara
    %
    % Last update: 19/02/2026

    properties (Access = public)
        UIFigure          matlab.ui.Figure
        LeftPanel         matlab.ui.container.Panel
        Axes              matlab.ui.control.UIAxes

        ProblemDropDown   matlab.ui.control.DropDown

        % Visibility
        VelCheckBox       matlab.ui.control.CheckBox
        StreamCheckBox    matlab.ui.control.CheckBox
        PathCheckBox      matlab.ui.control.CheckBox
        StreakCheckBox    matlab.ui.control.CheckBox

        % Color pickers
        VelColorPicker    matlab.ui.control.ColorPicker
        StreamColorPicker matlab.ui.control.ColorPicker
        PathColorPicker   matlab.ui.control.ColorPicker
        StreakColorPicker matlab.ui.control.ColorPicker

        % Parameters
        V0Field           matlab.ui.control.NumericEditField
        V0Label           matlab.ui.control.Label
        X0Field           matlab.ui.control.NumericEditField
        Y0Field           matlab.ui.control.NumericEditField
        OmegaField        matlab.ui.control.NumericEditField
        OmegaLabel        matlab.ui.control.Label
        X0Label           matlab.ui.control.Label
        Y0Label           matlab.ui.control.Label

        % Buttons
        RunButton         matlab.ui.control.Button
        StopButton        matlab.ui.control.Button
        ClearButton       matlab.ui.control.Button
    end

    properties (Access = private)
        % Simulation settings
        FPS double = 60                  % Fixed simulation frame rate [Hz]
        TMAX double = 12                 % Reference max time (used only for axis sizing in uniform flows)
        DTSTREAK double = 0.2            % Time interval between streak particle releases [s]
        DTPATHLINE double = 1/60         % Time step used for pathline precomputation [s]
        numPointsPathline double = 1000  % Maximum number of points stored for a pathline polyline
        maxStreakPts double = 1200       % Maximum number of streak particles retained
        tickCount uint64 = uint64(0)     % Total number of simulation ticks since start
        numStreamlines int32 = 6         % Number of streamlines to display
        numPointsStreamline double = 200 % Grid resolution for streamline generation via ψ contours
        numContoursStreamline int32 = 6  % Number of streamline contours to display (when using ψ contour method)

        % Timer / simulation state
        Tmr = timer.empty                % Timer object driving the simulation loop
        t double = 0                     % Current simulation time [s]
        isRunning logical = false        % True when the simulation is running

        % Re-entrancy guards
        isResetting logical = false      % Prevents timer callback execution during reset
        isTicking logical = false        % Prevents reset execution during timer callback

        % Graphics handles
        qh                               % Quiver plot handle (velocity field)
        srcDot                           % Source point marker handle
        streamH = gobjects(0)            % Streamline line object handles
        pathLine                         % Pathline polyline handle
        pDot                             % Pathline particle marker handle
        streakLine                       % Streakline polyline handle
        streakDots                       % Streak particle marker handles

        % Grid for quiver visualization
        xg double                        % Quiver grid x-coordinates
        yg double                        % Quiver grid y-coordinates

        % Axis limits
        xMin double                      % Minimum x-axis limit
        xMax double                      % Maximum x-axis limit
        yMin double                      % Minimum y-axis limit
        yMax double                      % Maximum y-axis limit

        % Histories
        pathX double = 0                 % Pathline history: x positions
        pathY double = 0                 % Pathline history: y positions
        pathT double = 0                 % Pathline history: time samples
        pathBox double = [0 0 0 0]       % Pathline bounding box (used for validity / clipping checks)
        streakX double = 0               % Streakline history: x positions
        streakY double = 0               % Streakline history: y positions
        lastReleaseT double = 0          % Time of last streak particle release

        % UI styling
        fontsizeLabels = 18              % Font size for axis labels and UI controls

        % View / camera behaviour
        qNx double = 14                  % Number of quiver grid points in x
        qNy double = 10                  % Number of quiver grid points in y
        viewDirty logical = false        % Indicates pending view update (axis/grid change)
        lastViewApplyT double = 0        % Time when the view was last applied
        viewApplyMinPeriod double = 1/20 % Minimum interval between view updates [s] (max 20 Hz)

        % Axis listeners
        xLimListener event.listener = event.listener.empty  % Listener for x-axis limit changes
        yLimListener event.listener = event.listener.empty  % Listener for y-axis limit changes

        % Flow 
        Flows struct = struct( ...
            'Label', {}, ...          % Dropdown label
            'Velocity', {}, ...       % @(app,x,y,t) -> [u,v]
            'AxisLimits', {}, ...     % @(app) -> [xMin xMax yMin yMax]
            'Streamlines', {}, ...    % @(app,tFix,n) -> cell arrays {x_i},{y_i} (optional)
            'UsesV0', {}, ...
            'UsesOmega', {}, ...
            'IsUnsteady', {} );
        flowIndex uint16 = 1;         % active flow index into Flows
    end


    methods (Access = public)
        function app = KinematicsLab
            % Constructor
            createComponents(app);
            startup(app);
        end

        function delete(app)
            % Destructor
            stopTimer(app);
            if isvalid(app.UIFigure)
                delete(app.UIFigure);
            end

        end

        function addFlow(app, flow)
            % Add a new flow definition to the Flow library
            %
            % Args:
            %     flow (struct): Struct defining the flow, with at least the required fields
            %
            % Note:
            %     * The required fields are Label and Velocity, which define the dropdown label and the velocity function, respectively.
            %     * The velocity function should have the signature @(app,x,y,t) -> [u,v], where (x,y) are coordinates and t is time.
            %     * Optional fields include AxisLimits (for automatic axis sizing), Streamlines (for custom streamline generation), and flags indicating whether the flow uses V0, omega, or is unsteady.
            %
            % Example:
            %     app.addFlow(struct( 
            %         'Label', 'Steady uniform: u=V0, v=0', ...
            %         'UsesV0', true, 'UsesOmega', false, 'IsUnsteady', false, ...
            %         'Velocity', @(app,x,y,t) deal(app.V0Field.Value + 0*x, 0 + 0*y), ...
            %         'AxisLimits',  @(app) axisLimitsSteadyUniform(app), ...
            %         'Streamlines', @(app,tFix,n) streamlinesUniform(app,tFix,n) ...
            %     ));

            arguments
                app
                flow struct
            end

            mustHave = {'Label','Velocity'};
            for k = 1:numel(mustHave)
                if ~isfield(flow, mustHave{k})
                    error('addFlow:MissingField', 'Flow must define "%s".', mustHave{k});
                end
            end

            % Defaults
            if ~isfield(flow,'UsesV0'),     flow.UsesV0 = false; end
            if ~isfield(flow,'UsesOmega'),  flow.UsesOmega = false; end
            if ~isfield(flow,'IsUnsteady'), flow.IsUnsteady = false; end
            if ~isfield(flow,'Streamlines'), flow.Streamlines = []; end
            if ~isfield(flow,'AxisLimits')
                flow.AxisLimits = @(app) defaultAxisLimits(app); % fallback
            end

            app.Flows(end+1) = flow;

            % Keep dropdown synced
            if ~isempty(app.ProblemDropDown) && isvalid(app.ProblemDropDown)
                app.ProblemDropDown.Items = {app.Flows.Label};
            end

        end

    end

    methods (Access = private)

        function startup(app)
            % Initialize defaults and reset the app.
            %
            % Args:
            %     app (KinematicsLab): App instance

            app.Flows = struct('Label',{},'Velocity',{},'AxisLimits',{},'Streamlines',{}, ...
                            'UsesV0',{},'UsesOmega',{},'IsUnsteady',{});

            % 1) Unsteady uniform: u=V0, v=V0 sin(ωt)
            app.addFlow(struct( ...
                'Label', 'Unsteady uniform: u=V0, v=V0 sin(ωt)', ...
                'UsesV0', true, 'UsesOmega', true, 'IsUnsteady', true, ...
                'Velocity', @(app,x,y,t) deal( ...
                    app.V0Field.Value + 0*x, ...
                    app.V0Field.Value * sin(app.OmegaField.Value*t) + 0*y), ...
                'AxisLimits',  @(app) axisLimitsUnsteadyUniform(app), ...
                'Streamlines', @(app,tFix,n) streamlinesUniform(app,tFix,n) ...
            ));

            % 2) Steady uniform: u=V0, v=0
            app.addFlow(struct( ...
                'Label', 'Steady uniform: u=V0, v=0', ...
                'UsesV0', true, 'UsesOmega', false, 'IsUnsteady', false, ...
                'Velocity', @(app,x,y,t) deal( ...
                    app.V0Field.Value + 0*x, ...
                    0 + 0*y), ...
                'AxisLimits',  @(app) axisLimitsSteadyUniform(app), ...
                'Streamlines', @(app,tFix,n) streamlinesUniform(app,tFix,n) ...
            ));

            % 3) Steady saddle: u=x, v=-y
            app.addFlow(struct( ...
                'Label', 'Steady saddle: u=x, v=-y', ...
                'UsesV0', false, 'UsesOmega', false, 'IsUnsteady', false, ...
                'Velocity', @(app,x,y,t) deal(x, -y), ...
                'Streamlines', @(app,tFix,n) streamlinesSaddle(app,n) ...
            ));

            % 4) Steady source: u=x, v=y
            app.addFlow(struct( ...
                'Label', 'Steady source: u=x, v=y', ...
                'UsesV0', false, 'UsesOmega', false, 'IsUnsteady', false, ...
                'Velocity', @(app,x,y,t) deal(x, y), ...
                'Streamlines', @(app,tFix,n) streamlinesSource(app,n) ...
            ));

            % 5) Steady rotation: u=y, v=-x
            app.addFlow(struct( ...
                'Label', 'Steady rotation: u=y, v=-x', ...
                'UsesV0', false, 'UsesOmega', false, 'IsUnsteady', false, ...
                'Velocity', @(app,x,y,t) deal(y, -x), ...
                'Streamlines', @(app,tFix,n) streamlinesRotation(app,n) ...
            ));

            % 6) Steady shear: u=V0, v=V0 x
            app.addFlow(struct( ...
                'Label', 'Steady shear: u=V0, v=V0 x', ...
                'UsesV0', true, 'UsesOmega', false, 'IsUnsteady', false, ...
                'Velocity', @(app,x,y,t) deal( ...
                    app.V0Field.Value + 0*x, ...
                    app.V0Field.Value .* x), ...
                'Streamlines', @(app,tFix,n) streamlinesShear(app,n) ...
            ));

            % 7) Steady cellular: u=cos(ωx), v=sin(ωy)
            app.addFlow(struct( ...
                'Label', 'Steady cellular: u=cos(ωx), v=sin(ωy)', ...
                'UsesV0', false, 'UsesOmega', true, 'IsUnsteady', false, ...
                'Velocity', @(app,x,y,t) deal( ...
                    cos(app.OmegaField.Value .* x), ...
                    sin(app.OmegaField.Value .* y)), ...
                'AxisLimits',  @(app) axisLimitsCellular(app), ...
                'Streamlines', @(app,tFix,n) streamlinesCellular(app,tFix,n) ...
            ));

            % 8) Steady Taylor–Green vortex (periodic cellular vortices)
            app.addFlow(struct( ...
                'Label', 'Steady Taylor–Green: u=sin(ωx)cos(ωy), v=-cos(ωx)sin(ωy)', ...
                'UsesV0', false, 'UsesOmega', true, 'IsUnsteady', false, ...
                'Velocity', @(app,x,y,t) taylorGreenVel(app,x,y,t), ...
                'AxisLimits', @(app) axisLimitsTaylorGreen(app), ...
                'Streamlines', @(app,tFix,n) streamlinesTaylorGreen(app,tFix,n) ... % contours of ψ
            ));

            % 9) Unsteady rotating saddle (saddle rotated by angle ωt)
            app.addFlow(struct( ...
                'Label', 'Unsteady rotating saddle (ωt-rotated u=x, v=-y)', ...
                'UsesV0', false, 'UsesOmega', true, 'IsUnsteady', true, ...
                'Velocity', @(app,x,y,t) rotatingSaddleVel(app,x,y,t), ...
                'Streamlines', @(app,tFix,n) streamlinesRotatingSaddle(app,tFix,n) ... % contours of ψ
            ));

            % 10) Steady spiral sink (focus): u=-a x - V0 y, v=V0 x - a y
            app.addFlow(struct( ...
                'Label', 'Steady spiral sink (focus): u=-a x - V0 y, v=V0 x - a y', ...
                'UsesV0', true, 'UsesOmega', false, 'IsUnsteady', false, ...
                'Velocity', @(app,x,y,t) spiralSinkVel(app,x,y,t) ...
            ));

            % Defaults
            app.ProblemDropDown.Items = {app.Flows.Label};
            app.ProblemDropDown.Value = app.Flows(1).Label;
            app.flowIndex = uint16(1);
            
            app.VelCheckBox.Value    = true;
            app.StreamCheckBox.Value = true;
            app.PathCheckBox.Value   = true;
            app.StreakCheckBox.Value = true;

            app.V0Field.Value    = 1.0;
            app.X0Field.Value    = 0.0;
            app.Y0Field.Value    = 0.0;
            app.OmegaField.Value = 2;

            updateOmegaEnable(app);
            updateV0Enable(app);
            resetAll(app);

            % NESTED FUNCTIONS
            function lim = axisLimitsUnsteadyUniform(app)
                V0 = app.V0Field.Value;
                om = app.OmegaField.Value;
                x0 = app.X0Field.Value;
                y0 = app.Y0Field.Value;

                Vabs = abs(V0);
                xMin = x0 - Vabs * max(app.TMAX, 0) - 1;
                xMax = x0 + Vabs * max(app.TMAX, 0) + 1;

                if abs(om) > 1e-12
                    A = 2 * Vabs / abs(om);
                    if ~isfinite(A), A = 1; end
                else
                    A = 1;
                end

                yMin = y0 - 1.2 * A - 1;
                yMax = y0 + 1.2 * A + 1;

                lim = [xMin xMax yMin yMax];
            end

            function lim = axisLimitsSteadyUniform(app)
                V0 = app.V0Field.Value;
                x0 = app.X0Field.Value;
                y0 = app.Y0Field.Value;

                Vabs = abs(V0);
                xMin = x0 - Vabs * max(app.TMAX, 0) - 1;
                xMax = x0 + Vabs * max(app.TMAX, 0) + 1;

                yMin = y0 - 4;
                yMax = y0 + 4;

                lim = [xMin xMax yMin yMax];
            end

            function lim = axisLimitsCellular(app)
                om = app.OmegaField.Value;
                x0 = app.X0Field.Value;
                y0 = app.Y0Field.Value;

                K = max(abs(om), 0.5);
                L = 2*pi / K;

                lim = [x0 - L, x0 + L, y0 - L, y0 + L];
            end

            function lim = axisLimitsTaylorGreen(app)
                om = app.OmegaField.Value;
                x0 = app.X0Field.Value;
                y0 = app.Y0Field.Value;

                k = max(abs(om), 0.5);
                L = 2*pi / k;          % one spatial period
                half = 0.55 * L;       % show ~one period with a bit of margin

                lim = [x0 - half, x0 + half, y0 - half, y0 + half];
            end

            function [XS, YS] = streamlinesUniform(app, tFix, n)
                x = linspace(app.xMin, app.xMax, 200);

                [u0, v0] = app.velocity(0, 0, tFix);
                slope = 0;
                if abs(u0) > 1e-12
                    slope = v0 / u0;
                end

                xRef = x(1);
                C = linspace(app.yMin, app.yMax, n);

                XS = cell(n,1);
                YS = cell(n,1);
                for i = 1:n
                    XS{i} = x;
                    YS{i} = C(i) + slope * (x - xRef);
                end
            end

            function [XS, YS] = streamlinesSaddle(app, n)
                x = linspace(app.xMin, app.xMax, 200);

                XS = cell(n,1);
                YS = cell(n,1);

                if n >= 1
                    XS{1} = x; YS{1} = 0*x; % y=0
                end
                if n >= 2
                    yAxis = linspace(app.yMin, app.yMax, 800);
                    XS{2} = 0*yAxis; YS{2} = yAxis; % x=0
                end

                m = max(n-2, 0);
                if m > 0
                    Cmax = 0.35 * (max(abs([app.xMin, app.xMax])) * max(abs([app.yMin, app.yMax])) + 1);
                    Cvals = linspace(-Cmax, Cmax, m);
                    Cvals(abs(Cvals) < 1e-9) = 0.15 * Cmax;

                    xx = x;
                    epsx = 1e-6 * max(1, max(abs(x)));
                    xx(abs(xx) < epsx) = epsx;

                    for k = 1:m
                        i = k + 2;
                        XS{i} = xx;
                        YS{i} = Cvals(k) ./ xx; % x*y=C
                    end
                end

                for i = 1:n
                    if isempty(XS{i})
                        XS{i} = nan; YS{i} = nan;
                    end
                end
            end

            function [XS, YS] = streamlinesSource(app, n)
                theta = linspace(0, pi, n);
                R = 1.5 * max(app.xMax - app.xMin, app.yMax - app.yMin);
                r = linspace(-R, R, 1200);

                XS = cell(n,1);
                YS = cell(n,1);
                for i = 1:n
                    XS{i} = r * cos(theta(i));
                    YS{i} = r * sin(theta(i));
                end
            end

            function [XS, YS] = streamlinesRotation(app, n)
                rMax = 0.95 * min(max(abs([app.xMin, app.xMax])), max(abs([app.yMin, app.yMax])));
                rMin = max(0.15 * rMax, 0.2);
                rVals = linspace(rMin, rMax, n);
                th = linspace(0, 2*pi, 600);

                XS = cell(n,1);
                YS = cell(n,1);
                for i = 1:n
                    XS{i} = rVals(i) * cos(th);
                    YS{i} = rVals(i) * sin(th);
                end
            end

            function [XS, YS] = streamlinesShear(app, n)
                x = linspace(app.xMin, app.xMax, 200);
                C = linspace(app.yMin, app.yMax, n);

                XS = cell(n,1);
                YS = cell(n,1);
                for i = 1:n
                    XS{i} = x;
                    YS{i} = 0.5 * x.^2 + C(i); % y = 0.5 x^2 + C
                end
            end

            function [XS, YS] = streamlinesCellular(app, tFix, n)
                om = app.OmegaField.Value;
            
                % ω -> 0 limit: u=1, v=0 => streamlines are y = const
                if abs(om) < 1e-12
                    psiFcn = @(app,X,Y,tFix) Y;
                    [XS, YS] = streamlinesFromPsi(app, psiFcn, tFix, n);
                    return
                end
            
                % "Psi" here is a FIRST INTEGRAL (invariant), not a true streamfunction.
                % Streamlines satisfy: tan(ωy/2) / (sec(ωx)+tan(ωx)) = const
                % We contour atan(ratio) to compress dynamic range and avoid blow-ups.
                psiFcn = @(app,X,Y,tFix) cellularInvariant(app, X, Y, om);
            
                [XS, YS] = streamlinesFromPsi(app, psiFcn, tFix, n);
            
                function psi = cellularInvariant(app, X, Y, om) %#ok<INUSL>
                    wx = om .* X;
                    wy = om .* Y;
            
                    denom = (1 ./ cos(wx)) + tan(wx);
                    ratio = tan(wy/2) ./ denom; % invariant C
            
                    ratio(~isfinite(ratio)) = NaN;
            
                    % Compress range: preserves ordering/sign (ratio = tan(psi))
                    psi = atan(ratio);
                end
            end

            % Streamlines via ψ contours
            function [XS, YS] = streamlinesFromPsi(app, psiFcn, tFix, n)
                % Streamlines as contours of streamfunction psi

                % Definitions
                Nx = app.numPointsStreamline;
                Ny = app.numPointsStreamline;
                xvec = linspace(app.xMin, app.xMax, Nx);
                yvec = linspace(app.yMin, app.yMax, Ny);
                [X,Y] = meshgrid(xvec, yvec);
                numContoursStreamline = app.numContoursStreamline;

                psi = psiFcn(app, X, Y, tFix);
                psi(~isfinite(psi)) = NaN;

                XS = cell(n,1);
                YS = cell(n,1);

                p = psi(:);
                p = p(isfinite(p));
                if isempty(p)
                    for i = 1:n, XS{i} = nan; YS{i} = nan; end
                    return
                end

                % Robust magnitude range (avoid corner outliers)
                absMax = prctile(abs(p), 90);
                if ~isfinite(absMax) || absMax < 1e-12
                    for i = 1:n, XS{i} = nan; YS{i} = nan; end
                    return
                end

                % Ensures a separatrix level even when n is even
                levels = linspace(-absMax, absMax, n);
                levels(round((n+1)/2)) = 0; 

                for i = 1:n
                    lev = levels(i);
                    C = contourc(xvec, yvec, psi, [lev lev]);

                    [xs, ys] = collectTopSegments(C, numContoursStreamline);
                    XS{i} = xs;
                    YS{i} = ys;
                end

                % NESTED FUNCTIONS
                function [xs, ys] = collectTopSegments(C, numContoursStreamline)
                    xs = nan; ys = nan;
                    if isempty(C), return; end

                    segs = {};
                    lens = [];

                    idx = 1;
                    while idx < size(C,2)
                        npts = C(2,idx);
                        if idx + npts > size(C,2), break; end
                        seg = C(:, idx+1:idx+npts);
                        segs{end+1} = seg;
                        lens(end+1) = npts;
                        idx = idx + npts + 1;
                    end

                    if isempty(segs), return; end

                    % Keep largest segments
                    [~, ord] = sort(lens, 'descend');
                    ord = ord(1:min(numContoursStreamline, numel(ord)));

                    xx = [];
                    yy = [];
                    for k = 1:numel(ord)
                        s = segs{ord(k)};
                        xx = [xx, s(1,:), NaN];
                        yy = [yy, s(2,:), NaN];
                    end

                    % Remove trailing NaN
                    if ~isempty(xx)
                        xx(end) = [];
                        yy(end) = [];
                    end

                    xs = xx;
                    ys = yy;
                end

            end

            function [u, v] = taylorGreenVel(app, x, y, t)
                k = max(abs(app.OmegaField.Value), 0.5);
                u =  sin(k.*x) .* cos(k.*y);
                v = -cos(k.*x) .* sin(k.*y);
            end

            function [XS, YS] = streamlinesTaylorGreen(app, tFix, n)
                k = max(abs(app.OmegaField.Value), 0.5);
                psiFcn = @(app, X, Y, tFix) sin(k.*X) .* sin(k.*Y);
                [XS, YS] = streamlinesFromPsi(app, psiFcn, tFix, n);
            end

            function [u, v] = rotatingSaddleVel(app, x, y, t)
                th = app.OmegaField.Value .* t;
                c = cos(th); s = sin(th);

                % Rotate coordinates: [X;Y] = R^T [x;y]
                X = c.*x + s.*y;
                Y = -s.*x + c.*y;

                % Saddle in rotated coords: U=X, V=-Y; transform back: [u;v]=R[U;V]
                u = c.*X + s.*Y;
                v = s.*X - c.*Y;
            end

            function [XS, YS] = streamlinesRotatingSaddle(app, tFix, n)
                th = app.OmegaField.Value .* tFix;
                c = cos(th); s = sin(th);

                psiFcn = @(app, Xg, Yg, tFix) ( (c.*Xg + s.*Yg) .* (-s.*Xg + c.*Yg) );

                [XS, YS] = streamlinesFromPsi(app, psiFcn, tFix, n);
            end

            function [u, v] = spiralSinkVel(app, x, y, t)
                V0 = app.V0Field.Value;

                a = 0.25; % Decay rate
                u = -a.*x - V0.*y;
                v =  V0.*x - a.*y;
            end

        end

        function lim = defaultAxisLimits(app)
            % Default axis limits (used when flow doesn't specify its own limits)
            [x0,y0] = getX0Y0(app);
            lim = [x0 - 4, x0 + 4, y0 - 4, y0 + 4];
        end

        %% -------------------- UI --------------------
        function createComponents(app)
            % Create UI layout and wire callbacks.
            %
            % Args:
            %     app (KinematicsLab): App instance

            app.UIFigure = uifigure('Name', 'Kinematics', 'Color', [0.94, 0.94, 0.94]);
            app.UIFigure.Position(3:4) = [900, 400];
            app.UIFigure.CloseRequestFcn = @(~, ~) onClose(app);

            gl = uigridlayout(app.UIFigure, [1, 2]);
            gl.ColumnWidth = {240, '1x'};
            gl.RowHeight = {'1x'};
            gl.Padding = [10, 10, 10, 10];
            gl.ColumnSpacing = 10;

            app.LeftPanel = uipanel(gl, ...
                'Title', 'Setup', ...
                'BackgroundColor', [0.97, 0.97, 0.97], ...
                'BorderType', 'none', ...
                'FontSize', app.fontsizeLabels - 4);

            app.Axes = uiaxes(gl);
            configureAxesAppearance(app);

            % Left-panel grid
            lp = uigridlayout(app.LeftPanel, [12, 3]);
            lp.Scrollable = 'on';
            lp.RowHeight = {22, 22, 22, 22, 22, 22, 22, 22, 22, 12, 30, 30};
            lp.ColumnWidth = {100, '1x', 50};
            lp.Padding = [10, 10, 10, 10];
            lp.RowSpacing = 6;
            lp.ColumnSpacing = 8;

            % Problem row
            lbl = uilabel(lp, 'Text', 'Problem', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            lbl.Layout.Row = 1;
            lbl.Layout.Column = 1;

            app.ProblemDropDown = uidropdown(lp, 'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.ProblemDropDown.Layout.Row = 1;
            app.ProblemDropDown.Layout.Column = [2, 3];

            % Velocity
            app.VelCheckBox = uicheckbox(lp, 'Text', 'Velocity field', 'ValueChangedFcn', @(~, ~) toggleVisibility(app), ...
                'FontSize', app.fontsizeLabels - 4);
            app.VelCheckBox.Layout.Row = 2;
            app.VelCheckBox.Layout.Column = [1, 2];

            app.VelColorPicker = uicolorpicker(lp, 'Value', [204, 204, 204] / 255, 'ValueChangedFcn', @(~, ~) colorChanged(app));
            app.VelColorPicker.Layout.Row = 2;
            app.VelColorPicker.Layout.Column = 3;

            % Streamline
            app.StreamCheckBox = uicheckbox(lp, 'Text', 'Streamline', 'ValueChangedFcn', @(~, ~) toggleVisibility(app), ...
                'FontSize', app.fontsizeLabels - 4);
            app.StreamCheckBox.Layout.Row = 3;
            app.StreamCheckBox.Layout.Column = [1, 2];

            app.StreamColorPicker = uicolorpicker(lp, 'Value', [0.91, 0.51, 0.49], 'ValueChangedFcn', @(~, ~) colorChanged(app));
            app.StreamColorPicker.Layout.Row = 3;
            app.StreamColorPicker.Layout.Column = 3;

            % Pathline
            app.PathCheckBox = uicheckbox(lp, 'Text', 'Pathline', 'ValueChangedFcn', @(~, ~) toggleVisibility(app), ...
                'FontSize', app.fontsizeLabels - 4);
            app.PathCheckBox.Layout.Row = 4;
            app.PathCheckBox.Layout.Column = [1, 2];

            app.PathColorPicker = uicolorpicker(lp, 'Value', [89, 156, 155] / 255, 'ValueChangedFcn', @(~, ~) colorChanged(app));
            app.PathColorPicker.Layout.Row = 4;
            app.PathColorPicker.Layout.Column = 3;

            % Streakline
            app.StreakCheckBox = uicheckbox(lp, 'Text', 'Streakline', 'ValueChangedFcn', @(~, ~) toggleVisibility(app), ...
                'FontSize', app.fontsizeLabels - 4);
            app.StreakCheckBox.Layout.Row = 5;
            app.StreakCheckBox.Layout.Column = [1, 2];

            app.StreakColorPicker = uicolorpicker(lp, 'Value', [152, 209, 208] / 255, 'ValueChangedFcn', @(~, ~) colorChanged(app));
            app.StreakColorPicker.Layout.Row = 5;
            app.StreakColorPicker.Layout.Column = 3;

            % V0
            app.V0Label = uilabel(lp, 'Text', '$V_0$ [m/s]', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            app.V0Label.Interpreter = 'latex';
            app.V0Label.Layout.Row = 6;
            app.V0Label.Layout.Column = 1;

            app.V0Field = uieditfield(lp, 'numeric', 'Value', 1.0, 'Limits', [-Inf, Inf], ...
                'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.V0Field.Layout.Row = 6;
            app.V0Field.Layout.Column = [2, 3];

            % x0
            app.X0Label = uilabel(lp, 'Text', '$x_0$ [m]', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            app.X0Label.Interpreter = 'latex';
            app.X0Label.Layout.Row = 7;
            app.X0Label.Layout.Column = 1;

            app.X0Field = uieditfield(lp, 'numeric', 'Value', 0.0, 'Limits', [-Inf, Inf], ...
                'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.X0Field.Layout.Row = 7;
            app.X0Field.Layout.Column = [2, 3];

            % y0
            app.Y0Label = uilabel(lp, 'Text', '$y_0$ [m]', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            app.Y0Label.Interpreter = 'latex';
            app.Y0Label.Layout.Row = 8;
            app.Y0Label.Layout.Column = 1;

            app.Y0Field = uieditfield(lp, 'numeric', 'Value', 0.0, 'Limits', [-Inf, Inf], ...
                'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.Y0Field.Layout.Row = 8;
            app.Y0Field.Layout.Column = [2, 3];

            % Omega
            app.OmegaLabel = uilabel(lp, 'Text', '$\omega$ [1/s]', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            app.OmegaLabel.Interpreter = 'latex';
            app.OmegaLabel.Layout.Row = 9;
            app.OmegaLabel.Layout.Column = 1;

            app.OmegaField = uieditfield(lp, 'numeric', 'Value', 2, 'Limits', [-Inf, Inf], ...
                'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.OmegaField.Layout.Row = 9;
            app.OmegaField.Layout.Column = [2, 3];

            % Buttons
            btnGrid = uigridlayout(lp, [1, 2]);
            btnGrid.Layout.Row = 11;
            btnGrid.Layout.Column = [1, 3];
            btnGrid.ColumnWidth = {'1x', '1x'};
            btnGrid.RowHeight = {30};
            btnGrid.Padding = [0, 0, 0, 0];
            btnGrid.ColumnSpacing = 10;

            app.RunButton = uibutton(btnGrid, 'Text', 'run', 'ButtonPushedFcn', @(~, ~) onRun(app), 'FontSize', app.fontsizeLabels - 4);
            app.StopButton = uibutton(btnGrid, 'Text', 'stop', 'ButtonPushedFcn', @(~, ~) onStop(app), 'FontSize', app.fontsizeLabels - 4);

            app.ClearButton = uibutton(lp, 'Text', 'clear', 'ButtonPushedFcn', @(~, ~) resetAll(app), 'FontSize', app.fontsizeLabels - 4);
            app.ClearButton.Layout.Row = 12;
            app.ClearButton.Layout.Column = [1, 3];
        end

        function configureAxesAppearance(app)
            % Configure axes aesthetics and LaTeX labels.
            %
            % Args:
            %     app (KinematicsLab): App instance

            app.Axes.Box = 'off';
            app.Axes.FontSize = app.fontsizeLabels + 2;
            app.Axes.XLabel.FontSize = app.fontsizeLabels;
            app.Axes.YLabel.FontSize = app.fontsizeLabels;
            app.Axes.Title.FontSize = app.fontsizeLabels + 4;
            app.Axes.LineWidth = 1.8;

            app.Axes.TickLabelInterpreter = 'latex';
            app.Axes.XLabel.Interpreter = 'latex';
            app.Axes.YLabel.Interpreter = 'latex';
            app.Axes.Title.Interpreter = 'latex';

            app.Axes.XLabel.String = '$x$ [m]';
            app.Axes.YLabel.String = '$y$ [m]';

            app.Axes.XGrid = 'off';
            app.Axes.YGrid = 'off';
            app.Axes.XMinorGrid = 'off';
            app.Axes.YMinorGrid = 'off';

            updateTitle(app);
        end

        function onClose(app)
            stopTimer(app);
            delete(app);
        end

        function paramChanged(app)
            % Callback for parameter changes (problem / V0 / x0 / y0 / omega).
            % Stops simulation, resets all graphics & histories, keeps running if it was running.
            wasRunning = app.isRunning;

            app.flowIndex = uint16(find(strcmp(app.ProblemDropDown.Value, {app.Flows.Label}), 1, 'first'));

            if isempty(app.flowIndex),
                app.flowIndex = 1;
            end

            updateOmegaEnable(app);
            updateV0Enable(app);
            resetAll(app);

            if wasRunning
                startTimer(app);
            end

        end

        function colorChanged(app)
            updateColors(app);
        end

        %% ------------------ Enable/Disable parameters ------------------
        function updateOmegaEnable(app)
            % Enable/disable omega input depending on the selected problem

            if problemUsesOmega(app)
                app.OmegaField.Enable = 'on';
                return
            end

            app.OmegaField.Enable = 'off';
        end

        function updateV0Enable(app)
            % Enable/disable V0 input depending on the selected problem

            if problemUsesV0(app)
                app.V0Field.Enable = 'on';
                return
            end
            app.V0Field.Enable = 'off';
        end

        function tf = problemUsesOmega(app)
            % Return true if the current problem definition uses omega
            %
            % Returns:
            %     tf (logical): True if omega is used by the selected velocity field
            tf = app.Flows(app.flowIndex).UsesOmega;
        end

        function tf = problemUsesV0(app)
            % Return true if the current problem definition uses V0
            %
            % Returns:
            %     tf (logical): True if V0 is used by the selected velocity field
            tf = app.Flows(app.flowIndex).UsesV0;
        end

        function tf = isUnsteady(app)
            % Return true if the selected velocity field is time-dependent
            tf = app.Flows(app.flowIndex).IsUnsteady;
        end

        function [x0, y0] = getX0Y0(app)
            % Get initial position vector
            %
            % Returns:
            %     x0 (double): Initial x-position [m]
            %     y0 (double): Initial y-position [m]
            x0 = app.X0Field.Value;
            y0 = app.Y0Field.Value;
        end

        %% -------------------- Flow definition --------------------
        function [u, v] = velocity(app, x, y, t)
            % Evaluate the velocity field at (x, y, t)
            %
            % Args:
            %     x (double): x-position
            %     y (double): y-position
            %     t (double): time [s]
            %
            % Returns:
            %     u (double): x-velocity component at (x, y, t)
            %     v (double): y-velocity component at (x, y, t)

            F = app.Flows(app.flowIndex);
            [u, v] = F.Velocity(app, x, y, t);
        end


        %% -------------------- RK4 helpers --------------------
        function [xn, yn] = rk4(app, x, y, t, dt)
            % Advance one particle with classic RK4
            %
            % Args:
            %     x, y (double): current particle position
            %     t (double): current time
            %     dt (double): time step
            %
            % Returns:
            %     xn, yn (double): updated particle position after dt

            [k1x, k1y] = app.velocity(x, y, t);
            [k2x, k2y] = app.velocity(x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, t + 0.5 * dt);
            [k3x, k3y] = app.velocity(x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, t + 0.5 * dt);
            [k4x, k4y] = app.velocity(x + dt * k3x, y + dt * k3y, t + dt);

            xn = x + dt * (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
            yn = y + dt * (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        end

        %% -------------------- Pathline (precompute + marker) --------------------
        function precomputePathline(app)
            % Precompute a pathline from t=0 using a fixed number of samples
            %
            % Notes:
            %   * This draws the pathline "to infinity" in practice by choosing
            %     numPointsPathline and DTPATHLINE sufficiently large/small.
            %   * The marker position is later interpolated from (pathT,pathX,pathY).

            [x0, y0] = getX0Y0(app);

            N = max(2, round(app.numPointsPathline));
            dt = max(1e-4, app.DTPATHLINE);

            app.pathX = nan(1, N);
            app.pathY = nan(1, N);
            app.pathT = nan(1, N);

            app.pathX(1) = x0;
            app.pathY(1) = y0;
            app.pathT(1) = 0;

            % Bounding box used only to decide if a zoom-out should trigger recompute.
            xl = app.Axes.XLim; yl = app.Axes.YLim;
            cx = mean(xl); cy = mean(yl);
            rx = 1.5 * 0.5 * diff(xl);
            ry = 1.5 * 0.5 * diff(yl);
            app.pathBox = [cx - rx, cx + rx, cy - ry, cy + ry];

            for k = 2:N
                [xn, yn] = rk4(app, app.pathX(k-1), app.pathY(k-1), app.pathT(k-1), dt);
                app.pathX(k) = xn;
                app.pathY(k) = yn;
                app.pathT(k) = app.pathT(k-1) + dt;
            end
        end

        function [xp, yp] = pathPosAtTime(app, tNow)
            % Interpolate precomputed pathline to get marker position at time tNow

            if isempty(app.pathT) || numel(app.pathT) < 2
                [x0, y0] = getX0Y0(app);
                xp = x0; yp = y0;
                return
            end

            dt = app.DTPATHLINE;
            N  = numel(app.pathX);

            s = tNow / dt;
            i = floor(s) + 1;
            if i < 1, xp = app.pathX(1); yp = app.pathY(1); return; end
            if i >= N, xp = app.pathX(end); yp = app.pathY(end); return; end

            a = s - floor(s);  % in [0,1)
            xp = (1-a)*app.pathX(i) + a*app.pathX(i+1);
            yp = (1-a)*app.pathY(i) + a*app.pathY(i+1);
        end

        %% -------------------- Streamlines numeric helpers --------------------
        function [dx, dy] = unitDir(app, x, y, tFix)
            % Unit tangent direction for streamline integration:
            %     dX/ds = u/||u||, dY/ds = v/||u||

            [u, v] = app.velocity(x, y, tFix);
            s = hypot(u, v);
            if s < 1e-12
                dx = 0;
                dy = 0;
            else
                dx = u / s;
                dy = v / s;
            end
        end

        function [XS, YS] = streamlinesUniform(app, tFix, n)
            x = linspace(app.xMin, app.xMax, 200);
            [u0, v0] = app.velocity(0, 0, tFix);

            slope = 0;
            if abs(u0) > 1e-12
                slope = v0 / u0;
            end

            xRef = x(1);
            C = linspace(app.yMin, app.yMax, n);

            XS = cell(n,1);
            YS = cell(n,1);
            for i = 1:n
                XS{i} = x;
                YS{i} = C(i) + slope*(x - xRef);
            end
        end

        function [XS, YS] = streamlinesSaddle(app, n)
            x = linspace(app.xMin, app.xMax, 200);

            XS = cell(n,1);
            YS = cell(n,1);

            if n >= 1
                XS{1} = x; YS{1} = 0*x;
            end
            if n >= 2
                yAxis = linspace(app.yMin, app.yMax, 800);
                XS{2} = 0*yAxis; YS{2} = yAxis;
            end

            m = max(n-2,0);
            if m > 0
                Cmax = 0.35*(max(abs([app.xMin, app.xMax]))*max(abs([app.yMin, app.yMax])) + 1);
                Cvals = linspace(-Cmax, Cmax, m);
                Cvals(abs(Cvals) < 1e-9) = 0.15*Cmax;

                xx = x;
                epsx = 1e-6*max(1, max(abs(x)));
                xx(abs(xx) < epsx) = epsx;

                for k = 1:m
                    i = k + 2;
                    XS{i} = xx;
                    YS{i} = Cvals(k) ./ xx;
                end
            end

            % Fill any empties if n small/odd cases
            for i = 1:n
                if isempty(XS{i}), XS{i} = nan; YS{i} = nan; end
            end
        end


        function [xs, ys] = integrateStreamline(app, x0, y0, tFix)
            % Numerically integrate a streamline from a seed point
            %
            % Streamlines satisfy:
            %     dy/dx = v/u  (or in arclength form: dX/ds = u/||u||).
            %
            % Args:
            %     x0, y0 (double): seed point
            %     tFix (double): time at which to compute instantaneous streamlines
            %
            % Returns:
            %     xs, ys (double): polyline describing the streamline

            % Definitions
            dom = max(app.xMax - app.xMin, app.yMax - app.yMin);
            ds   = max(0.03 * dom, 0.02);
            maxN = 200;
            vTol = 1e-10;
            margin = 0.5;

            % Cache bounds and constants
            xLo = app.xMin - margin;  xHi = app.xMax + margin;
            yLo = app.yMin - margin;  yHi = app.yMax + margin;

            % Preallocate
            Xb = zeros(maxN, 1); Yb = zeros(maxN, 1);
            Xf = zeros(maxN, 1); Yf = zeros(maxN, 1);

            Xb(1) = x0; Yb(1) = y0; nb = 1;
            Xf(1) = x0; Yf(1) = y0; nf = 1;

            % Integrate backward and forward
            [Xb, Yb, nb] = integrateDir(x0, y0, -ds, Xb, Yb, nb);
            [Xf, Yf, nf] = integrateDir(x0, y0, +ds, Xf, Yf, nf);

            % Merge (avoid duplicating seed point)
            xs = [flipud(Xb(1:nb)); Xf(2:nf)];
            ys = [flipud(Yb(1:nb)); Yf(2:nf)];

            % NESTED FUNCTIONS
            function [X, Y, n] = integrateDir(x, y, h, X, Y, n)
                for k = 2:maxN
                    [x, y, ok] = rk4Step(x, y, h);
                    if ~ok
                        break;
                    end
                    n = n + 1;
                    X(n) = x;
                    Y(n) = y;
                end
            end

            function [xn, yn, ok] = rk4Step(x, y, h)
                % One RK4 step: dX/ds = unitDir(X), with step h (±ds).

                [d1x, d1y] = app.unitDir(x, y, tFix);

                if ~(isfinite(d1x) && isfinite(d1y))
                    xn = x; yn = y; ok = false; return;
                end

                [d2x, d2y] = app.unitDir(x + 0.5*h*d1x, y + 0.5*h*d1y, tFix);
                [d3x, d3y] = app.unitDir(x + 0.5*h*d2x, y + 0.5*h*d2y, tFix);
                [d4x, d4y] = app.unitDir(x + h*d3x, y + h*d3y, tFix);

                xn = x + h * (d1x + 2*(d2x + d3x) + d4x) / 6;
                yn = y + h * (d1y + 2*(d2y + d3y) + d4y) / 6;

                % Domain check
                ok = (xn >= xLo) && (xn <= xHi) && (yn >= yLo) && (yn <= yHi);
                if ~ok
                    return;
                end

                % Stop if velocity magnitude is effectively zero
                [ux, uy] = app.velocity(xn, yn, tFix);
                ok = hypot(ux, uy) >= vTol;
            end
        end

        %% -------------------- Axis limits --------------------
        function computeAxisLimits(app)
            % Compute axis limits for the selected problem and (x0,y0)
            %
            % Notes:
            %   * For uniform flows, x-range depends on V0 and TmaxConst
            %   * For other cases, a default window around (x0,y0) is used

            F = app.Flows(app.flowIndex);
            lim = F.AxisLimits(app); % [xMin xMax yMin yMax]

            app.xMin = lim(1);
            app.xMax = lim(2);
            app.yMin = lim(3);
            app.yMax = lim(4);
        end

        %% -------------------- View-dependent quiver --------------------
        function installViewListeners(app)
            if ~isempty(app.xLimListener) && isvalid(app.xLimListener), return; end
            app.xLimListener = addlistener(app.Axes, 'XLim', 'PostSet', @(~,~) markViewDirty(app));
            app.yLimListener = addlistener(app.Axes, 'YLim', 'PostSet', @(~,~) markViewDirty(app));
        end

        function onViewChanged(app)
            if app.isResetting || ~isgraphics(app.Axes)
                return
            end

            xl = app.Axes.XLim; yl = app.Axes.YLim;
            app.xMin = xl(1); app.xMax = xl(2);
            app.yMin = yl(1); app.yMax = yl(2);

            rebuildQuiverGrid(app);
            refreshQuiver(app);

            % Streamlines depend on xMin/xMax/yMin/yMax
            updateStreamlines(app);

            % If user zooms out a lot, recompute the precomputed pathline so it remains representative.
            if ~isPathBoxCoveringView(app)
                precomputePathline(app);
                if ~isempty(app.pathLine) && isgraphics(app.pathLine)
                    set(app.pathLine, 'XData', app.pathX, 'YData', app.pathY);
                end
            end

        end

        function tf = isPathBoxCoveringView(app)
            if numel(app.pathBox) ~= 4
                tf = false;
                return
            end
            xl = app.Axes.XLim; yl = app.Axes.YLim;
            tf = (xl(1) >= app.pathBox(1)) && (xl(2) <= app.pathBox(2)) && ...
                 (yl(1) >= app.pathBox(3)) && (yl(2) <= app.pathBox(4));
        end

        function rebuildQuiverGrid(app)
            % Always build grid based on current view limits
            xl = app.Axes.XLim; yl = app.Axes.YLim;
            [app.xg, app.yg] = meshgrid( ...
                linspace(xl(1), xl(2), app.qNx), ...
                linspace(yl(1), yl(2), app.qNy));
        end

        function refreshQuiver(app)
            % Create/update quiver using the current grid
            [u, v] = app.velocity(app.xg, app.yg, app.t);

            if isempty(app.qh) || ~isgraphics(app.qh)
                app.qh = quiver(app.Axes, app.xg, app.yg, u, v, 'AutoScale', 'on');
                app.qh.DisplayName = 'Velocity field';
                app.qh.HandleVisibility = 'on';
                app.qh.HitTest = 'off';
                app.qh.PickableParts = 'none';

                return
            end

            set(app.qh, 'XData', app.xg, 'YData', app.yg, 'UData', u, 'VData', v);
        end

        %% -------------------- Plot setup --------------------
        function resetAll(app)
            % Full reset: stop timer, clear axes, recreate graphics, reset histories
            %
            % Notes:
            %   * This is the ONLY operation that resets time and traces

            if app.isResetting
                return
            end
            app.isResetting = true;
            c = onCleanup(@() setResetFalse(app)); %#ok<NASGU>

            stopTimer(app);
            while app.isTicking
                pause(0.01);
            end

            % Delete all children (including HandleVisibility='off')
            try
                delete(allchild(app.Axes));
            catch
                try
                    delete(app.Axes.Children);
                catch
                end
            end

            drawnow;

            % Reset handle caches
            app.qh = [];
            app.srcDot = [];
            app.streamH = gobjects(0);
            app.pathLine = [];
            app.pDot = [];
            app.streakLine = [];
            app.streakDots = [];

            % Reset time + histories using (x0, y0)
            [x0, y0] = getX0Y0(app);

            app.isRunning = false;
            app.t = 0;
            app.tickCount = uint64(0);

            app.pathX = x0;
            app.pathY = y0;
            app.pathT = 0;

            app.streakX = x0;
            app.streakY = y0;
            app.lastReleaseT = 0;

            configureAxesAppearance(app);
            hold(app.Axes, 'on');

            computeAxisLimits(app);
            axis(app.Axes, [app.xMin, app.xMax, app.yMin, app.yMax]);

            % Keep quiver in sync with zoom/pan
            installViewListeners(app);

            % Quiver grid
            rebuildQuiverGrid(app);
            refreshQuiver(app);

            % Source / release point
            app.srcDot = plot(app.Axes, x0, y0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
            app.srcDot.HandleVisibility = 'off';

            % Precompute pathline polyline once (marker will move with time)
            precomputePathline(app);

            ensureStreamObjects(app);
            ensurePathObjects(app);
            ensureStreakObjects(app);

            updateStreamlines(app);
            updateColors(app);
            toggleVisibility(app);

            % Update plot
            updatePlot(app, 0);

            % Update title
            updateTitle(app);

            function setResetFalse(appObj)
                if isvalid(appObj)
                    appObj.isResetting = false;
                end
            end
        end

        function ensureStreamObjects(app)
            % Create streamline line objects (one per streamline)
            if ~isempty(app.streamH) && all(isgraphics(app.streamH))
                return
            end

            app.streamH = gobjects(app.numStreamlines, 1);
            for i = 1:app.numStreamlines
                app.streamH(i) = plot(app.Axes, nan, nan, '-', 'LineWidth', 1.2);
                app.streamH(i).HandleVisibility = 'off';
                app.streamH(i).HitTest = 'off';
                app.streamH(i).PickableParts = 'none';
            end
        end

        function ensurePathObjects(app)
            % Create pathline polyline and particle marker
            if ~isempty(app.pathLine) && isgraphics(app.pathLine) && ...
               ~isempty(app.pDot) && isgraphics(app.pDot)
                return
            end

            app.pathLine = plot(app.Axes, app.pathX, app.pathY, '-', 'LineWidth', 2.2, 'DisplayName', 'Pathline');
            app.pDot = plot(app.Axes, app.pathX(1), app.pathY(1), 'o', ...
                'MarkerFaceColor', app.PathColorPicker.Value, 'MarkerSize', 7);
            app.pDot.HandleVisibility = 'off';
            app.pathLine.HitTest = 'off';  app.pathLine.PickableParts = 'none';
            app.pDot.HitTest = 'off';      app.pDot.PickableParts = 'none';
        end

        function ensureStreakObjects(app)
            % Create streakline polyline and scatter points
            if ~isempty(app.streakLine) && isgraphics(app.streakLine) && ...
               ~isempty(app.streakDots) && isgraphics(app.streakDots)
                return
            end

            app.streakLine = plot(app.Axes, app.streakX, app.streakY, '-', 'LineWidth', 2.2, 'DisplayName', 'Streakline');
            app.streakDots = plot(app.Axes, app.streakX, app.streakY, 'o', ...
                'LineStyle','none', 'MarkerSize', 4, 'MarkerFaceColor', app.StreakColorPicker.Value, ...
                'MarkerEdgeColor', app.StreakColorPicker.Value);
            app.streakDots.HandleVisibility = 'off';
            app.streakLine.HitTest = 'off'; app.streakLine.PickableParts = 'none';
            app.streakDots.HitTest = 'off'; app.streakDots.PickableParts = 'none';
        end

        %% -------------------- Streamlines (analytic where possible) --------------------
        function updateStreamlines(app)
            % Update streamline polylines at current time
            %
            % Notes:
            %   * Uniform cases: explicit y(x) formula
            %   * Several steady canonical fields: analytic families
            %   * Fallback: numeric streamline integration

            if isempty(app.streamH) || ~all(isgraphics(app.streamH))
                return
            end

            % Definitions
            n = numel(app.streamH);
            F = app.Flows(app.flowIndex);

            if ~isempty(F.Streamlines)
                [XS, YS] = F.Streamlines(app, app.t, n);
                for i = 1:n
                    setLine(i, XS{i}, YS{i});
                end
            else
                % Fallback numeric seeding (your existing fallback)
                xSeed = app.xMin + 0.08 * (app.xMax - app.xMin);
                ySeeds = linspace(app.yMin, app.yMax, n);

                for i = 1:n
                    [xs, ys] = integrateStreamline(app, xSeed, ySeeds(i), app.t);
                    setLine(i, xs, ys);
                end
            end

            % SUB-PASS FUNCTIONS
            function setLine(i, x, y)
                % Helper to set streamline line data with axis clipping (NaN outside limits)
                mask = (x < app.xMin) | (x > app.xMax) | (y < app.yMin) | (y > app.yMax);
                x(mask) = NaN;
                y(mask) = NaN;
                set(app.streamH(i), 'XData', x, 'YData', y);
            end

        end

        %% -------------------- Colors --------------------
        function updateColors(app)
            % Apply UI color pickers to plotted objects

            % Get colors from pickers
            colorVelocity    = app.VelColorPicker.Value;
            colorStreamline  = app.StreamColorPicker.Value;
            colorPathline    = app.PathColorPicker.Value;
            colorStreakline  = app.StreakColorPicker.Value;

            if ~isempty(app.qh) && isgraphics(app.qh)
                app.qh.Color = colorVelocity;
            end

            if ~isempty(app.streamH) && all(isgraphics(app.streamH))
                for i = 1:numel(app.streamH)
                    app.streamH(i).Color = colorStreamline;
                end
            end

            if ~isempty(app.pathLine) && isgraphics(app.pathLine)
                app.pathLine.Color = colorPathline;
            end

            if ~isempty(app.pDot) && isgraphics(app.pDot)
                app.pDot.Color = colorPathline;
                app.pDot.MarkerFaceColor = colorPathline;
            end

            if ~isempty(app.streakLine) && isgraphics(app.streakLine)
                app.streakLine.Color = colorStreakline;
            end

            if ~isempty(app.streakDots) && isgraphics(app.streakDots)
                app.streakDots.MarkerFaceColor = colorStreakline;
                app.streakDots.MarkerEdgeColor = colorStreakline;
            end

        end

        %% -------------------- Visibility --------------------
        function toggleVisibility(app)
            % Toggle object visibility from checkboxes

            setOnOff(app.qh, app.VelCheckBox.Value);
            setOnOff(app.streamH, app.StreamCheckBox.Value);
            setOnOff([app.pathLine, app.pDot], app.PathCheckBox.Value);
            setOnOff([app.streakLine, app.streakDots], app.StreakCheckBox.Value);

            function setOnOff(h, tf)
                if isempty(h), return; end
                for k = 1:numel(h)
                    if ~isgraphics(h(k)), continue; end
                    h(k).Visible = tern(tf, 'on', 'off');
                end
            end

            function out = tern(cond, a, b)
                if cond, out = a; else, out = b; end
            end
        end

        function updateTitle(app)
            % Update title with current time and velocity info
            timeLine = sprintf('$t = %.2f\\,\\mathrm{s}$', app.t);
            app.Axes.Title.String = timeLine;
        end

        %% -------------------- Update per frame --------------------
        function updatePlot(app, dt)
            % Advance visualization by one frame.
            %
            % Args:
            %     dt (double): time step used for advection (0 for a "refresh")

            if app.isResetting || ~isgraphics(app.Axes)
                return
            end

            % Apply view changes at most every viewApplyMinPeriod
            if app.viewDirty && (app.t - app.lastViewApplyT) >= app.viewApplyMinPeriod
                app.viewDirty = false;
                app.lastViewApplyT = app.t;
                onViewChanged(app);
            end

            % Update quiver
            if ~isempty(app.qh) && isgraphics(app.qh) && app.VelCheckBox.Value && app.isUnsteady()
                [u, v] = app.velocity(app.xg, app.yg, app.t);
                set(app.qh, 'UData', u, 'VData', v);
            end

            % Update streamlines only if unsteady
            if app.isUnsteady() && app.StreamCheckBox.Value
                updateStreamlines(app);
            end

            % NOTE: The pathline curve is precomputed from t=0 and stays fixed;
            %       only the marker position depends on the current app.t.
            [xp, yp] = pathPosAtTime(app, app.t);

            if dt > 0
                app.tickCount = app.tickCount + 1;

                t0 = app.t - dt;
                [x0, y0] = getX0Y0(app);

                % Streakline: release new particles at fixed time intervals, then advect all
                nNew = floor((app.t - app.lastReleaseT) / app.DTSTREAK);
                if nNew > 0
                    [x0,y0] = getX0Y0(app);
                    app.streakX = [app.streakX, repmat(x0, 1, nNew)];
                    app.streakY = [app.streakY, repmat(y0, 1, nNew)];
                    app.lastReleaseT = app.lastReleaseT + nNew * app.DTSTREAK;

                    ns = numel(app.streakX);
                    if ns > app.maxStreakPts
                        app.streakX = app.streakX(ns - app.maxStreakPts + 1 : ns);
                        app.streakY = app.streakY(ns - app.maxStreakPts + 1 : ns);
                    end
                end

                % Streak advect
                [app.streakX, app.streakY] = rk4(app, app.streakX, app.streakY, t0, dt);
            end

            % Push data to graphics
            if ~isempty(app.pDot) && isgraphics(app.pDot) && app.PathCheckBox.Value
                set(app.pDot, 'XData', xp, 'YData', yp);
            end

            if ~isempty(app.streakLine) && isgraphics(app.streakLine) && app.StreakCheckBox.Value
                set(app.streakLine, 'XData', app.streakX, 'YData', app.streakY);
            end

            if ~isempty(app.streakDots) && isgraphics(app.streakDots) && app.StreakCheckBox.Value
                set(app.streakDots, 'XData', app.streakX, 'YData', app.streakY);
            end

            % Update title
            updateTitle(app);

            % Update graphics
            drawnow limitrate nocallbacks
        end

        %% -------------------- Timer --------------------
        function onRun(app),  startTimer(app); end
        function onStop(app), stopTimer(app);  end

        function startTimer(app)
            % Start (or restart) the fixed-rate timer
            %
            % Notes:
            %   * Period is quantized to milliseconds to avoid MATLAB timer warnings

            stopTimer(app);
            app.isRunning = true;

            dt = 1 / app.FPS;
            dt = round(dt * 1000) / 1000;
            dt = max(dt, 0.001);

            app.Tmr = timer( ...
                'ExecutionMode', 'fixedRate', ...
                'Period', dt, ...
                'BusyMode', 'drop', ...
                'TimerFcn', @(~, ~) onTick(app, dt));
            start(app.Tmr);
        end

        function stopTimer(app)
            % Stop and delete timer safely
            app.isRunning = false;

            if ~isempty(app.Tmr) && isvalid(app.Tmr)
                try, stop(app.Tmr); catch, end
                if app.isTicking
                    return
                end

                try, delete(app.Tmr); catch, end
            end

            app.Tmr = timer.empty;
        end

        function onTick(app, dt)
            % Timer tick: advance time and update plots
            %
            % Args:
            %     dt (double): timer period (time increment)
            if app.isTicking
                return
            end

            app.isTicking = true;
            c = onCleanup(@() setTickFalse(app));

            if app.isResetting || ~app.isRunning || ~isgraphics(app.Axes)
                return
            end

            app.t = app.t + dt;
            updatePlot(app, dt);

            function setTickFalse(appObj)
                if ~isvalid(appObj)
                    return
                end

                appObj.isTicking = false;
            end

        end

        function markViewDirty(app)
            app.viewDirty = true;
        end

    end
end
