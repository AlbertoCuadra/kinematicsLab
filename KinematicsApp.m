classdef KinematicsApp < matlab.apps.AppBase
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
    % Last update: 18/02/2026

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
        % Fixed simulation settings
        fpsConst double = 30          % Fixed FPS
        TmaxConst double = 12         % Used only for axis sizing in uniform cases

        % Streamline
        numStream int32 = 9           % Number of streamlines

        % Performance limits (keep only most recent samples)
        maxPathPts double = 6000      % Pathline polyline points
        maxStreakPts double = 1200    % Streak particles kept
        pathStride double = 1         % Store every k ticks (1 = store all)
        tickCount uint64 = uint64(0)  % Number of ticks since the start of the simulation

        % Timer
        Tmr = timer.empty             % Timer object for the simulation loop
        t double = 0                  % Simulation time (seconds)
        isRunning logical = false     % Flag for simulation state

        % Flags for clear
        isResetting logical = false   % Flag to avoid running the timer callback during reset
        isTicking logical = false     % Flag to avoid running the reset callback during timer tick

        % Graphics handles
        qh
        srcDot
        streamH = gobjects(0)         % Streamlines lines
        pathLine                      % Pathline
        pDot                          % Particle marker
        streakLine                    % Streakline
        streakDots                    % Streakdots

        % Cached geometry
        Cvals double                  % Streamlines
        dtRel double = 0.12           % Streak release interval (seconds)
        xg double                     % Grid of x values for streamlines
        yg double                     % Grid of y values for streamlines

        % Axis limits
        xMin double                  % Axis limits (computed from the velocity field and TmaxConst)
        xMax double                  % Axis limits (computed from the velocity field and TmaxConst)
        yMin double                  % Axis limits (computed from the velocity field and TmaxConst)
        yMax double                  % Axis limits (computed from the velocity field and TmaxConst)

        % Histories (RK)
        pathX double = 0             % Pathline history of x positions
        pathY double = 0             % Pathline history of y positions
        streakX double = 0           % Streakline history of x positions
        streakY double = 0           % Streakline history of y positions
        lastReleaseT double = 0      % Last time a streak particle was released

        % Interface
        fontsizeLabels = 18          % Font size for labels and buttons

        % Flags
        showLegend = false           % Flag to show legend
    end

    methods (Access = public)
        function app = KinematicsApp
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

    end

    methods (Access = private)

        %% -------------------- Startup --------------------
        function startup(app)
            % Initialize defaults and reset the app.
            %
            % Args:
            %     app (KinematicsApp): App instance

            app.ProblemDropDown.Items = { ...
                'Unsteady uniform: u=V0, v=V0 sin(ωt)', ...
                'Steady uniform: u=V0, v=0', ...
                'Steady saddle: u=x, v=-y', ...
                'Steady source: u=x, v=y', ...
                'Steady rotation: u=y, v=-x', ...
                'Steady shear: u=V0, v=V0 x', ...
                'Steady cellular: u=cos(ωx), v=sin(ωy)', ...
                'Steady polynomial: u=(x^2+x-2)(x^2+x-12), v=y-3x+7' ...
                };
            app.ProblemDropDown.Value = app.ProblemDropDown.Items{1};

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
        end

        %% -------------------- UI --------------------
        function createComponents(app)
            % Create UI layout and wire callbacks.
            %
            % Args:
            %     app (KinematicsApp): App instance
            
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

            % Left-panel grid (scrollable)
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

            % Velocity field row
            app.VelCheckBox = uicheckbox(lp, ...
                'Text', 'Velocity field', ...
                'ValueChangedFcn', @(~, ~) toggleVisibility(app), ...
                'FontSize', app.fontsizeLabels - 4);
            app.VelCheckBox.Layout.Row = 2;
            app.VelCheckBox.Layout.Column = [1, 2];

            app.VelColorPicker = uicolorpicker(lp, ...
                'Value', [204, 204, 204] / 255, ...
                'ValueChangedFcn', @(~, ~) colorChanged(app));
            app.VelColorPicker.Layout.Row = 2;
            app.VelColorPicker.Layout.Column = 3;

            % Streamline row
            app.StreamCheckBox = uicheckbox(lp, ...
                'Text', 'Streamline', ...
                'ValueChangedFcn', @(~, ~) toggleVisibility(app), ...
                'FontSize', app.fontsizeLabels - 4);
            app.StreamCheckBox.Layout.Row = 3;
            app.StreamCheckBox.Layout.Column = [1, 2];

            app.StreamColorPicker = uicolorpicker(lp, ...
                'Value', [0.91, 0.51, 0.49], ...
                'ValueChangedFcn', @(~, ~) colorChanged(app));
            app.StreamColorPicker.Layout.Row = 3;
            app.StreamColorPicker.Layout.Column = 3;

            % Pathline row
            app.PathCheckBox = uicheckbox(lp, ...
                'Text', 'Pathline', ...
                'ValueChangedFcn', @(~, ~) toggleVisibility(app), ...
                'FontSize', app.fontsizeLabels - 4);
            app.PathCheckBox.Layout.Row = 4;
            app.PathCheckBox.Layout.Column = [1, 2];

            app.PathColorPicker = uicolorpicker(lp, ...
                'Value', [89, 156, 155] / 255, ...
                'ValueChangedFcn', @(~, ~) colorChanged(app));
            app.PathColorPicker.Layout.Row = 4;
            app.PathColorPicker.Layout.Column = 3;

            % Streakline row
            app.StreakCheckBox = uicheckbox(lp, ...
                'Text', 'Streakline', ...
                'ValueChangedFcn', @(~, ~) toggleVisibility(app), ...
                'FontSize', app.fontsizeLabels - 4);
            app.StreakCheckBox.Layout.Row = 5;
            app.StreakCheckBox.Layout.Column = [1, 2];

            app.StreakColorPicker = uicolorpicker(lp, ...
                'Value', [152, 209, 208] / 255, ...
                'ValueChangedFcn', @(~, ~) colorChanged(app));
            app.StreakColorPicker.Layout.Row = 5;
            app.StreakColorPicker.Layout.Column = 3;

            % V0 row
            app.V0Label = uilabel(lp, 'Text', '$V_0$ [m/s]', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            app.V0Label.Interpreter = 'latex';
            app.V0Label.Layout.Row = 6;
            app.V0Label.Layout.Column = 1;

            app.V0Field = uieditfield(lp, 'numeric', 'Value', 1.0, 'Limits', [-Inf, Inf], ...
                'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.V0Field.Layout.Row = 6;
            app.V0Field.Layout.Column = [2, 3];

            % x0 row
            app.X0Label = uilabel(lp, 'Text', '$x_0$ [m]', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            app.X0Label.Interpreter = 'latex';
            app.X0Label.Layout.Row = 7;
            app.X0Label.Layout.Column = 1;

            app.X0Field = uieditfield(lp, 'numeric', 'Value', 0.0, 'Limits', [-Inf, Inf], ...
                'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.X0Field.Layout.Row = 7;
            app.X0Field.Layout.Column = [2, 3];

            % y0 row
            app.Y0Label = uilabel(lp, 'Text', '$y_0$ [m]', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            app.Y0Label.Interpreter = 'latex';
            app.Y0Label.Layout.Row = 8;
            app.Y0Label.Layout.Column = 1;

            app.Y0Field = uieditfield(lp, 'numeric', 'Value', 0.0, 'Limits', [-Inf, Inf], ...
                'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.Y0Field.Layout.Row = 8;
            app.Y0Field.Layout.Column = [2, 3];

            % Omega row
            app.OmegaLabel = uilabel(lp, 'Text', '$\omega$ [1/s]', 'HorizontalAlignment', 'left', 'FontSize', app.fontsizeLabels - 4);
            app.OmegaLabel.Interpreter = 'latex';
            app.OmegaLabel.Layout.Row = 9;
            app.OmegaLabel.Layout.Column = 1;

            app.OmegaField = uieditfield(lp, 'numeric', 'Value', 2, 'Limits', [-Inf, Inf], ...
                'ValueChangedFcn', @(~, ~) paramChanged(app), 'FontSize', app.fontsizeLabels - 4);
            app.OmegaField.Layout.Row = 9;
            app.OmegaField.Layout.Column = [2, 3];

            % Buttons row (nested grid)
            btnGrid = uigridlayout(lp, [1, 2]);
            btnGrid.Layout.Row = 11;
            btnGrid.Layout.Column = [1, 3];
            btnGrid.ColumnWidth = {'1x', '1x'};
            btnGrid.RowHeight = {30};
            btnGrid.Padding = [0, 0, 0, 0];
            btnGrid.ColumnSpacing = 10;

            app.RunButton = uibutton(btnGrid, 'Text', 'run', 'ButtonPushedFcn', @(~, ~) onRun(app), 'FontSize', app.fontsizeLabels - 4);
            app.StopButton = uibutton(btnGrid, 'Text', 'stop', 'ButtonPushedFcn', @(~, ~) onStop(app), 'FontSize', app.fontsizeLabels - 4);

            % Clear row
            app.ClearButton = uibutton(lp, 'Text', 'clear', 'ButtonPushedFcn', @(~, ~) resetAll(app), 'FontSize', app.fontsizeLabels - 4);
            app.ClearButton.Layout.Row = 12;
            app.ClearButton.Layout.Column = [1, 3];
        end

        function configureAxesAppearance(app)
            % Configure axes aesthetics and LaTeX labels.
            %
            % Args:
            %     app (KinematicsApp): App instance

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
            updateTitle(app);

            app.Axes.XGrid = 'off';
            app.Axes.YGrid = 'off';
            app.Axes.XMinorGrid = 'off';
            app.Axes.YMinorGrid = 'off';
        end

        function onClose(app)
            % UIFigure CloseRequestFcn callback. Stops timer and deletes the app.
            stopTimer(app);
            delete(app);
        end

        function paramChanged(app)
            % Callback for parameter changes (problem / V0 / x0 / y0 / omega).
            % Stops simulation, resets all graphics & histories, keeps running if it was running.
            wasRunning = app.isRunning;
            updateOmegaEnable(app);
            updateV0Enable(app);
            resetAll(app);
            
            if wasRunning
                startTimer(app);
            end

        end

        function colorChanged(app)
            % Callback for color pickers. Updates the colors of all graphics objects.
            updateColors(app);
            updateLegend(app);
        end

        %% ------------------ Enable/Disable parameters ------------------
        function updateOmegaEnable(app)
            % Enable/disable omega input depending on the selected problem

            % Definitions
            needsOmega = problemUsesOmega(app);

            % Update UI
            if needsOmega
                app.OmegaField.Enable = 'on';
                return
            end
            
            app.OmegaField.Enable = 'off';
        end

        function tf = problemUsesOmega(app)
            % Return true if the current problem definition uses omega
            %
            % Returns:
            %     tf (logical): True if omega is used by the selected velocity field
            val = app.ProblemDropDown.Value;
            tf = contains(val, 'ω');
        end

        function updateV0Enable(app)
            % Enable/disable V0 input depending on the selected problem

            % Definitions
            needsV0 = problemUsesV0(app);

            % Update UI
            if needsV0
                app.V0Field.Enable = 'on';
                return
            end

            app.V0Field.Enable = 'off';
        end

        function tf = problemUsesV0(app)
            % Return true if the current problem definition uses V0
            %
            % Returns:
            %     tf (logical): True if V0 is used by the selected velocity field
            val = app.ProblemDropDown.Value;
            tf = contains(val, 'V0');
        end

        function tf = isUnsteady(app)
            % Return true if the selected velocity field is time-dependent.
            tf = strcmp(app.ProblemDropDown.Value, 'Unsteady uniform: u=V0, v=V0 sin(ωt)');
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

            V0 = app.V0Field.Value;
            om = app.OmegaField.Value;
            problemCase = app.ProblemDropDown.Value;

            switch problemCase
                case 'Unsteady uniform: u=V0, v=V0 sin(ωt)'
                    u = V0 + 0 * x;
                    v = V0 * sin(om * t) + 0 * y;

                case 'Steady uniform: u=V0, v=0'
                    u = V0 + 0 * x;
                    v = 0 + 0 * y;

                case 'Steady saddle: u=x, v=-y'
                    u = x;
                    v = -y;

                case 'Steady source: u=x, v=y'
                    u = x;
                    v = y;

                case 'Steady rotation: u=y, v=-x'
                    u = y;
                    v = -x;

                case 'Steady shear: u=V0, v=V0 x'
                    u = V0 + 0 * x;
                    v = V0 * x;

                case 'Steady cellular: u=cos(ωx), v=sin(ωy)'
                    u = cos(om * x);
                    v = sin(om * y);

                case 'Steady polynomial: u=(x^2+x-2)(x^2+x-12), v=y-3x+7'
                    u = (x .^ 2 + x - 2) .* (x .^ 2 + x - 12);
                    v = y - 3 * x + 7;

                otherwise
                    u = 0 + 0 * x;
                    v = 0 + 0 * y;
            end
        end

        %% -------------------- RK4 helpers --------------------
        function [xn, yn] = rk4Particle(app, x, y, t, dt)
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

        function [xn, yn] = rk4ParticleVec(app, x, y, t, dt)
            % Vectorized RK4 for arrays x, y (streak particles)
            %
            % Args:
            %     x, y (double): arrays of particle positions
            %     t (double): current time
            %     dt (double): time step
            %
            % Returns:
            %     xn, yn (double): arrays of updated particle positions

            [k1x, k1y] = app.velocity(x, y, t);
            [k2x, k2y] = app.velocity(x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, t + 0.5 * dt);
            [k3x, k3y] = app.velocity(x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, t + 0.5 * dt);
            [k4x, k4y] = app.velocity(x + dt * k3x, y + dt * k3y, t + dt);

            xn = x + dt * (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
            yn = y + dt * (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        end

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

            ds = 0.03 * max(app.xMax - app.xMin, app.yMax - app.yMin);
            ds = max(ds, 0.02);
            maxN = 700;

            % Backward integration
            xb = x0;
            yb = y0;
            Xb = nan(maxN, 1);
            Yb = nan(maxN, 1);
            Xb(1) = xb;
            Yb(1) = yb;
            nb = 1;

            for k = 2:maxN
                [xb, yb, ok] = stepStream(app, xb, yb, tFix, -ds);
                if ~ok
                    break;
                end
                nb = nb + 1;
                Xb(nb) = xb;
                Yb(nb) = yb;
            end

            % Forward integration
            xf = x0;
            yf = y0;
            Xf = nan(maxN, 1);
            Yf = nan(maxN, 1);
            Xf(1) = xf;
            Yf(1) = yf;
            nf = 1;

            for k = 2:maxN
                [xf, yf, ok] = stepStream(app, xf, yf, tFix, +ds);
                if ~ok
                    break;
                end
                nf = nf + 1;
                Xf(nf) = xf;
                Yf(nf) = yf;
            end

            % Merge (avoid duplicating seed point)
            Xb = flipud(Xb(1:nb));
            Yb = flipud(Yb(1:nb));
            Xf = Xf(2:nf);
            Yf = Yf(2:nf);

            xs = [Xb; Xf];
            ys = [Yb; Yf];

            function [xn, yn, ok] = stepStream(appObj, x, y, tF, h)
                % One RK4 step in arclength parameter s for the unit direction field
                [d1x, d1y] = appObj.unitDir(x, y, tF);
                [d2x, d2y] = appObj.unitDir(x + 0.5 * h * d1x, y + 0.5 * h * d1y, tF);
                [d3x, d3y] = appObj.unitDir(x + 0.5 * h * d2x, y + 0.5 * h * d2y, tF);
                [d4x, d4y] = appObj.unitDir(x + h * d3x, y + h * d3y, tF);

                xn = x + h * (d1x + 2 * d2x + 2 * d3x + d4x) / 6;
                yn = y + h * (d1y + 2 * d2y + 2 * d3y + d4y) / 6;

                margin = 0.5;
                ok = xn >= (appObj.xMin - margin) && xn <= (appObj.xMax + margin) && ...
                     yn >= (appObj.yMin - margin) && yn <= (appObj.yMax + margin);

                if ~ok
                    return;
                end

                [ux, uy] = appObj.velocity(xn, yn, tF);
                if hypot(ux, uy) < 1e-10
                    ok = false;
                end
            end
        end

        %% -------------------- Axis limits --------------------
        function computeAxisLimits(app)
            % Compute axis limits for the selected problem and (x0,y0)
            %
            % Notes:
            %   * For uniform flows, x-range depends on V0 and TmaxConst
            %   * For other cases, a default window around (x0,y0) is used

            V0 = app.V0Field.Value;
            om = app.OmegaField.Value;
            val = app.ProblemDropDown.Value;
            [x0, y0] = getX0Y0(app);

            switch val
                case 'Unsteady uniform: u=V0, v=V0 sin(ωt)'
                    xMin0 = -1;
                    xMax0 = x0 + V0 * max(app.TmaxConst, 0) + 1;

                    if om ~= 0
                        A = 2 * V0 / om;
                        if ~isfinite(A)
                            A = 1;
                        end
                    else
                        A = 1;
                    end
                    yMin0 = -1.2 * abs(A) - 1;
                    yMax0 =  1.2 * abs(A) + 1;

                    app.xMin = min(xMin0, x0 - 1);
                    app.xMax = max(xMax0, x0 + 1);
                    app.yMin = min(yMin0, y0 - 1);
                    app.yMax = max(yMax0, y0 + 1);

                case 'Steady uniform: u=V0, v=0'
                    xMin0 = -1;
                    xMax0 = x0 + V0 * max(app.TmaxConst, 0) + 1;
                    yMin0 = -4;
                    yMax0 = 4;

                    app.xMin = min(xMin0, x0 - 1);
                    app.xMax = max(xMax0, x0 + 1);
                    app.yMin = min(yMin0, y0 - 1);
                    app.yMax = max(yMax0, y0 + 1);

                case 'Steady cellular: u=cos(ωx), v=sin(ωy)'
                    K = max(abs(om), 0.5);
                    L = 2 * pi / K;
                    app.xMin = x0 - L;
                    app.xMax = x0 + L;
                    app.yMin = y0 - L;
                    app.yMax = y0 + L;

                case 'Steady polynomial: u=(x^2+x-2)(x^2+x-12), v=y-3x+7'
                    app.xMin = -6;
                    app.xMax = 6;
                    app.yMin = -22;
                    app.yMax = 6;

                    app.xMin = min(app.xMin, x0 - 1);
                    app.xMax = max(app.xMax, x0 + 1);
                    app.yMin = min(app.yMin, y0 - 1);
                    app.yMax = max(app.yMax, y0 + 1);

                otherwise
                    app.xMin = x0 - 4;
                    app.xMax = x0 + 4;
                    app.yMin = y0 - 4;
                    app.yMax = y0 + 4;
            end
        end

        %% -------------------- Plot setup --------------------
        function resetAll(app)
            % Full reset: stop timer, clear axes, recreate graphics, reset histories
            %
            % Notes:
            %   * This is the ONLY operation that resets time and traces

            if app.isResetting
                return;
            end
            app.isResetting = true;
            c = onCleanup(@() setResetFalse(app)); %#ok<NASGU>

            stopTimer(app);
            while app.isTicking
                pause(0.01);
            end

            % Delete all children (including HandleVisibility='off')
            legend(app.Axes, 'off');

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

            app.streakX = x0;
            app.streakY = y0;
            app.lastReleaseT = 0;

            configureAxesAppearance(app);
            hold(app.Axes, 'on');

            computeAxisLimits(app);
            axis(app.Axes, [app.xMin, app.xMax, app.yMin, app.yMax]);

            % Quiver grid
            [app.xg, app.yg] = meshgrid( ...
                linspace(app.xMin, app.xMax, 14), ...
                linspace(app.yMin, app.yMax, 10));
            [u, v] = app.velocity(app.xg, app.yg, app.t);

            app.qh = quiver(app.Axes, app.xg, app.yg, u, v, 'AutoScale', 'on');
            app.qh.DisplayName = 'Velocity field';
            app.qh.HandleVisibility = 'on';

            % Source / release point
            app.srcDot = plot(app.Axes, x0, y0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
            app.srcDot.HandleVisibility = 'off';

            % Streamlines setup
            app.Cvals = linspace(app.yMin, app.yMax, app.numStream);

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
                return;
            end

            app.streamH = gobjects(app.numStream, 1);
            for i = 1:app.numStream
                app.streamH(i) = plot(app.Axes, nan, nan, '-', 'LineWidth', 1.2);

                if i == 1
                    app.streamH(i).DisplayName = 'Streamlines';
                    app.streamH(i).HandleVisibility = 'on';
                else
                    app.streamH(i).HandleVisibility = 'off';
                end

            end

        end

        function ensurePathObjects(app)
            % Create pathline polyline and particle marker
            if ~isempty(app.pathLine) && isgraphics(app.pathLine) && ...
               ~isempty(app.pDot) && isgraphics(app.pDot)
                return;
            end

            app.pathLine = plot(app.Axes, app.pathX, app.pathY, '-', 'LineWidth', 2.2, 'DisplayName', 'Pathline');
            app.pDot = plot(app.Axes, app.pathX(end), app.pathY(end), 'o', ...
                'MarkerFaceColor', app.PathColorPicker.Value, 'MarkerSize', 7);
            app.pDot.HandleVisibility = 'off';
        end

        function ensureStreakObjects(app)
            % Create streakline polyline and scatter points
            if ~isempty(app.streakLine) && isgraphics(app.streakLine) && ...
               ~isempty(app.streakDots) && isgraphics(app.streakDots)
                return;
            end

            app.streakLine = plot(app.Axes, app.streakX, app.streakY, '-', 'LineWidth', 2.2, 'DisplayName', 'Streakline');
            app.streakDots = scatter(app.Axes, app.streakX, app.streakY, 22, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);
            app.streakDots.HandleVisibility = 'off';
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
                return;
            end

            % Definitions
            n = numel(app.streamH);
            val = app.ProblemDropDown.Value;
            x = linspace(app.xMin, app.xMax, 1500);

            switch val

                % Uniform cases: y = C + (v/u) x
                case { ...
                    'Unsteady uniform: u=V0, v=V0 sin(ωt)', ...
                    'Steady uniform: u=V0, v=0' }

                    [u0, v0] = app.velocity(0, 0, app.t);
                    slope = 0;
                    if abs(u0) > 1e-12
                        slope = v0 / u0;
                    end

                    xRef = x(1);
                    C = linspace(app.yMin, app.yMax, n);
                    for i = 1:n
                        y = C(i) + slope * (x - xRef);
                        setLine(i, x, y);
                    end

                % Source: u=x, v=y -> y = C x (rays)
                case 'Steady source: u=x, v=y'
                    theta = linspace(0, pi, n);
                    R = 1.5 * max(app.xMax - app.xMin, app.yMax - app.yMin);
                    r = linspace(-R, R, 1200);
                    for i = 1:n
                        xi = r * cos(theta(i));
                        yi = r * sin(theta(i));
                        setLine(i, xi, yi);
                    end

                % Saddle: u=x, v=-y -> x y = C, plus axes
                case 'Steady saddle: u=x, v=-y'
                    if n >= 1
                        setLine(1, x, 0 * x); % y = 0
                    end
                    
                    if n >= 2
                        yAxis = linspace(app.yMin, app.yMax, 800);
                        setLine(2, 0 * yAxis, yAxis); % x = 0
                    end

                    m = max(n - 2, 0);
                    if m > 0
                        Cmax = 0.35 * (max(abs([app.xMin, app.xMax])) * max(abs([app.yMin, app.yMax])) + 1);
                        Cvals = linspace(-Cmax, Cmax, m);
                        Cvals(abs(Cvals) < 1e-9) = 0.15 * Cmax;

                        xx = x;
                        epsx = 1e-6 * max(1, max(abs(x)));
                        idx0 = abs(xx) < epsx;
                        xx(idx0) = epsx;

                        for k = 1:m
                            i = k + 2;
                            y = Cvals(k) ./ xx;
                            setLine(i, xx, y);
                        end
                    end

                % Rotation: u=y, v=-x -> circles
                case 'Steady rotation: u=y, v=-x'
                    rMax = 0.95 * min(max(abs([app.xMin, app.xMax])), max(abs([app.yMin, app.yMax])));
                    rMin = max(0.15 * rMax, 0.2);
                    rVals = linspace(rMin, rMax, n);
                    th = linspace(0, 2 * pi, 600);
                    for i = 1:n
                        xi = rVals(i) * cos(th);
                        yi = rVals(i) * sin(th);
                        setLine(i, xi, yi);
                    end

                % Shear: u=V0, v=V0 x -> dy/dx = x -> y = 0.5 x^2 + C
                case 'Steady shear: u=V0, v=V0 x'
                    C = linspace(app.yMin, app.yMax, n);
                    for i = 1:n
                        y = 0.5 * x .^ 2 + C(i);
                        setLine(i, x, y);
                    end

                % Cellular: u=cos(ωx), v=sin(ωy)
                % Implicit analytic: tan(ωy/2) = C (sec(ωx) + tan(ωx))
                case 'Steady cellular: u=cos(ωx), v=sin(ωy)'
                    om = app.OmegaField.Value;

                    if abs(om) < 1e-12

                        for i = 1:n
                            setLine(i, x, app.Cvals(min(i, numel(app.Cvals))) + 0 * x);
                        end

                        return;
                    end

                    y0Targets = linspace(app.yMin, app.yMax, n);
                    Cvals = tan(om * y0Targets / 2);

                    wx = om * x;
                    cwx = cos(wx);
                    f = (1 ./ cwx) + tan(wx);

                    bad = abs(cwx) < 0.03;
                    f(bad) = NaN;

                    for i = 1:n
                        y = (2 / om) * atan(Cvals(i) * f);
                        setLine(i, x, y);
                    end

                % Fallback: numeric streamline integration (steady non-uniform)
                otherwise
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

        %% -------------------- Hide/Show (no reset) --------------------
        function toggleVisibility(app)
            % Toggle object visibility from checkboxes (no reset)

            setVisible(app.qh, app.VelCheckBox.Value);
            setVisible(app.streamH, app.StreamCheckBox.Value);
            setVisible(app.pathLine, app.PathCheckBox.Value);
            setVisible(app.pDot, app.PathCheckBox.Value);
            setVisible(app.streakLine, app.StreakCheckBox.Value);
            setVisible(app.streakDots, app.StreakCheckBox.Value);

            
            updateLegend(app);

            function setVisible(h, tf)
                if isempty(h)
                    return;
                end

                if iscell(h)
                    h = [h{:}];
                end

                for k = 1:numel(h)

                    if ~isgraphics(h(k))
                        continue;
                    end

                    if tf
                        h(k).Visible = 'on';
                    else
                        h(k).Visible = 'off';
                    end

                end

            end
        end

        function updateTitle(app)
            % Update title with current time and velocity info      
            timeLine = sprintf('$t = %.2f\\,\\mathrm{s}$', app.t);
            app.Axes.Title.String = timeLine;
        end

        function updateLegend(app)
            % Update legend (optional).
            %
            % Notes:
            %   * Controlled by app.showLegend.Value

            if ~app.showLegend
                return
            end

            hs = gobjects(0);
            names = {};

            if ~isempty(app.qh) && isgraphics(app.qh) && strcmp(app.qh.Visible, 'on')
                hs(end + 1) = app.qh; 
                names{end + 1} = 'Velocity field';
            end

            if ~isempty(app.streamH) && isgraphics(app.streamH(1)) && strcmp(app.streamH(1).Visible, 'on')
                hs(end + 1) = app.streamH(1);
                names{end + 1} = 'Streamlines';
            end

            if ~isempty(app.pathLine) && isgraphics(app.pathLine) && strcmp(app.pathLine.Visible, 'on')
                hs(end + 1) = app.pathLine;
                names{end + 1} = 'Pathline';
            end

            if ~isempty(app.streakLine) && isgraphics(app.streakLine) && strcmp(app.streakLine.Visible, 'on')
                hs(end + 1) = app.streakLine;
                names{end + 1} = 'Streakline';
            end

            if isempty(hs)
                legend(app.Axes, 'off');
            else
                lg = legend(app.Axes, hs, names, 'Location', 'northwest', 'Box', 'on');
                lg.AutoUpdate = 'off';
                lg.Interpreter = 'latex';
            end
        end

        %% -------------------- Update per frame --------------------
        function updatePlot(app, dt)
            % Advance visualization by one frame.
            %
            % Args:
            %     dt (double): time step used for advection (0 for a "refresh")

            if app.isResetting || ~isgraphics(app.Axes)
                return;
            end

            % Update quiver
            if ~isempty(app.qh) && isgraphics(app.qh)
                [u, v] = app.velocity(app.xg, app.yg, app.t);
                set(app.qh, 'UData', u, 'VData', v);
            end

            % Update streamlines only if unsteady
            if app.isUnsteady()
                updateStreamlines(app);
            end

            if dt > 0
                app.tickCount = app.tickCount + 1;

                t0 = app.t - dt;
                [x0, y0] = getX0Y0(app);

                % Pathline: advance single particle with RK4, append to path every pathStride ticks
                [xn, yn] = rk4Particle(app, app.pathX(end), app.pathY(end), t0, dt);

                if mod(double(app.tickCount), app.pathStride) == 0
                    app.pathX(end + 1) = xn;
                    app.pathY(end + 1) = yn;

                    np = numel(app.pathX);
                    if np > app.maxPathPts
                        app.pathX = app.pathX(np - app.maxPathPts + 1:np);
                        app.pathY = app.pathY(np - app.maxPathPts + 1:np);
                    end
                else
                    % keep marker moving even if not storing this tick
                    app.pathX(end) = xn;
                    app.pathY(end) = yn;
                end

                % Streakline: release new particles at fixed time intervals, then advect all
                while (app.t - app.lastReleaseT) >= app.dtRel
                    app.streakX(end + 1) = x0;
                    app.streakY(end + 1) = y0;
                    app.lastReleaseT = app.lastReleaseT + app.dtRel;

                    ns = numel(app.streakX);
                    if ns > app.maxStreakPts
                        app.streakX = app.streakX(ns - app.maxStreakPts + 1:ns);
                        app.streakY = app.streakY(ns - app.maxStreakPts + 1:ns);
                    end

                end

                [app.streakX, app.streakY] = rk4Particle(app, app.streakX, app.streakY, t0, dt);
            end

            % Push data to graphics
            if ~isempty(app.pathLine) && isgraphics(app.pathLine)
                set(app.pathLine, 'XData', app.pathX, 'YData', app.pathY);
            end

            if ~isempty(app.pDot) && isgraphics(app.pDot)
                set(app.pDot, 'XData', app.pathX(end), 'YData', app.pathY(end));
            end

            if ~isempty(app.streakLine) && isgraphics(app.streakLine)
                set(app.streakLine, 'XData', app.streakX, 'YData', app.streakY);
            end

            if ~isempty(app.streakDots) && isgraphics(app.streakDots)
                set(app.streakDots, 'XData', app.streakX, 'YData', app.streakY);
            end

            % Update title
            updateTitle(app);

            % Update graphics
            drawnow limitrate nocallbacks
        end

        %% -------------------- Timer --------------------
        function onRun(app)
            % Run button callback
            startTimer(app);
        end

        function onStop(app)
            % Stop button callback
            stopTimer(app);
        end

        function startTimer(app)
            % Start (or restart) the fixed-rate timer
            %
            % Notes:
            %   * Period is quantized to milliseconds to avoid MATLAB timer warnings
            
            stopTimer(app);
            app.isRunning = true;

            dt = 1 / app.fpsConst;
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
                try
                    stop(app.Tmr);
                catch
                end

                if app.isTicking
                    return;
                end
                try
                    delete(app.Tmr);
                catch
                end

            end

            app.Tmr = timer.empty;
        end

        function onTick(app, dt)
            % Timer tick: advance time and update plots
            %
            % Args:
            %     dt (double): timer period (time increment)

            if app.isTicking
                return;
            end

            app.isTicking = true;
            c = onCleanup(@() setTickFalse(app));

            if app.isResetting || ~app.isRunning || ~isgraphics(app.Axes)
                return;
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

    end
end
