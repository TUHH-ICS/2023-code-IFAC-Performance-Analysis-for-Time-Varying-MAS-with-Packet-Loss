%---------------------------------------------------------------------------------------------------
% For Paper
% "Robust Performance Analysis for Time-Varying Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

classdef ParforProgressMeter < handle
    %PARFORPROGRESSMETER Class to display a progress meter in the command line when using parfor
    %loops to parallelize computations.
    %   Tracking progress in Matlabs parfor loops is not straight forward, because one cannot
    %   increment a shared variable to indicate progress in the usual way. Instead, we need to make
    %   use of a parallel computing primitive to synchronize the threads and print the achieved
    %   progress. This is handled by this class in an encapsulated manner.

    properties(GetAccess = public, SetAccess = private)
        progress    % Tracks the number of completed tasks
    end

    properties(GetAccess = public, SetAccess = immutable)
        elements    % Total number of tasks
        interval    % Increments in which to print the progress
    end

    properties(GetAccess = private, SetAccess = private)
        timer       % Timer for the runtime, to predict the remaining time
        lastPrint   % Percentage of last status print
    end

    properties(GetAccess = private, SetAccess = immutable)
        queue       % Synchronized primitive for status update
    end

    methods
        function obj = ParforProgressMeter(elements, interval)
            %PARFORPROGRESSMETER Construct an instance of this class
            %   Constructing the progress meter does not start the counter. This is only done when
            %   explicitely calling start() or notify(i)
            %
            %   Arguments:
            %       elements -> Number of events until 100% progress
            %       interval -> [Optional] Interval of update, given in percent. Default: 0.05

            % Initialize numeric fields
            obj.elements  = elements;
            obj.progress  = 0;
            obj.timer     = [];
            obj.lastPrint = 0;

            % Handle optional argument
            if nargin <= 1
                obj.interval = 0.05; % Default: Update every 5%
            else
                obj.interval = interval;
            end

            % Initialize queue
            obj.queue = parallel.pool.DataQueue;
            afterEach(obj.queue, @(~)obj.update());
        end

        function start(obj)
            %START Start the timer for predicting the remaining time
            obj.timer = tic;
        end

        function notify(obj, i)
            %NOTIFY Increment the progress by one
            send(obj.queue, i)
        end
    end

    methods(Access = private)
        function update(obj)
            %UPDATE Increment progress and print status if required

            % If the timer has not been started yet, this is now done after the fact
            if isempty(obj.timer)
                obj.start()
            end

            % Update the progress
            obj.progress    = obj.progress + 1;
            percentFinished = obj.progress / obj.elements;

            % If the last update was sufficiently long ago, print a new status update
            if percentFinished - obj.lastPrint > obj.interval
                obj.lastPrint = percentFinished;
                timeElapsed   = toc(obj.timer);
                timeRemaining = (1-percentFinished)/percentFinished * timeElapsed;

                fprintf('Calculation %5.1f%% finished. Time remaining: %s\n', 100*percentFinished, format_duration(timeRemaining))
            end
        end
    end
end
