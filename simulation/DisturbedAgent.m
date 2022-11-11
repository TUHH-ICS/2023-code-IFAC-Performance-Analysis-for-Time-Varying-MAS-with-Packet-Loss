%---------------------------------------------------------------------------------------------------
% For Paper
% "Robust Performance Analysis for Time-Varying Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

classdef DisturbedAgent < IntegratorAgent
    %DISTURBEDAGENT Simple LTI agent that can optionally be disturbed by a
    %discrete impulse. The agent tries to maintain a formation using a
    %distributed consensus protocol.

    properties(Constant, GetAccess = public)
        kappa = 0.1 % Gain of the consensus protocol
    end

    properties(GetAccess = public, SetAccess = private)
        disturbance % Disturbance signal of the agents
        pattern_idx % Index of the topology pattern
    end

    properties(GetAccess = public, SetAccess = immutable)
        ref       % Constant formation reference
        adjacency % Adjacency vectors
        pattern   % Switching pattern of the adjacency
    end

    properties(Dependent, GetAccess = public, SetAccess = immutable)
        neighbours % Current neighbours of the agent
    end

    methods
        function obj = DisturbedAgent(id, initialPos, ref, adjacency, pattern, disturbed)
            %DISTURBEDAGENT Construct an instance of this class
            %   Sets up the correct agent dynamics and initializes the
            %   agent to the given initial position.

            obj@IntegratorAgent(id, 1, initialPos);
            obj.ref         = ref;
            obj.adjacency   = adjacency;
            obj.pattern     = pattern;
            obj.pattern_idx = 1;

            % Set disturbance status
            obj.disturbance = zeros(obj.dim, 1);
            if disturbed
                obj.disturbance(disturbed) = 1;
            end
        end

        function value = get.neighbours(obj)
            %GET.NEIGHBOURS Implementation of the dependent neighbours property.
            %   Selects the appropriate column of the adjancy and finds the
            %   current neighbour IDs.

            index = obj.pattern(obj.pattern_idx);
            value = find(obj.adjacency(:, index));
        end

        function step(obj)
            % STEP Perform one step of simulation

            % Receive messages from the network
            messages = obj.receive();

            % Filter to only neighbours
            if ~isempty(messages)
                mask = any(obj.neighbours == [messages.sender], 1);
                messages = messages(mask);
            end

            if ~isempty(messages)
                % Calculate formation error and distributed controller
                % state. We use the standard Laplacian, therefore we sum
                % up the data.
                data = [messages.data];
                err  = sum(obj.position - [data.position] - obj.ref, 2);
            else
                err = zeros(obj.dim,1);
            end

            % Evaluate agent dynamics
            u = -DisturbedAgent.kappa * err;
            obj.move(u + obj.disturbance);

            % Shift to next communication pattern
            obj.pattern_idx = obj.pattern_idx + 1;
            if obj.pattern_idx > length(obj.pattern)
                obj.pattern_idx = 1;
            end

            % Apply disturbance only once
            obj.disturbance(:) = 0;

            % Send message to network, include only the position
            data = struct;
            data.position  = obj.position - obj.ref;
            obj.send(data)
        end
    end
end
