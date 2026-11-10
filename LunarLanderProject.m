% Lunar Lander Project 

% Lunar Constants
moonRadius = 1736; % km 
moonMu     = 4902.8; % km^3/s^2

% Define Low Lunar Polar Orbit 
altitudeLLO = 100; % km 
radiusLLO   = moonRadius + altitudeLLO; 
orbitPeriod = 2 * pi * sqrt(radiusLLO^3/moonMu);
timeSpan    = linspace(0, orbitPeriod, 1000);

% initial state vector in MCI 
stateLL0 = [0 ; 0 ; radiusLLO; -sqrt(moonMu/radiusLLO); 0; 0]; 

% Propagate lunar obit 
twoBodyDynamics = @(t, state) [state(4:6); -moonMu*state(1:3)/norm(state(1:3))^3];
[time, stateHistory] = ode45(twoBodyDynamics, timeSpan, stateLL0);

% Choose initial condition from LLO to descent from 
initialDescendState = stateHistory(400,:)'; 

% Landing site in MCI 
moonLandingSite = [0;0;-moonRadius]; 

% Choose final conditions for landing 
finalDescendState = [moonLandingSite; 0; 0; 0]; 

% Final time also constrained 
tf = orbitPeriod/5; % Seconds 

% Initialize G value to something large
G = ones(6,1); 

% Initial adjoint state guess 
adjointStates = 1e-5 * [1;1;1;1;1;1];

% Iteration counter 
iter = 0; 

% Integratino tolerance
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% If G is not small enough then continue iterating
while norm(G) > .1 && iter < 30

    % Propagate state and adjoints with ode45
    initialDescentStateWithAdjoints = [initialDescendState; adjointStates];
    descentTimeSpan = linspace(0, tf, 1000);
    [descentTime, descentHistory] = ode45(@(t,X) stateAndAdjointProp(t,X,moonMu), ...
        descentTimeSpan, initialDescentStateWithAdjoints, opts);

    % Determine how well my constraints were met
    Gdesired =  finalDescendState;
    Gachieved = descentHistory(end,1:6)';

    G = Gdesired - Gachieved;

    % Need partial of G wrt P0 with forward finitie difference 
    deltaP0 = 1e-8; 

    % Loop each of the 6 adjoint states 
    for j = 1:6
        Adjointpert = adjointStates;
        Adjointpert(j) = adjointStates(j) + deltaP0;

        % Intergate with this pert adjoint 
        perturbedInitialState = [initialDescendState; Adjointpert];

        [~, perturbedHistory] = ode45(@(t,X) stateAndAdjointProp(t,X,moonMu), ...
            descentTimeSpan, perturbedInitialState, opts);

        Gpert = finalDescendState - perturbedHistory(end,1:6)';

        GwrtP0(:,j) = (Gpert - G) / deltaP0;
    end

    % Update P0 and run it again
    adjointStates = adjointStates - .5 * (GwrtP0\G); 

    iter = iter + 1; 

end



% Plot orbit with moon 
figure;
plot3(stateHistory(:,1), stateHistory(:,2), stateHistory(:,3), 'b', 'LineWidth', 1.5);
hold on;
% Plot the descent 
plot3(descentHistory(:,1), descentHistory(:,2), descentHistory(:,3))

quiver3(stateHistory(1:10:end,1), stateHistory(1:10:end,2), stateHistory(1:10:end,3), ...
    stateHistory(1:10:end,4), stateHistory(1:10:end,5), stateHistory(1:10:end,6), ...
    2, 'r', 'LineWidth', 1);
plot3(moonLandingSite(1), moonLandingSite(2), moonLandingSite(3), 'r.', 'MarkerSize', 25);
[xMoon, yMoon, zMoon] = sphere(40);
surf(moonRadius*xMoon, moonRadius*yMoon, moonRadius*zMoon, ...
    'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
axis equal;
grid on;
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
title('Low Lunar Polar Orbit');



function [XDot] = stateAndAdjointProp(t, X, moonMu)
% Set up integration for pos/vel and adjoints given PMP result 

% States
r = X(1:3); 
v = X(4:6); 

% Adjoints 
Pr = X(7:9);
Pv = X(10:12);

% Math 
rNorm = norm(r);
rHat = r/rNorm;

% State Dynamics 
xDot = v; 
vDot = -Pv - moonMu/(rNorm^3) * r; 
% Adjoint Dynamics
PrDot = -moonMu/(rNorm^3) * (3*rHat*rHat' - eye(3)) * Pv; 
PvDot = -Pr; 

% Build up state derivative
XDot = [xDot; vDot; PrDot; PvDot];

end


