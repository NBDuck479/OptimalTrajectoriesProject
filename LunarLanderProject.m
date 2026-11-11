% Lunar Lander Project 

% Lunar Constants
moonRadius = 1736; % km 
moonMu     = 4902.8; % km^3/s^2

% Define Low Lunar Polar Orbit 
altitudeLLO = 100; % km 
radiusLLO   = moonRadius + altitudeLLO; 
orbitPeriod = 2 * pi * sqrt(radiusLLO^3/moonMu);
timeSpan    = linspace(0, orbitPeriod, 1000);

% Define Rocket parameters for lander 
engineParams.Isp  = 300; % Sec 
engineParams.Ve   = engineParams.Isp * 9.80665/1000; 
gMoon             = 1.62/1000;
initialMass       = 10000; % kg
Thrust2Weight     = 2; 
Tmax              = initialMass * gMoon * Thrust2Weight;  
engineParams.bMax = Tmax / engineParams.Ve ; 

 
% initial state vector in MCI 
stateLL0 = [0 ; 0 ; radiusLLO; -sqrt(moonMu/radiusLLO); 0; 0]; 

% Propagate lunar obit 
twoBodyDynamics = @(t, state) [state(4:6); -moonMu*state(1:3)/norm(state(1:3))^3];
[time, stateHistory] = ode45(twoBodyDynamics, timeSpan, stateLL0);

% Choose initial condition from LLO to descent from 
idxState = 400; 
initialDescendState = stateHistory(idxState,:)'; 

% Landing site in MCI 
moonLandingSite = [0;0;-moonRadius]; 

% Choose final conditions for landing 
finalDescendState = [moonLandingSite; 0; 0; 0]; 

% Final time also constrained 
tf = orbitPeriod/5; % Seconds 

% Initialize G value to something large
G = ones(6,1); 

% Initial adjoint state guess 
adjointStates = 1e-5 * [1;1;1;1;1;1;1];

% Iteration counter 
iter = 0; 

% Integratino tolerance
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% If G is not small enough then continue iterating
while norm(G) > .1 && iter < 30

    % Propagate state and adjoints with ode45
    initialDescentStateWithAdjoints = [initialDescendState; initialMass; adjointStates];
    descentTimeSpan = linspace(0, tf, 1000);
    [descentTime, descentHistory] = ode45(@(t,X) stateAndAdjointProp(t,X,moonMu,engineParams), ...
        descentTimeSpan, initialDescentStateWithAdjoints, opts);

    % Determine how well my constraints were met
    Gdesired =  finalDescendState;
    Gachieved = descentHistory(end,1:6)';

    G = Gdesired - Gachieved;

    % Need partial of G wrt P0 with forward finitie difference 
    deltaP0 = 1e-8; 

    % Loop each of the 7 adjoint states 
    for j = 1:7
        Adjointpert = adjointStates;
        Adjointpert(j) = adjointStates(j) + deltaP0;

        % Intergate with this pert adjoint 
        perturbedInitialState = [initialDescendState; initialMass; Adjointpert];

        [~, perturbedHistory] = ode45(@(t,X) stateAndAdjointProp(t,X,moonMu,engineParams), ...
            descentTimeSpan, perturbedInitialState, opts);

        Gpert = finalDescendState - perturbedHistory(end,1:6)';

        GwrtP0(:,j) = (Gpert - G) / deltaP0;
    end

    % Update P0 and run it again
    adjointStates = adjointStates - .5 * (GwrtP0\G); 

    iter = iter + 1; 

end

% Primer vector gives the control history throughout the trajecotry, u =
% -Pv 
u = - descentHistory(:, 10:12); 



% Plot orbit with moon 
figure;
plot3(stateHistory(1:idxState,1), stateHistory(1:idxState,2), stateHistory(1:idxState,3), 'b', 'LineWidth', 1.5);
hold on;
% Plot the descent 
plot3(descentHistory(:,1), descentHistory(:,2), descentHistory(:,3))

% quiver3(stateHistory(1:10:idxState,1), stateHistory(1:10:idxState,2), stateHistory(1:10:idxState,3), ...
%     stateHistory(1:10:idxState,4), stateHistory(1:10:idxState,5), stateHistory(1:10:idxState,6), ...
%     1, 'r', 'LineWidth', 1);

% Plot the control input for descent 
quiver3(descentHistory(1:10:end,1), descentHistory(1:10:end,2), descentHistory(1:10:end,3), ...
    u(1:10:end,1), u(1:10:end,2), u(1:10:end,3), 5, 'r', 'LineWidth', 1); 
plot3(moonLandingSite(1), moonLandingSite(2), moonLandingSite(3), 'r.', 'MarkerSize', 25);
[xMoon, yMoon, zMoon] = sphere(40);
surf(moonRadius*xMoon, moonRadius*yMoon, moonRadius*zMoon, ...
    'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
axis equal;
grid on;
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
title('Low Lunar Polar Orbit');



function [XDot] = stateAndAdjointProp(t, X, moonMu, engineParams)
% Set up integration for pos/vel and adjoints given PMP result 

% Set up engine parameters 
bMax = engineParams.bMax; 
Ve   = engineParams.Ve; 

% States
r = X(1:3); 
v = X(4:6); 
m = X(7); 

% Adjoints 
Pr = X(8:10);
Pv = X(11:13);
Pm = X(14); 

% Determine switching function 
Kappa = Ve/m * Pv + Pm;

if Kappa < 0
    % Minimizing the Hamiltonian will be no thrust
    b = 0;

elseif Kappa > 0
    % Minimizing th hamiltonian will be max thrust (mass flow)
    b = bMax;

else
    error("idk what to do Kappa == 0");
end

% Math 
rNorm = norm(r);
rHat = r/rNorm;

% Control direction 
uDir = -Pv/norm(Pv); 

% State Dynamics 
xDot = v; 
vDot = b*Ve/m * uDir - moonMu/(rNorm^3) * r; 
mDot = -b; 

% Adjoint Dynamics
PrDot = -moonMu/(rNorm^3) * (3*rHat*rHat' - eye(3)) * Pv; 
PvDot = -Pr; 
PmDot = -b*Ve/m^2 * norm(Pv); 

% Build up state derivative
XDot = [xDot; vDot; mDot; PrDot; PvDot; PmDot];

end


