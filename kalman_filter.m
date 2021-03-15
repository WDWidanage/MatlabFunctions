% Implemnting a Kalman filter
%
% Read, understood and now to test
%
% W.D. Widanage 22/06/11 (Disconnected)

function [x_post, P_post]=kalman_filter(Model,z,R,Q,x_p,P_p,u)

% Inputs:
% Model - structure variable with intial values for each of the state
% space matrices
%       Model.A - A matrix size na x na
%       Model.B - B matrix size na x 1
%       Model.C - C matrix size 1 x na
% z - Measured/Observed state at time instant t
% R - Process noise variance at time instant t
% Q - Measured noise variance at time instant t
% x_p - Estimate of state for prediction
% P_p - Estimate of prediction error covaraince
% u - Control/Input at time instant t 
%
% Outputs:
% x_post - Posteriori estimate of state after correction
% P_post - Posteriori estimate of preiction error covariance


% Prediction/Time update step, priori estimates
x_pre = Model.A*x_p + Model.B*u; 
P_pre = Model.A*P_p*Model.A'+R;

% Correction step, posteriori estimates
%Kalman gain
K = P_pre*Model.C'*inv(Model.C*P_pre*Model.C'+Q);
x_post = x_pre+K*(z-Model.C*x_pre);
P_post = (eye(size(K*Model.C))-K*Model.C)*P_pre;


