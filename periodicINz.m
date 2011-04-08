function [Et,Hr,Hz] = periodicINz(Et,Hr,Hz)

% lower boundaries
Et(1,:) = Et(end-1,:);
Hr(1,:) = Hr(end-1,:);
Hz(1,:) = Hz(end-1,:);

% upper boundaries
Et(end,:) = Et(2,:);
Hr(end,:) = Hr(2,:);
Hz(end,:) = Hz(2,:);