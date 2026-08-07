function loss = pupilSequenceLoss(prediction, target)
%PUPILSEQUENCELOSS NMSE of amplitude plus first temporal derivative.

amplitudeDenominator = mean((target - mean(target, 1)).^2, 'all') + eps;
amplitudeLoss = mean((prediction - target).^2, 'all') ./ amplitudeDenominator;

predictionDerivative = diff(prediction, 1, 1);
targetDerivative = diff(target, 1, 1);
derivativeDenominator = mean((targetDerivative - ...
    mean(targetDerivative, 1)).^2, 'all') + eps;
derivativeLoss = mean((predictionDerivative - targetDerivative).^2, 'all') ...
    ./ derivativeDenominator;
loss = amplitudeLoss + derivativeLoss;
end
