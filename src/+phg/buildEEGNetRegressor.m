function network = buildEEGNetRegressor(nChannels, nTime, nTargetTime, options)
%BUILDEEGNETREGRESSOR Compact EEGNet-style pupil time-course regressor.

arguments
    nChannels (1,1) double {mustBeInteger,mustBePositive}
    nTime (1,1) double {mustBeInteger,mustBePositive}
    nTargetTime (1,1) double {mustBeInteger,mustBePositive}
    options.temporalFilters (1,1) double {mustBeInteger,mustBePositive} = 8
    options.depthMultiplier (1,1) double {mustBeInteger,mustBePositive} = 2
    options.dropout (1,1) double {mustBeGreaterThanOrEqual(options.dropout,0),mustBeLessThan(options.dropout,1)} = 0.5
end

F1 = options.temporalFilters;
D = options.depthMultiplier;
F2 = F1 * D;
temporalKernel = min(nTime, max(8, round(nTime / 8)));
separableKernel = min(nTime, 16);

layers = [
    imageInputLayer([nChannels nTime 1], 'Normalization', 'none', ...
        'Name', 'ieeg_input')
    convolution2dLayer([1 temporalKernel], F1, 'Padding', 'same', ...
        'BiasLearnRateFactor', 0, 'Name', 'temporal_convolution')
    batchNormalizationLayer('Name', 'temporal_batchnorm')
    groupedConvolution2dLayer([nChannels 1], F2, F1, ...
        'Name', 'spatial_depthwise')
    batchNormalizationLayer('Name', 'spatial_batchnorm')
    eluLayer('Name', 'spatial_elu')
    averagePooling2dLayer([1 4], 'Stride', [1 4], ...
        'Name', 'pool_1')
    dropoutLayer(options.dropout, 'Name', 'dropout_1')
    groupedConvolution2dLayer([1 separableKernel], F2, F2, ...
        'Padding', 'same', 'Name', 'separable_depthwise')
    convolution2dLayer([1 1], F2, 'Name', 'separable_pointwise')
    batchNormalizationLayer('Name', 'separable_batchnorm')
    eluLayer('Name', 'separable_elu')
    averagePooling2dLayer([1 8], 'Stride', [1 8], ...
        'Name', 'pool_2')
    dropoutLayer(options.dropout, 'Name', 'dropout_2')
    flattenLayer('Name', 'flatten')
    fullyConnectedLayer(nTargetTime, 'Name', 'pupil_time_course')
    ];
network = dlnetwork(layers);
end
