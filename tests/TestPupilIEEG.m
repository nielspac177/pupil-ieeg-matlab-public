classdef TestPupilIEEG < matlab.unittest.TestCase
    methods (Test)
        function loadsAuditedDerivedDataset(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            cfg = default_config(repoRoot);
            [tables, audit] = phg.loadDerivedTables(cfg);
            testCase.verifyEqual(height(tables.channel), 913);
            testCase.verifyEqual(numel(unique(string(tables.channel.PtID))), 18);
            selected = audit.value(audit.metric == "legacy_selected_channels");
            testCase.verifyEqual(selected, 177);
            testCase.verifyTrue(all(isfinite(tables.channel.XYZMNI), 'all'));
        end

        function storedResponseIsSemicontinuous(testCase)
            % The whole analysis design rests on this: the stored signed
            % response is exactly zero precisely when no suprathreshold
            % excursion was found. If that stopped being true, the hurdle
            % decomposition in runDerivedAnalyses would be modelling the wrong
            % thing and a Gaussian model might become appropriate again.
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            cfg = default_config(repoRoot);
            tables = phg.loadDerivedTables(cfg);
            isZeroArea = tables.channel.RespAreaNet == 0;
            isZeroContiguity = tables.channel.RespSig == 0;
            testCase.verifyEqual(isZeroArea, isZeroContiguity);
            testCase.verifyEqual(sum(~isZeroArea), 285);
        end

        function adjustsPValuesWithBenjaminiHochberg(testCase)
            p = [0.001; 0.008; 0.039; 0.041; 0.042; 0.60];
            q = phg.benjaminiHochberg(p);
            testCase.verifyEqual(q(1), 0.006, 'AbsTol', 1e-12);
            testCase.verifyEqual(q(2), 0.024, 'AbsTol', 1e-12);
            % Monotonicity must be enforced from the largest value downward.
            testCase.verifyTrue(all(diff(q) >= -1e-12));
            testCase.verifyTrue(all(q <= 1));
            testCase.verifyTrue(all(q >= p - 1e-12));
            withMissing = phg.benjaminiHochberg([0.01; NaN; 0.02]);
            testCase.verifyTrue(isnan(withMissing(2)));
        end

        function parsesElectrodeShaftFromContactLabel(testCase)
            labels = ["LHIP1"; "LHIP12"; "RAMG3"; "LPOST32"; "  LOFC7 "; ""];
            expected = ["LHIP"; "LHIP"; "RAMG"; "LPOST"; "LOFC"; "UnknownLead"];
            testCase.verifyEqual(phg.parseLeadLabel(labels), expected);
        end

        function detectsOverlappingTextInFigureAudit(testCase)
            % The previous audit could not see label collisions at all. This
            % builds a figure with two deliberately coincident labels and
            % checks that the audit reports the overlap.
            fig = figure('Visible', 'off', 'Color', 'w', ...
                'Units', 'centimeters', 'Position', [1 1 8 6]);
            cleanup = onCleanup(@() close(fig));
            ax = axes(fig); %#ok<LAXES>
            plot(ax, 1:10);
            text(ax, 0.5, 0.5, 'COLLIDING LABEL', 'Units', 'normalized');
            text(ax, 0.5, 0.5, 'COLLIDING LABEL', 'Units', 'normalized');
            imagePath = string(tempname) + ".png";
            imageCleanup = onCleanup(@() deleteIfPresent(imagePath));
            print(fig, char(imagePath), '-dpng', '-r100');
            qc = phg.auditFigureLayout(fig, imagePath, "overlapTest");
            testCase.verifyGreaterThanOrEqual(qc.text_overlap_count, 1);
            testCase.verifyFalse(qc.layout_ok);
            clear imageCleanup cleanup
        end

        function detectsSyntheticRipple(testCase)
            rng(20260802, 'twister');
            fs = 1000;
            time = (0:1/fs:5-1/fs)';
            signal = 2 .* randn(size(time));
            burst = time >= 2 & time < 2.10;
            signal(burst) = signal(burst) + ...
                20 .* sin(2*pi*100*time(burst));
            events = phg.detectRipples(signal, fs, ...
                'maximumZ', 100, 'minimumZ', 2.5);
            accepted = events(events.accepted, :);
            testCase.verifyGreaterThanOrEqual(height(accepted), 1);
            overlaps = accepted.onset_seconds <= 2.10 & ...
                accepted.offset_seconds >= 2.00;
            testCase.verifyTrue(any(overlaps));
            testCase.verifyEqual(accepted.peak_frequency_hz(find(overlaps,1)), ...
                100, 'AbsTol', 8);
        end

        function assignsExposureBalancedPupilStates(testCase)
            pupilTime = (0:0.01:60)';
            pupilSize = sin(2*pi*pupilTime/20) + pupilTime/600;
            events = table([5;25;45], [5.05;25.05;45.05], ...
                'VariableNames', {'onset_seconds', 'offset_seconds'});
            [assigned, exposure] = phg.assignPupilStates( ...
                events, pupilTime, pupilSize, 6);
            testCase.verifyTrue(all(isfinite(assigned.pupil_state)));
            testCase.verifyEqual(height(exposure), 6);
            testCase.verifyLessThan(max(exposure.duration_seconds) - ...
                min(exposure.duration_seconds), 0.1);
        end

        function buildsPrivateBidsScaffold(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            cfg = default_config(repoRoot);
            temporaryRoot = string(tempname);
            mkdir(temporaryRoot);
            cleanup = onCleanup(@() rmdir(temporaryRoot, 's'));
            cfg.bidsRoot = fullfile(temporaryRoot, 'bids');
            phg.createBidsScaffold(cfg);
            datasetDescription = fullfile(cfg.bidsRoot, ...
                'dataset_description.json');
            testCase.verifyTrue(isfile(datasetDescription));
            decoded = jsondecode(fileread(datasetDescription));
            testCase.verifyEqual(string(decoded.BIDSVersion), cfg.bidsVersion);
            testCase.verifyTrue(isfile(fullfile(cfg.bidsRoot, ...
                'participants.tsv')));
            testCase.verifyTrue(isfile(fullfile(cfg.bidsRoot, ...
                'derivatives', 'pupil-ieeg', 'dataset_description.json')));
            clear cleanup
        end

        function buildsEEGNetWhenToolboxIsAvailable(testCase)
            testCase.assumeTrue(exist('dlnetwork', 'file') == 2 && ...
                exist('groupedConvolution2dLayer', 'file') == 2, ...
                'Deep Learning Toolbox is not installed.');
            network = phg.buildEEGNetRegressor(8, 256, 64);
            testCase.verifyClass(network, 'dlnetwork');
            names = string({network.Layers.Name});
            testCase.verifyTrue(any(names == "spatial_depthwise"));
            testCase.verifyTrue(any(names == "separable_depthwise"));
        end

        function computesKnownSha256(testCase)
            temporaryFile = string(tempname);
            cleanup = onCleanup(@() deleteIfPresent(temporaryFile));
            phg.writeText(temporaryFile, "abc", true);
            testCase.verifyEqual(phg.sha256File(temporaryFile), ...
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
            clear cleanup
        end
    end
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
