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

        function appliesPrespecifiedReplicationRule(testCase)
            % The replication verdict is the one number in this project most
            % exposed to motivated reasoning, so the rule is tested against
            % the wording of replication_plan_ebrains.md section 2 rather
            % than trusted. The case that matters is the fourth: an odds
            % ratio well below 1 whose interval crosses 1 is explicitly NOT
            % a replication, however tempting it looks.
            discovery = [0.021 0.207];
            verdict = @(o, l, h) phg.replicationVerdict(o, l, h, discovery);

            testCase.verifyEqual(verdict(0.065, 0.021, 0.207), "replication");
            testCase.verifyEqual(verdict(0.50, 0.30, 0.85), "replication");

            % OR inside the discovery interval, interval crosses 1.
            testCase.verifyEqual(verdict(0.10, 0.01, 1.40), ...
                "directionally_consistent_underpowered");
            testCase.verifyEqual(verdict(0.207, 0.02, 2.10), ...
                "directionally_consistent_underpowered");

            % Below 1 and crossing 1, but outside the discovery interval:
            % consistent in direction only, and not claimable as either.
            testCase.verifyEqual(verdict(0.60, 0.20, 1.80), "inconclusive");
            testCase.verifyEqual(verdict(0.010, 0.0001, 1.5), "inconclusive");

            % Failure in both of its forms.
            testCase.verifyEqual(verdict(1.00, 0.40, 2.50), "failure");
            testCase.verifyEqual(verdict(3.20, 1.40, 7.30), "failure");

            testCase.verifyEqual(verdict(NaN, NaN, NaN), "not_estimable");
            testCase.verifyEqual(verdict(0.05, NaN, Inf), "not_estimable");
        end

        function periPeakMeasurementRecoversKnownResponses(testCase)
            % The ported measurement must be sign-symmetric and must return
            % exactly zero when nothing is locked to the peak times, because
            % the hurdle model depends on that zero being structural.
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            cfg = default_config(repoRoot);
            rng(20260807, 'twister');

            fs = cfg.replication.pupilFs;
            t = (0:1 / fs:2000 - 1 / fs)';
            background = cumsum(randn(size(t))) / sqrt(fs);
            background = background - movmean(background, 60 * fs);
            peaks = (60:7:1940)';

            kernelTime = (-2:1 / fs:6)';
            kernel = 3 * exp(-((kernelTime - 1.5) / 1.2) .^ 2);
            dilating = background;
            for k = 1:numel(peaks)
                first = round((peaks(k) + kernelTime(1)) * fs) + 1;
                index = first:(first + numel(kernel) - 1);
                inside = index >= 1 & index <= numel(dilating);
                dilating(index(inside)) = dilating(index(inside)) + ...
                    kernel(inside);
            end

            dilation = phg.measurePeriPeakResponse(dilating, t, peaks, cfg);
            testCase.verifyGreaterThan(dilation.RespAreaNet, 0);
            testCase.verifyGreaterThan(dilation.RespSig, 0);

            null = phg.measurePeriPeakResponse(background, t, peaks, cfg);
            testCase.verifyEqual(null.RespAreaNet, 0);
            testCase.verifyEqual(null.RespSig, 0);

            % The identity the hurdle decomposition rests on.
            testCase.verifyEqual(null.RespAreaNet == 0, null.RespSig == 0);
            testCase.verifyEqual(dilation.RespAreaNet == 0, ...
                dilation.RespSig == 0);
        end
    end
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
