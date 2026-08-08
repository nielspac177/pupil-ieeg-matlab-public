function configurations = replicationConfigurations()
%REPLICATIONCONFIGURATIONS The four analyses fixed before any was run.
%   The primary is the gaze-regressed pupil under epoch rule A. The other three
%   are sensitivity analyses named in docs/replication_amendment_01.md sections
%   3 and 5, and are reported alongside the primary regardless of how they come
%   out. Defined in one place so that the staged runs and the assembler cannot
%   disagree about what was run.

configurations = [
    struct('label', "primary",   'variant', "gaze_regressed", 'rule', "A")
    struct('label', "raw_pupil", 'variant', "raw",            'rule', "A")
    struct('label', "rule_b",    'variant', "gaze_regressed", 'rule', "B")
    struct('label', "rule_c",    'variant', "gaze_regressed", 'rule', "C")
    ];
end
