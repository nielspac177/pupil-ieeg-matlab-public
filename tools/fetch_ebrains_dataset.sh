#!/usr/bin/env bash
# Resumable download of the EBRAINS pupillometry dataset.
#
# The archive is large and the connection will drop. This script is built to be
# re-run: every file is fetched with `wget -c`, which resumes a partial file
# from wherever it stopped rather than starting over, and files already present
# at full size are skipped. Killing it and running it again is the intended
# recovery path, not a failure mode.
#
# Requires an EBRAINS access token. Get one while logged in at
#   https://data-proxy.ebrains.eu/api/v1/buckets/<dataset-id>
# or from the "Get token" control in the EBRAINS Knowledge Graph UI, then:
#
#   export EBRAINS_TOKEN='eyJ...'
#   bash tools/fetch_ebrains_dataset.sh
#
# Tokens are short-lived. When one expires the script stops with a clear
# message; export a fresh token and re-run, and it picks up where it left off.

set -uo pipefail

DATASET_ID="${DATASET_ID:-f2f41456-82d6-40cc-88e2-12a8d5b375f8}"
TARGET="${TARGET:-/Volumes/Niels_3/ebrains_pupil_ieeg}"
API="https://data-proxy.ebrains.eu/api/v1/buckets"

if [[ -z "${EBRAINS_TOKEN:-}" ]]; then
    echo "EBRAINS_TOKEN is not set. See the header of this script." >&2
    exit 1
fi

mkdir -p "$TARGET"
cd "$TARGET" || exit 1

echo "Target:    $TARGET"
echo "Free space: $(df -h . | awk 'NR==2 {print $4}')"
echo

# ---------------------------------------------------------------- file listing
# The bucket API pages through objects with a marker. The listing is written to
# disk so a re-run does not need to re-enumerate, and so you can inspect what
# the script believes it should be fetching.
LISTING="$TARGET/_file_listing.json"
if [[ ! -s "$LISTING" ]]; then
    echo "Listing bucket contents..."
    : > "$LISTING"
    marker=""
    while :; do
        page=$(curl -sf -H "Authorization: Bearer $EBRAINS_TOKEN" \
            "$API/$DATASET_ID?limit=1000${marker:+&marker=$marker}") || {
            echo "Listing failed. Token expired, or access not yet granted." >&2
            rm -f "$LISTING"
            exit 1
        }
        count=$(printf '%s' "$page" | python3 -c \
            'import json,sys; print(len(json.load(sys.stdin).get("objects",[])))')
        [[ "$count" -eq 0 ]] && break
        printf '%s\n' "$page" >> "$LISTING"
        marker=$(printf '%s' "$page" | python3 -c \
            'import json,sys; o=json.load(sys.stdin)["objects"]; print(o[-1]["name"])')
        echo "  ... $count objects (through $marker)"
    done
fi

mapfile -t FILES < <(python3 - "$LISTING" <<'PY'
import json, sys
names = []
for line in open(sys.argv[1]):
    line = line.strip()
    if not line:
        continue
    for obj in json.loads(line).get("objects", []):
        names.append(f"{obj['name']}\t{obj.get('bytes', 0)}")
print("\n".join(names))
PY
)

echo
echo "${#FILES[@]} files to consider."
echo

# ------------------------------------------------------------------- download
downloaded=0
skipped=0
failed=0

for entry in "${FILES[@]}"; do
    name="${entry%%$'\t'*}"
    size="${entry##*$'\t'}"

    if [[ -f "$name" ]]; then
        have=$(stat -f%z "$name" 2>/dev/null || stat -c%s "$name" 2>/dev/null || echo 0)
        if [[ "$have" == "$size" ]]; then
            skipped=$((skipped + 1))
            continue
        fi
    fi

    mkdir -p "$(dirname "$name")"

    # Ask the proxy for a temporary direct URL rather than streaming through it,
    # because a redirect-following download cannot be resumed reliably.
    url=$(curl -sf -H "Authorization: Bearer $EBRAINS_TOKEN" \
        "$API/$DATASET_ID/${name}?redirect=false" | python3 -c \
        'import json,sys; print(json.load(sys.stdin).get("url",""))')

    if [[ -z "$url" ]]; then
        echo "  ! no URL for $name (token expired?)" >&2
        failed=$((failed + 1))
        continue
    fi

    echo "→ $name"
    if wget -c -q --show-progress --tries=20 --waitretry=15 \
            --read-timeout=120 -O "$name" "$url"; then
        downloaded=$((downloaded + 1))
    else
        echo "  ! failed: $name (re-run to resume)" >&2
        failed=$((failed + 1))
    fi
done

echo
echo "Downloaded $downloaded, already complete $skipped, failed $failed."
if [[ "$failed" -gt 0 ]]; then
    echo "Re-run this script to resume the failed files."
    exit 1
fi
echo "Complete. Analysis protocol: docs/replication_plan_ebrains.md"
