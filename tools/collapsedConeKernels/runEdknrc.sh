#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${EGS_HOME:-}" || -z "${HEN_HOUSE:-}" ]]; then
    echo "EGSnrc environment not set (EGS_HOME, HEN_HOUSE)." >&2
    exit 1
fi

inpDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/egsinp"
outDir="${1:?usage: runEdknrc.sh <outputDir> [pegsFile]}"
pegs="${2:-521icru}"
appDir="${EGS_HOME%/}/edknrc"

mkdir -p "$outDir"
outDir="$(cd "$outDir" && pwd)"

runOne() {
    local name="$1"
    local start end
    cp "$inpDir/$name.egsinp" "$appDir/"
    start=$(date +%s)
    (cd "$appDir" && nice -n 19 edknrc -i "$name" -p "$pegs" > "$outDir/$name.log" 2>&1)
    end=$(date +%s)
    for ext in egslst egsdat errors angplot radplot; do
        [[ -f "$appDir/$name.$ext" ]] && mv "$appDir/$name.$ext" "$outDir/"
    done
    rm -f "$appDir/$name.egsinp" "$appDir/$name.keV"
    echo "$name $((end - start))" >> "$outDir/walltime.txt"
}

for f in "$inpDir"/*.egsinp; do
    runOne "$(basename "$f" .egsinp)" &
done
wait
