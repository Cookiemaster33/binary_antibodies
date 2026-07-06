#!/bin/bash
# ============================================================
# push_results.sh
# Copy pipeline outputs into pipeline_results/<folder> and push to GitHub.
#
# Environment:
#   RESULTS_GITHUB_DIR  — target subfolder (default: pipeline_results/v5_two_round_refine)
#   GITHUB_BRANCH       — branch to push (default: cursor/conditional-nanobody-design-992c)
#   GITHUB_TOKEN        — required on the Lambda instance for git push
#   GITHUB_REPO         — default: Cookiemaster33/binary_antibodies
# ============================================================
set -eo pipefail

PIPELINE=${PIPELINE_DIR:-/home/ubuntu/pipeline}
RESULTS_SUBDIR=${RESULTS_GITHUB_DIR:-pipeline_results/v5_two_round_refine}
BRANCH=${GITHUB_BRANCH:-cursor/conditional-nanobody-design-992c}
REPO_SLUG=${GITHUB_REPO:-Cookiemaster33/binary_antibodies}
REPO_DIR=${GITHUB_REPO_DIR:-/home/ubuntu/binary_antibodies}

if [ -z "${GITHUB_TOKEN:-}" ]; then
    echo "WARNING: GITHUB_TOKEN not set — skipping GitHub push."
    exit 0
fi

DEST="$REPO_DIR/$RESULTS_SUBDIR"
mkdir -p "$DEST/structures"

echo "Collecting results → $RESULTS_SUBDIR"

cp "$PIPELINE/final/final_results.json" "$DEST/" 2>/dev/null || true
cp "$PIPELINE/final/round1_final_results.json" "$DEST/" 2>/dev/null || true
cp "$PIPELINE/final/validated_minibinders.fasta" "$DEST/" 2>/dev/null || true
cp "$PIPELINE/outputs/mpnn/minibinder_sequences.fasta" "$DEST/mpnn_minibinder_sequences.fasta" 2>/dev/null || true
cp "$PIPELINE/top_structures/"*.tsv "$DEST/structures/" 2>/dev/null || true
cp "$PIPELINE/top_structures/"*.cif "$DEST/structures/" 2>/dev/null || true
cp "$PIPELINE/full_pipeline.log" "$DEST/pipeline.log" 2>/dev/null || true

# Run metadata for reproducibility
cat > "$DEST/run_config.json" << EOF
{
  "results_dir": "$RESULTS_SUBDIR",
  "rfd3_rounds": "${RFD3_ROUNDS:-1}",
  "n_designs_round1": "${N_DESIGNS:-200}",
  "mb_length_range": "${MB_LENGTH_RANGE:-35-70}",
  "rfd3_round2_templates": "${RFD3_ROUND2_TEMPLATES:-5}",
  "rfd3_round2_designs_per_template": "${RFD3_ROUND2_DESIGNS_PER_TEMPLATE:-8}",
  "rfd3_partial_t": "${RFD3_PARTIAL_T:-2.0}",
  "completed_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF

# Clone or update repo
if [ ! -d "$REPO_DIR/.git" ]; then
    git clone "https://x-access-token:${GITHUB_TOKEN}@github.com/${REPO_SLUG}.git" "$REPO_DIR"
fi

cd "$REPO_DIR"
git config user.email "cursor-agent@users.noreply.github.com"
git config user.name "Cursor Agent"
git fetch origin "$BRANCH" 2>/dev/null || true
git checkout "$BRANCH" 2>/dev/null || git checkout -b "$BRANCH" "origin/$BRANCH" 2>/dev/null || git checkout -b "$BRANCH"
git pull --rebase origin "$BRANCH" 2>/dev/null || true

git add "$RESULTS_SUBDIR"
if git diff --cached --quiet; then
    echo "No new results to commit."
    exit 0
fi

git commit -m "results: ${RESULTS_SUBDIR} — RFD3_ROUNDS=${RFD3_ROUNDS:-1} pipeline run"
git push "https://x-access-token:${GITHUB_TOKEN}@github.com/${REPO_SLUG}.git" "$BRANCH"
echo "Pushed results to $BRANCH:$RESULTS_SUBDIR"
