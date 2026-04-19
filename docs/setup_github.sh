#!/bin/bash
# =============================================================================
# GitHub Setup Script — Parkinson RNA-Seq Analysis
# Run this script once from your Ubuntu terminal
# Usage: bash docs/setup_github.sh
# =============================================================================

set -e  # stop on any error

echo "============================================="
echo "  Parkinson RNA-Seq — GitHub Setup Script"
echo "============================================="

# ── Step 1: Check git is installed ───────────────────────────────────────────
if ! command -v git &> /dev/null; then
  echo "Installing git..."
  sudo apt update && sudo apt install git -y
fi
echo "✅ git is ready"

# ── Step 2: Configure git (edit these) ───────────────────────────────────────
echo ""
echo "Setting up git config..."
read -p "Enter your GitHub username: " GH_USER
read -p "Enter your GitHub email: " GH_EMAIL

git config --global user.name  "$GH_USER"
git config --global user.email "$GH_EMAIL"
echo "✅ git config set"

# ── Step 3: Init repo ─────────────────────────────────────────────────────────
cd "$(dirname "$0")/.."   # go to repo root
git init
echo "✅ git init done"

# ── Step 4: Create placeholder files so empty dirs are tracked ───────────────
touch data/.gitkeep
touch results/plots/.gitkeep
echo "✅ placeholder files created"

# ── Step 5: Initial commit ────────────────────────────────────────────────────
git add .
git commit -m "feat: initial project structure and analysis pipeline

- Complete R analysis script (EDA, PCA 2D/3D, heatmaps, variance)
- Professional README with results summary and usage instructions
- .gitignore excluding large data and generated plots
- MIT License"

echo "✅ initial commit done"

# ── Step 6: Instructions for GitHub remote ───────────────────────────────────
echo ""
echo "============================================="
echo "  NEXT STEPS — Do these manually:"
echo "============================================="
echo ""
echo "1. Go to https://github.com/new"
echo "   Repository name : parkinson-rnaseq-analysis"
echo "   Visibility       : Public"
echo "   DON'T add README (we already have one)"
echo ""
echo "2. Then run these commands:"
echo ""
echo "   git remote add origin https://github.com/$GH_USER/parkinson-rnaseq-analysis.git"
echo "   git branch -M main"
echo "   git push -u origin main"
echo ""
echo "3. To push your data files (optional — they're in .gitignore by default):"
echo "   git add -f data/Parkinson_exp.txt data/Parkinson_phenotable.txt"
echo "   git commit -m 'data: add expression matrix and phenotype table'"
echo "   git push"
echo ""
echo "✅ All done! Your repo is ready 🎉"
