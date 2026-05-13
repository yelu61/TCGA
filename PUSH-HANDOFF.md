# Push hand-off — public + private split

Everything is prepared. The public repo lives at `/tmp/tcga-toolkit-clean/`;
the private repo is the current project directory.

Before running anything, replace the `<YOUR-HANDLE>` placeholder in the
README / hand-off:

```bash
HANDLE=your-github-username

# Replace in public README
sed -i '' "s|<YOUR-HANDLE>|$HANDLE|g" /tmp/tcga-toolkit-clean/README.md

# Replace in private README
sed -i '' "s|<YOUR-HANDLE>|$HANDLE|g" \
  "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/TCGA/PRIVATE-README.md"
```

---

## A) Public repo — `<YOUR-HANDLE>/tcga-toolkit`

### 1. Create the empty GitHub repo

In the GitHub UI, create `tcga-toolkit` as **Public**. **Do not** initialise
with README / LICENSE / .gitignore — they already exist in the local copy.

### 2. Init & push

```bash
cd /tmp/tcga-toolkit-clean

git init -b main
git add .
git status                            # sanity-check what is staged
git commit -m "Initial public release of tcga-toolkit v0.3.0"

git remote add origin git@github.com:<YOUR-HANDLE>/tcga-toolkit.git
git push -u origin main
```

### 3. Verify

- README renders correctly on GitHub
- `bash tcga_toolkit/tests/run_all.sh` still passes locally
- No `0-Data/` / `GDCdata/` / `tcga_runs/` / personal templates leaked

---

## B) Private repo — `<YOUR-HANDLE>/tcga-research` (or your preferred name)

### 1. Create the empty GitHub repo

In the GitHub UI, create `tcga-research` as **Private**. Do not initialise.

### 2. Decide: submodule or sibling clone?

**Option A — git submodule (recommended)**. The toolkit appears at
`tcga_toolkit/` inside this repo but tracks the public repo's commits.

```bash
cd "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/TCGA"

# Move the working toolkit aside first (safety: backup, do not lose any
# local-only edits)
mv tcga_toolkit tcga_toolkit_backup_$(date +%Y%m%d)

git init -b main
git submodule add git@github.com:<YOUR-HANDLE>/tcga-toolkit.git tcga_toolkit

# Compare backup to fresh submodule; commit any local-only edits to the
# public repo first if you find diffs:
diff -r tcga_toolkit_backup_$(date +%Y%m%d) tcga_toolkit | head -50
rm -rf tcga_toolkit_backup_$(date +%Y%m%d)   # once you are sure
```

**Option B — sibling clone (simpler)**. Keep the toolkit working copy in
this repo, but `.gitignore` it. Edit `.gitignore` and uncomment the line:

```
# tcga_toolkit/
```

Then proceed without `git submodule add`.

### 3. First commit

```bash
cd "/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/TCGA"

git add .gitignore PRIVATE-README.md
git add -A AGENTS.md CLAUDE.md notebooks/ tcga_runs/ 2-Output/
git add -A .claude/settings.json   # local settings stays
git add Human\ DNA\ Repair\ Genes.xlsx   # opt-in — comment out if too sensitive
git status                          # sanity-check

git commit -m "Initial private research workspace"
git remote add origin git@github.com:<YOUR-HANDLE>/tcga-research.git
git push -u origin main
```

### 4. Verify

```bash
git status                          # should show clean working tree
du -sh .                            # the repo should be at most a few MB
                                    # if it shows GB, something leaked
git ls-files | xargs -I {} stat -f "%z %N" {} | sort -nr | head -20
                                    # top 20 largest tracked files
```

---

## C) Day-to-day workflow after the split

### When you change the toolkit code

```bash
cd tcga_toolkit             # inside the submodule
# edit / add a new task
bash tests/run_all.sh
git add . && git commit -m "Add task X" && git push   # pushes to public repo

cd ..
git add tcga_toolkit
git commit -m "Bump toolkit to include task X"
git push                    # records the new toolkit pointer in private repo
```

### When you change a research config

```bash
# directly in the private repo root
git add templates/ tcga_runs/
git commit -m "Run lasso-cox on FA pathway"
git push
```

---

## D) Data backup (handled outside git)

Your 21 GB of TCGA / GTEx / GDC data **never** goes into either repo.
Use:

- **Time Machine** for daily incremental backup.
- **External SSD / HDD** for quarterly snapshots.
- **Backblaze B2 / S3 Glacier** for yearly archives.

To recover on a new machine:

```bash
git clone --recurse-submodules git@github.com:<YOUR-HANDLE>/tcga-research.git TCGA
cd TCGA
# Restore 0-Data/, GDCdata/, GTEX/, 1-Input/ from your backup
bash tcga_toolkit/tests/run_all.sh
Rscript tcga_toolkit/scripts/run_task.R \
  --config tcga_toolkit/templates/audit_tcga.json
```
