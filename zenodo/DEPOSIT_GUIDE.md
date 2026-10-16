# Archiving this repository on Zenodo

Two routes are available. **Route A is strongly preferred** because Zenodo then re-archives the repository automatically every time you cut a GitHub release, and each release receives its own version DOI. This is the arrangement journals and reviewers expect.

Only one route is needed — but whichever you use, the metadata is already prepared:

- `.zenodo.json` (repository root) — read automatically by the GitHub integration
- `zenodo/metadata.json` — ready to paste into the web form, or to POST to the API

---

## Before you start: two placeholders to replace

1. In `CITATION.cff`, replace `REPLACE_WITH_GITHUB_USERNAME` (two occurrences) with the GitHub account or organisation that will host the repository.
2. Confirm the creator list and affiliations in both `.zenodo.json` and `zenodo/metadata.json`. ORCIDs for the second and third authors should be added if available.

---

## Route A — GitHub release + Zenodo integration (recommended)

### A1. Create the GitHub repository

```bash
cd hnscc-tls-signatures
git init
git add .
git commit -m "Initial release: multi-signature TLS analysis pipeline and results"

git branch -M main
git remote add origin https://github.com/<YOUR_USERNAME>/hnscc-tls-signatures.git
git push -u origin main
```

Create the empty repository on GitHub first (no README, no licence — this repository already contains both).

> **Repository visibility.** Choose **Public** if the manuscript's data availability statement promises open code. If you prefer to keep it private, Zenodo's GitHub integration does not support private repositories and you must use Route B instead.

### A2. Enable the Zenodo integration

1. Sign in to [zenodo.org](https://zenodo.org) using your GitHub account (or link GitHub under *Account → Linked accounts*).
2. Go to **Account → Settings → GitHub**.
3. Find `hnscc-tls-signatures` in the repository list and flip the switch **On**.
   Zenodo reads `.zenodo.json` at this point. If any field is rejected, the metadata appears in the deposit draft with that field blank — check the deposit after the first release.

### A3. Cut a release

1. On GitHub, go to **Releases → Draft a new release**.
2. Create the tag `v1.0.0` (matching the `version` field in the metadata).
3. Title: `v1.0.0 — initial release`.
4. Publish the release.

Within a minute or two, Zenodo archives the release as a new record and mints a DOI. The **badge** for the DOI appears in the GitHub repository sidebar under *Settings → Webhooks → Zenodo*, and on the Zenodo record page you will find both:

| DOI type | What it identifies | Use it for |
|---|---|---|
| **Version DOI** (e.g. `10.5281/zenodo.XXXXXXX`) | This exact release | Citing a specific analysis version |
| **Concept DOI** (e.g. `10.5281/zenodo.YYYYYYY`) | All versions of the software | The manuscript's data availability statement |

Use the **concept DOI** in the manuscript. It always resolves to the latest version.

### A4. Later versions

To archive a revision: bump `version` in `.zenodo.json` and `CITATION.cff`, commit, then draft a new release with a new tag (`v1.0.1`, `v1.1.0`, …). Zenodo files it under the same concept DOI automatically.

---

## Route B — manual Zenodo upload

Use this if the repository stays private, or if you prefer not to link accounts.

### B1. Build the archive

```bash
cd repo
zip -r hnscc-tls-signatures-v1.0.0.zip hnscc-tls-signatures \
    -x "*/__pycache__/*" "*/.git/*" "*.pyc"
```

### B2. Create the record

1. Sign in to Zenodo and click **New upload**.
2. **Upload type:** Software.
3. **DOI:** leave blank and click *Reserve DOI* — Zenodo allocates one and shows the format up front.
4. **Basic information:** copy each field from `zenodo/metadata.json`. Note that the `description` field accepts HTML; paste it as-is.
5. **Licence:** the metadata specifies MIT for the code. Because the repository also contains CC BY 4.0 data, state the split explicitly in the description (it is already worded that way) — Zenodo allows only one licence per record.
6. **Related works:** already listed in the metadata (`isDerivedFrom` the GDC TCGA-HNSC project; `isSupplementTo` the TCGA-HNSC marker paper).
7. **Funding:** add the two grants recorded in the `notes` field if you want them indexed:
   - Shapingba District Technology Innovation Project — 2025147, 2026095
   - Science-Health Collaborative Medical Research Project — Grant 55
8. **Contributors:** optional; add if any author is not in the creator list.

### B3. Publish

Click **Publish**. The DOI becomes active immediately and cannot be deleted afterwards (only a new version can be created).

---

## Route C — deposit via the REST API

Useful for automating the deposit from a script or CI job.

```bash
# 1. Create a personal access token with deposit:write and deposit:actions scopes
#    Account -> Settings -> Applications -> Personal access tokens
export ZENODO_TOKEN="<your token>"

# 2. Create an empty deposition
curl -s -H "Authorization: Bearer $ZENODO_TOKEN" \
     -H "Content-Type: application/json" \
     -d '{}' \
     https://zenodo.org/api/deposit/depositions | tee dep.json

# 3. Upload the archive (replace <DEPOSITION_ID> from dep.json)
curl -s -H "Authorization: Bearer $ZENODO_TOKEN" \
     -F "file=@hnscc-tls-signatures-v1.0.0.zip" \
     https://zenodo.org/api/deposit/depositions/<DEPOSITION_ID>/files

# 4. Attach the metadata prepared in zenodo/metadata.json
curl -s -X PUT -H "Authorization: Bearer $ZENODO_TOKEN" \
     -H "Content-Type: application/json" \
     -d @zenodo/metadata.json \
     https://zenodo.org/api/deposit/depositions/<DEPOSITION_ID>

# 5. Publish
curl -s -X POST -H "Authorization: Bearer $ZENODO_TOKEN" \
     https://zenodo.org/api/deposit/depositions/<DEPOSITION_ID>/actions/publish
```

`zenodo/metadata.json` is already wrapped in the `{"metadata": {...}}` envelope the API expects.

---

## After the DOI is minted

Three places need the DOI once it exists:

1. **`README.md`** — replace the sentence *"The archived release DOI will be added here once the Zenodo deposit is published."* with the concept DOI.
2. **`CITATION.cff`** — add a `doi:` field at the top level.
3. **The manuscript** — in *Availability of data and materials*, replace
   `〔待补：Zenodo/Figshare DOI〕` and `〔待补：GitHub 仓库链接〕` with the concept DOI and the repository URL.

Zenodo also generates a citation snippet and a DOI badge on the record page; the badge Markdown can be pasted at the top of `README.md`.

---

## Checklist

- [ ] `REPLACE_WITH_GITHUB_USERNAME` replaced in `CITATION.cff`
- [ ] Author ORCIDs completed (Xi Tang: 0000-0003-2740-1368; Yulian Zhang and Yungang He pending)
- [ ] Repository pushed and visible at the intended URL
- [ ] Zenodo integration switched on (Route A) *or* archive uploaded (Route B)
- [ ] Release tagged `v1.0.0` matching the `version` field
- [ ] DOI obtained — record both the version DOI and the concept DOI
- [ ] Concept DOI written into `README.md`, `CITATION.cff`, and the manuscript
- [ ] DOI badge added to `README.md`
