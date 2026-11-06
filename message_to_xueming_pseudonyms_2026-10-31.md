Subject: indelsig.net: replacing Hartwig sample IDs with HMF IDs

Hi Xueming,

Hartwig's publication policy says not to use Hartwig's internal sample IDs (CPCT..., DRUP..., WIDE..., ACTN..., CORE...) but to use the HMF... IDs in publications, and that would include indelsig.net.

Please make the Shiny repo XmWU123/indelsig (formerly indel-signature-browser) private for the time being.

I re-rendered the vignette pages and data tables with Hartwig's preferred IDs (of the form HMF003964A) in their place.

Please pull the latest `Liu2026_code_and_data` and copy these files into your local clone of XmWU123/indelsig:

- `vignette/connection_table.tsv` → `indel_shiny_clean/data/connection_table.tsv`
- `vignette/prot_table_1.csv` → `indel_shiny_clean/data/prot_table_1.csv`
- `vignette/separate_pages/*.html` (all 45 signature pages and `overview_table.html`) → `indel_shiny_clean/www/separate_pages/`

```
# In your local clone of XmWU123/indelsig, after copying in the new files:
cd <your local clone of XmWU123/indelsig>

# Remove the combined vignette page. It has the old IDs, nothing on the
# site links to it, and the new version (109 MB) is over GitHub's 100 MB
# file limit.
git rm indel_shiny_clean/www/repertoire_of_indel_mutational_signatures.html

git add -A
git commit -m "Replace Hartwig sample IDs with HMF IDs"

# Look for old IDs in the most recent commit. If you find any, stop and
# let me know. (The digits in the pattern matter: without them it also
# matches random text inside the images embedded in the HTML pages.)
git grep -n -E '(CPCT|DRUP|WIDE|ACTN|CORE)[0-9]{6}' HEAD

# The repository's git history still contains the Hartwig IDs, and
# replacing the files does not remove them from older commits.
# The next steps make a copy of the current files without the commit
# history. Move to the folder that contains your clone:
cd ..
git -C <your local clone of XmWU123/indelsig> archive --prefix=indelsig-browser/ HEAD | tar -x
# The command above writes the files from the latest commit (HEAD) of
# your clone into a new folder, indelsig-browser/, in the current
# directory. Only tracked files are copied. The .git folder and its
# history, uncommitted edits, and untracked or ignored files are all left
# out. git archive packs the files of the commit, and | tar -x unpacks
# them straight away, without writing a tar file to disk.

cd indelsig-browser
git init -b main
git add -A
git commit -m "Initial commit"

# In indelsig-browser, double check again that no Hartwig IDs are left:
rg -n '(CPCT|DRUP|WIDE|ACTN|CORE)[0-9]{6}'
# rg is ripgrep, a fast and flexible successor to grep.
# If you find old IDs in indelsig-browser, let me know where they are
# and wait.

# Create the new empty repo (XmWU123/indelsig-browser) on GitHub first,
# then:
git remote add origin git@github.com:XmWU123/indelsig-browser.git
git push -u origin main
# Make sure indelsig-browser is public.
```

Then update the web site (indelsig.net) to be populated from the new repo.

Then, on Zenodo, make a new version of the existing record (concept DOI 10.5281/zenodo.20092288, which the paper cites), rather than a new record, with a .zip of the new indelsig-browser repo. In the new version, set the Identifiers URL to https://indelsig.net and the repository URLs to https://github.com/XmWU123/indelsig-browser. Please do not turn on Zenodo's GitHub integration for indelsig-browser, because its first release would get a new concept DOI. The files of version v1.0.0 (10.5281/zenodo.20092289) still have the old IDs, and only Zenodo support can remove them. Could you ask them to (https://zenodo.org/support)?

Keep the old repo XmWU123/indelsig private for now. We can decide later whether to delete it.

Please let me know if you have any questions.

Thanks,
Steve
