# Before the course: install R and RStudio

**Please do this at least two days before the course**, not the night before.
It takes about 20–30 minutes, most of it waiting.

If something goes wrong, email the organizer with the console output. Do not
turn up on the day with a broken setup — we will not have time to fix it.

---

## Step 1 — Install R

R is the language. RStudio (step 2) is the program you actually look at. You
need both, and R must go first.

Go to **<https://cloud.r-project.org>** and download R for your system.

- **Mac**: choose the `.pkg` file. There are two — pick **arm64** if your Mac
  is M1/M2/M3/M4, **x86_64** if it is an older Intel Mac. (Apple menu → About
  This Mac will tell you which chip you have.)
- **Windows**: click "Download R for Windows" → "base" → "Download R for
  Windows".

Open the downloaded file and follow the installer. Accept all the defaults.

> **Already have R?** Open RStudio and look at the first line in the Console.
> If it says anything below **4.4.0**, please install the current version
> anyway. Older versions have to build packages from scratch, which is where
> most installation failures happen.

### Windows only: also install Rtools

Windows needs one extra piece. Get it from
<https://cran.r-project.org/bin/windows/Rtools/> — pick the version matching
your R (Rtools45 for R 4.5.x, Rtools44 for R 4.4.x). Install with defaults.

---

## Step 2 — Install RStudio

Go to **<https://posit.co/download/rstudio-desktop/>** and download RStudio
Desktop (free version). Install it.

---

## Step 3 — Get the course files

Go to **<https://github.com/savannahmwesigwa/bulk-rnaseq-course>** and click
the green **Code** button → **Download ZIP**.

Unzip it somewhere you will find again — your Desktop is fine. You should
end up with a folder containing `00_setup.R` and a `data` folder.

---

## Step 4 — Run the setup script

1. Open **RStudio**
2. **File → Open File…** and choose `00_setup.R` from the folder you just
   unzipped
3. Click **Source** (top right of the editor pane), or press
   **Ctrl+Shift+S** (**Cmd+Shift+S** on Mac)
4. Wait. It will print a lot of text. The first run takes 10–20 minutes.
5. If it asks *"Do you want to install from sources the package which needs
   compilation?"* — answer **no** (type `n` and press Enter).

### What you should see at the end

```
========================================
  ALL GOOD - you are ready for the course.
========================================
```

If you see that, you are done. Close RStudio and come to the course.

### If you see FAIL

```
========================================
  SOMETHING DID NOT WORK.
========================================
```

Select everything in the Console (click in it, then Ctrl+A / Cmd+A), copy,
and email it to **gushijima@mdanderson.org**. Include what kind of computer
you have. We will sort it out before the day.

---

## Common problems

**"Installation of package 'X' had non-zero exit status"**
Usually means it tried to build from source. Run the setup script again —
it often works second time. If not, send us the output.

**Mac: "cannot be opened because the developer cannot be verified"**
Right-click the installer → Open → Open. Or System Settings → Privacy &
Security → "Open Anyway".

**Mac: errors mentioning `gfortran` or `clang`**
You need Xcode command line tools. Open Terminal and run:
```
xcode-select --install
```
Then run the setup script again.

**It has been running for 40 minutes**
That can be normal on a slow connection — the Bioconductor packages are
large. If nothing has printed for over 15 minutes, restart R
(Session → Restart R) and run it again.

**"package 'DESeq2' is not available"**
Your R is too old. Go back to step 1.

---

## What we will cover

- What a count matrix is and why raw counts are not comparable
- CPM, RPKM, TPM — what they do and which to use when
- Checking whether your samples are what you think they are
- Finding contamination from co-cultured cells
- Running a differential expression analysis
- Reading the results and knowing when not to trust them

You do not need any prior R experience. Every line of code is explained.
