# System of Linear Equations Calculator (web version, powered by the MATLAB code)

A web page for solving 3 to 10 linear equations with Doolittle, Crout, Cholesky, Cramer's Rule,
Gauss-Jordan, Gauss-Seidel and Gauss-Jacobi. **The browser does no math.** When you press
*Solve with the MATLAB code*, the page sends your system to the server, which runs the project's
MATLAB solver functions in GNU Octave and sends back the Command Window printout.

## What is in this folder

| File | What it does |
|---|---|
| `matlab/solver_core.m` | The solver functions, identical to the ones in `LinearSystemSolver.m` and `LinearSystemSolverGUI.m` |
| `matlab/web_solve.m` | Reads one request, calls the right solver function, saves the printout and the answer |
| `server.py` | Small Python server: shows the page and runs `web_solve.m` in Octave for each request |
| `static/index.html` | The web page (input form + display of the MATLAB printout) |
| `Dockerfile` | Builds a Linux image with GNU Octave and Python |
| `render.yaml` | Tells Render how to deploy it |

## Put it online with Render (free)

You need a GitHub account and a Render account (you can sign up to Render with GitHub).

1. **Create a GitHub repository.** On github.com click **New**, name it
   `linear-equations-calculator`, make it Public, and click **Create repository**.
2. **Upload the files.** On the new repository page click **uploading an existing file**, then drag
   in everything from this folder: `Dockerfile`, `render.yaml`, `server.py`, `README.md`, and the
   `matlab` and `static` folders. Click **Commit changes**.
3. **Create the service on Render.** At https://dashboard.render.com choose **New > Blueprint**,
   connect your GitHub account if asked, pick the `linear-equations-calculator` repository, and
   click **Apply**. Render reads `render.yaml` and creates a free Docker web service.

   (Alternative without the Blueprint: **New > Web Service**, pick the repository, set
   *Language* to **Docker** and *Instance Type* to **Free**, then click **Deploy Web Service**.)
4. **Wait for the first build.** Installing Octave takes several minutes the first time. When the
   log shows `Serving on port ... using GNU Octave`, the service is live.
5. **Open your link.** It looks like `https://linear-equations-calculator.onrender.com`. Share that
   link with anyone; no account is needed to use it.

Every time you upload a change to the repository, Render rebuilds and redeploys automatically.

**About the free plan:** the service goes to sleep after 15 minutes without visitors. The next
visit wakes it up, which takes about a minute; the page shows a "waking up" message meanwhile.
After that, each solve takes about a second.

## Lecture mode (for teaching)

After solving a system, press **Lecture mode** in the Result card. The solution then plays like
slides, one calculation per step, in large type for a projector:

- Each step shows the matrices with the value **being computed** in solid blue and the values
  **used** in that step outlined, plus the formula with the numbers substituted.
- **Pause before answers** (eye button, or **R**) hides each result behind a **?** so the class can
  work it out first. Press **Next** once to show the answer, and again to move on.
- **Teaching notes** (light-bulb button, or **N**) show a one-line explanation of each step.
- The **step list** on the left (or **O**) jumps to any step.
- Keys: **→ / Space / Page Down** next, **← / Page Up** back, **Home / End**, **+ / −** text size,
  **F** full screen, **Esc** close. Presentation clickers that send Page Down / Page Up also work.

Lecture mode does no calculations of its own: every number on its slides is copied from the
printout of the MATLAB solver code.

## Run it on your own computer

- **With Docker Desktop:**
  `docker build -t lss .` then `docker run -p 7860:7860 lss`, and open http://localhost:7860
- **Without Docker (Linux, macOS or WSL):** install GNU Octave and Python 3
  (Ubuntu: `sudo apt install octave python3`), then run `python3 server.py` in this folder and
  open http://localhost:7860.
  On Windows, install Octave from https://octave.org and make sure `octave-cli` is on your PATH.

## Server settings (optional)

| Environment variable | Default | Meaning |
|---|---|---|
| `PORT` | 7860 | Port to listen on (Render sets this for you) |
| `SOLVE_TIMEOUT` | 60 | Seconds before a solve is stopped |
| `MAX_OCTAVE` | 2 | How many Octave processes may run at the same time |
