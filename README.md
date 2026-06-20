# **PolyGEN**

Design automation of polycistronic tRNA-based genes containing custom RNAs for assembly in Type IIs restriction enzyme-driven Golden Gate experiments. The backbone of PolyGEN is based on
[iBioCAD](https://ibiocad.igb.illinois.edu/) by HamediRad et al. (2019), for which the code is openly available [here](https://github.com/scottweis1/iBioCAD).

The code takes as input an array of custom RNAs and will compute the finished PTG together with the necessary oligomers to produce all parts from a plasmid containing a gRNA-tRNA template. Currently, the produced PTGs can include sgRNAs, pegRNAs, crRNAs, tigRNAs and other small RNAs. By default, PolyGEN will use the following parameters, which can be varied

- primer melting temperature between 55 and 65 °C if possible
- Digestion by BsaI
- 'tgcc' and 'gttt' as restriction overlaps with the plasmid
- no additional restriction sites in the plasmid

To calculate the primer melting temperatures, PolyGEN uses the same method and parameters as Benchling: [SantaLucia (1998)](https://www.pnas.org/content/95/4/1460).

## Prime editing PBS design

For pegRNA design, PolyGEN tests PBS lengths from 8 to 17 nt and chooses the PBS with melting temperature closest to 30 °C. The GenBank output annotates the selected PBS length instead of assuming a fixed 13 nt PBS.

## TasA/TasH tigRNA design

PolyGEN includes a Tas guide design page for TIGR-Tas guide RNAs. It currently supports TasA and TasH mature tigRNA architectures.

- TasA targets are represented as two 9 bp programmable spacer regions.
- TasH multiplexing targets are represented as two 8 bp programmable spacer regions. The first 8 bp become spacer A and the reverse complement of the second 8 bp becomes spacer B.
- The mature tigRNA is assembled as edge repeat 5 prime, spacer A, loop repeat, spacer B and edge repeat 3 prime.
- The default TasA scaffold is AACCG, spacer A, AGTAACCCC, spacer B, AGTG. Adjacent TasA array units repeat as right edge followed by the next left edge.
- For adjacent asymmetric Tas units, the upstream right edge plus the downstream left edge is treated as required edge-junction context. For default TasA this is AGTG-AACCG; for default TasH this is AAGCGAGCCA-GAGCGAGTTA.
- Exact target windows are 18 bp for TasA and 16 bp for TasH.
- Exact target windows can be separated with new lines, spaces, commas, semicolons or vertical bars.
- Scaffold fields can be edited for locus-specific tigRNA scaffolds. DNA or RNA bases are accepted, but the final mature tigRNA must be 36 nt for TasA or the multiplexing unit must be 46 nt for TasH. The TasH default is the asymmetric unit `GAGCGAGTTA` + spacer A + `AAAACAATCA` + spacer B + `AAGCGAGCCA`.

Designed tigRNAs can be selected and sent directly to the PTG page using the `tigRNA` architecture. Only checked rows are transferred; if no rows are checked, the full generated set is transferred. In this architecture, sequences are provided as mature tigRNAs:

```text
tigRNA;sequence0|tigRNA;sequence1|tigRNA;sequence2
```

The `tigRNA` architecture does not add tRNAs and does not add a Cpf1 direct repeat. It treats each input as a mature tigRNA unit, then designs oligo-extension fragments for scarless Golden Gate assembly of multiplex tigRNA arrays.

For tigRNA oligo design, PolyGEN chooses internal Golden Gate overhangs from spacer A or spacer B, so fixed edge and loop repeat regions are not used as variable assembly overhang sources. Spacer-A splits are accepted only when enough left-side sequence remains for the oligo-extension fragment. Each fragment is produced from a forward/reverse oligo pair with an overlapping region that can be filled in by polymerase. PolyGEN dynamically balances tigRNA oligo arms up to 60 bp when required, so short TasA/TasH fragments keep enough shared overlap for assembly. The reported oligo melting temperature for tigRNA fragments is calculated only from that shared overlap, not from the full oligo sequence.

TasA and TasH overhang selection is Tm-aware. PolyGEN scores valid spacer-A and spacer-B splits by the predicted fill-in overlap Tm of the final oligo pair and prefers combinations where every tigRNA fragment is at least 45 °C. If no valid combination reaches that floor, the best available design is kept and the warning is written to the oligo page, CSV, raw JSON and GenBank outputs.

If no complete optimal overhang set can be found, PolyGEN enters rescue overhang mode. In rescue mode, it still chooses overhangs that appear in curated PolyGEN overhang tables, but it flags that the exact selected collection was not found as a validated optimal set. The warning is shown on the oligo page and written into the CSV, raw JSON and GenBank outputs.

The `Reuse PTG border oligos` option is applied only to PTG/Cas9 tRNA-processed assemblies. These reusable border oligos are specific to the selected restriction enzyme and 4 bp border linkers. CA and `tigRNA` assemblies do not use the tRNA architecture, so PolyGEN always designs new target-specific border oligos for those modes.

The GenBank output annotates each mature tigRNA with edge repeat 5 prime, spacer A, loop repeat, spacer B and edge repeat 3 prime features. Oligo primer-binding annotations show only the oligo sequence that matches the final assembled construct: the Type IIS recognition site is excluded and the 4 bp Golden Gate overhang is included.

## GenBank annotations

PolyGEN annotates generated GenBank files with readable labels and notes for each array unit, spacer, scaffold, tRNA, direct repeat, tigRNA subfeature, pegRNA PBS/RT-template segment and primer-binding site. These labels preserve the link between each RNA part and the sequence it targets after the array has been assembled.

For TasH multiplexing units, GenBank exports annotate the left edge repeat, spacer A, loop repeat, spacer B and right edge repeat. If an older identical-edge TasH scaffold is supplied, PolyGEN can also mark the processing cuts inside those identical edge repeats.

For TasH arrays with identical left and right edge repeats, adjacent units can share the identical edge repeat at each junction. The current asymmetric TasH default is assembled as complete repeated units, while its AAGCGAGCCA-GAGCGAGTTA inter-unit edge junction is preserved in the oligo overlap.

TasH oligo design follows the selected scaffold model: asymmetric defaults are split across complete units with the required right-edge/left-edge junction retained, while identical-edge custom scaffolds can still use the shared-edge path.

_____________

**Setup Linux/Mac**

First, install docker by running apt-get install docker

In the terminal, run through the following pipeline

- Clone this repository via `git clone https://git.hhu.de/urquizag/polygen`
- Navigate into the cloned repo `cd polygen`
- execute `docker-compose up` (requires docker desktop to be active)

In the browser, open localhost:5000

## PyPI engine package

The reusable Python engine is published separately as `polygen-engine`:

```bash
pip install polygen-engine
```

For local development from this repository:

```bash
python3.9 -m venv .venv
. .venv/bin/activate
pip install -e .
polygen-engine --version
```

The engine package exposes the design API, including `PTGbldr`, `scarless_gg`,
`tas_guide_design`, `TAS_SYSTEMS` and the TasA/TasH tigRNA helpers. The Flask
web app remains in `polygen_scripts` and is deployed with Docker/Gunicorn.

## Cloud Run / container deployment

The Docker image runs PolyGEN with Gunicorn and listens on the `PORT`
environment variable, which is required by Cloud Run and similar container
hosts. Locally, Docker maps container port 8080 to host port 5000:

```bash
docker build -t polygen:test .
docker run --rm -p 5000:8080 -e PORT=8080 polygen:test
```

For a constrained Cloud Run demo deployment, run from the repository root and
verify that `Dockerfile` is listed before deploying. If Cloud Run says it is
"Building using Buildpacks", rerun the command from the directory containing
`Dockerfile`.

```bash
ls Dockerfile
gcloud run deploy polygen-v1-demo --source . --region=europe-west3 --allow-unauthenticated --min-instances=0 --max-instances=1 --concurrency=10 --cpu=1 --memory=512Mi --timeout=120 --set-env-vars SECRET_KEY=<new_secret>
```

For a public test service, set a real `SECRET_KEY` environment variable in the
Cloud Run service settings. The app keeps generated results in per-browser
server-side session state, so run a single Gunicorn worker and keep Cloud Run
to one instance for the first demo if downloads must survive between requests.

**Setup Windows**

Activate Windows Subsystem for Linux (WSL2) by 

- open Control Panel
- open Turn Windows features on or off
- check features Virtual Machine Platform and Windows Subsystem for Linux
- confirm with OK
- restart the pc
- download the Update Setup from [here](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi)
- Execute the wsl_update_x64.msi file
- open a command prompt
- run wsl --set-default-version 2

In the BIOS, enable Visualization tools. Next, install [Docker](https://docs.docker.com/desktop/windows/install/) and clone this git repository. With active docker, open a command prompt and navigate into the cloned repository using dir \<location\>. Execute docker-compose up.

In the browser, open localhost:5000
  
____________
  
**Common problems**
  
If the install fails due to issues with docker-snap:
  
- sudo rm -rf /etc/docker
- sudo snap refresh
