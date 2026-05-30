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
- TasH targets are represented as two 8 bp programmable spacer regions.
- The mature tigRNA is assembled as edge repeat 5 prime, spacer A, loop repeat, spacer B and edge repeat 3 prime.
- Exact target windows are 18 bp for TasA and 16 bp for TasH.
- Exact target windows can be separated with new lines, spaces, commas, semicolons or vertical bars.
- Scaffold fields can be edited for locus-specific tigRNA scaffolds. DNA or RNA bases are accepted, but the final mature tigRNA must be 36 nt.

Designed tigRNAs can be selected and sent directly to the PTG page using the `tigRNA` architecture. Only checked rows are transferred; if no rows are checked, the full generated set is transferred. In this architecture, sequences are provided as mature tigRNAs:

```text
tigRNA;sequence0|tigRNA;sequence1|tigRNA;sequence2
```

The `tigRNA` architecture does not add tRNAs and does not add a Cpf1 direct repeat. It treats each input as a mature tigRNA unit, then designs oligo-extension fragments for scarless Golden Gate assembly of multiplex tigRNA arrays.

For tigRNA oligo design, PolyGEN chooses internal Golden Gate overhangs from spacer B only, so fixed edge and loop repeat regions are not used as variable assembly overhang sources. Each fragment is produced from a forward/reverse oligo pair with an overlapping region that can be filled in by polymerase. The reported oligo melting temperature for tigRNA fragments is calculated only from that shared overlap, not from the full oligo sequence.

The `Reuse PTG border oligos` option is applied only to PTG/Cas9 tRNA-processed assemblies. These reusable border oligos are specific to the selected restriction enzyme and 4 bp border linkers. CA and `tigRNA` assemblies do not use the tRNA architecture, so PolyGEN always designs new target-specific border oligos for those modes.

The GenBank output annotates each mature tigRNA with edge repeat 5 prime, spacer A, loop repeat, spacer B and edge repeat 3 prime features.

_____________

**Setup Linux/Mac**

First, install docker by running apt-get install docker

In the terminal, run through the following pipeline

- Clone this repository via `git clone https://git.hhu.de/urquizag/polygen`
- Navigate into the cloned repo `cd polygen`
- execute `docker-compose up` (requires docker desktop to be active)

In the browser, open localhost:5000

## Cloud Run / container deployment

The Docker image runs PolyGEN with Gunicorn and listens on the `PORT`
environment variable, which is required by Cloud Run and similar container
hosts. Locally, Docker maps container port 8080 to host port 5000:

```bash
docker build -t polygen:test .
docker run --rm -p 5000:8080 -e PORT=8080 polygen:test
```

For a Cloud Run test deployment from the repository root:

```bash
gcloud run deploy polygen-test --source . --allow-unauthenticated --max-instances=1
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
