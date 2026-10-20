/*
 * PolyGEN browser demo — worker.
 *
 * Loads Pyodide, then runs the unmodified PolyGEN engine from this repository.
 * Everything happens off the main thread so a long optimisation can never
 * freeze the page; the main thread can terminate this worker to cancel.
 */

let pyodide = null;

function progress(text){ self.postMessage({type:"progress", text}); }

async function fetchText(url){
  const r = await fetch(url, {cache:"force-cache"});
  if(!r.ok) throw new Error(`HTTP ${r.status} for ${url}`);
  return await r.text();
}

async function boot({pyodideUrl, pkgBase}){
  progress("Downloading Python runtime (about 10 MB, cached afterwards)…");
  const {loadPyodide} = await import(`${pyodideUrl}pyodide.mjs`);
  pyodide = await loadPyodide({indexURL: pyodideUrl});

  progress("Loading numpy, pandas and biopython…");
  await pyodide.loadPackage(["numpy", "pandas", "biopython"]);

  progress("Loading the PolyGEN engine…");
  const FS = pyodide.FS;
  FS.mkdirTree("/pg/polygen_engine/overhangsets");

  // engine.py imports SeqIO, SeqRecord and Bio.Alphabet.IUPAC but never uses
  // them. IUPAC was removed in Biopython 1.78, so the imports are stripped here
  // to keep this page working against an older checkout that still pins 1.77.
  const mods = ["__init__.py", "cli.py", "engine.py"];
  const texts = await Promise.all(mods.map(m => fetchText(pkgBase + m)));
  mods.forEach((m, i) => {
    let src = texts[i];
    if(m === "engine.py"){
      src = src
        .replace(/^from Bio import SeqIO[ \t]*$/m, "")
        .replace(/^from Bio\.SeqRecord import SeqRecord[ \t]*$/m, "")
        .replace(/^from Bio\.Alphabet import IUPAC[ \t]*$/m, "");
    }
    FS.writeFile(`/pg/polygen_engine/${m}`, src);
  });

  // all_curated_overhangs() lists this directory, so every table must be here.
  const sets = [];
  for(let n = 10; n <= 50; n++) sets.push(`setsof${n}.csv`);
  const csvs = await Promise.all(sets.map(f => fetchText(`${pkgBase}overhangsets/${f}`)));
  sets.forEach((f, i) => FS.writeFile(`/pg/polygen_engine/overhangsets/${f}`, csvs[i]));

  const version = pyodide.runPython(`
import sys
sys.path.insert(0, "/pg")
import polygen_engine
polygen_engine.__version__
`);
  self.postMessage({type:"ready", version});
}

function design(cfg){
  pyodide.globals.set("cfg_json", JSON.stringify(cfg));
  return JSON.parse(pyodide.runPython(`
import json
from polygen_engine import PTGbldr, scarless_gg, annotatePrimers, polyToJson, InvalidUsage

cfg = json.loads(cfg_json)
try:
    ptg = PTGbldr(cfg["name"], [list(p) for p in cfg["parts"]], cfg["poltype"])
    poly = scarless_gg(ptg,
                       poltype=cfg["poltype"], enzm=cfg["enzm"],
                       tm_range=cfg["tm_range"],
                       bb_linkers=cfg["bb_linkers"], ad_linkers=[])
    poly = annotatePrimers(poly,
                           oligo_prefix=cfg["prefix"], oligo_index=cfg["index"],
                           staticBorderPrimers=False, noBorderPrimers=False,
                           poltype=cfg["poltype"], enzm=cfg["enzm"],
                           bb_linkers=cfg["bb_linkers"], ad_linkers=[])
    out = {"ok": True, "data": polyToJson(poly)}
except InvalidUsage as e:
    out = {"ok": False, "error": str(getattr(e, "message", e))}
except MemoryError:
    out = {"ok": False, "error": "Ran out of memory while searching for a compatible overhang set. Try fewer sequences, or use the full tool."}
except Exception as e:
    out = {"ok": False, "error": f"{type(e).__name__}: {e}"}

json.dumps(out)
`));
}

self.onmessage = async (ev) => {
  const msg = ev.data;
  try{
    if(msg.type === "boot"){
      await boot(msg);
    }else if(msg.type === "design"){
      const res = design(msg.cfg);
      self.postMessage({type:"result", id: msg.id, ...res});
    }
  }catch(err){
    if(msg.type === "boot") self.postMessage({type:"bootError", error: err.message});
    else self.postMessage({type:"result", id: msg.id, ok:false, error: err.message});
  }
};
