#!/usr/bin/env python3
"""
ensembl.py

Flask backend for Ensembl annotation with UniProt & NCBI fallback support.

Serves:
  - GET  /health
  - GET  /version
  - GET  /config
  - POST /annotate
  - GET  /annotate
  - GET  /           -> serves index.html from project root
  - GET  /index.html -> serves index.html
"""
from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import requests, time, sys, os, json
from urllib.parse import quote
from datetime import datetime

# ----- config -----
ENSEMBL_REST = "https://rest.ensembl.org"
UNIPROT_REST = "https://rest.uniprot.org"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS = {"Accept": "application/json"}
SLEEP_BETWEEN = 0.08  # polite delay between Ensembl calls
SLEEP_BETWEEN_UNIPROT = 0.1  # polite delay between UniProt calls
SLEEP_BETWEEN_NCBI = 0.1  # polite delay between NCBI calls
MAX_RETRIES = 5
MAX_IDS = 200  # safety limit
VERSION = "v1.0.0"
ENABLE_UNIPROT_FALLBACK = True  # Enable/disable UniProt fallback
ENABLE_NCBI_FALLBACK = True  # Enable/disable NCBI fallback

app = Flask(__name__, static_folder=None)
CORS(app)  # allow cross-origin (safe for local dev)

# ----- HTTP with retry/backoff -----
def retry_get(url, params=None, headers=HEADERS, max_tries=MAX_RETRIES):
    backoff = 1.0
    for attempt in range(max_tries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=30)
            if r.status_code == 200:
                return r
            # rate-limited or service unavailable -> backoff and retry
            if r.status_code in (429, 503):
                time.sleep(backoff)
                backoff *= 2
                continue
            # other non-200 -> return response for caller to inspect
            return r
        except requests.RequestException:
            time.sleep(backoff)
            backoff *= 2
    raise RuntimeError(f"Failed to GET {url} after {max_tries} attempts")

def retry_post(url, data=None, headers=HEADERS, max_tries=MAX_RETRIES):
    backoff = 1.0
    for attempt in range(max_tries):
        try:
            r = requests.post(url, data=data, headers=headers, timeout=30)
            if r.status_code == 200:
                return r
            # rate-limited or service unavailable -> backoff and retry
            if r.status_code in (429, 503):
                time.sleep(backoff)
                backoff *= 2
                continue
            # other non-200 -> return response for caller to inspect
            return r
        except requests.RequestException:
            time.sleep(backoff)
            backoff *= 2
    raise RuntimeError(f"Failed to POST {url} after {max_tries} attempts")

# ----- Ensembl helpers -----
def get_gene_symbol(ensembl_id):
    """Return the display_name (gene symbol) for a gene Ensembl ID or None."""
    url = f"{ENSEMBL_REST}/lookup/id/{quote(ensembl_id)}"
    r = retry_get(url)
    if r is None or r.status_code != 200:
        return None
    try:
        j = r.json()
    except Exception:
        return None
    return j.get("display_name") or j.get("external_name") or None

def get_go_xrefs(ensembl_id):
    """Return list of (go_id, description) tuples from Ensembl xrefs for a gene."""
    url = f"{ENSEMBL_REST}/xrefs/id/{quote(ensembl_id)}"
    r = retry_get(url)
    if r is None or r.status_code != 200:
        return []
    try:
        items = r.json()
    except Exception:
        return []
    gos = []
    for it in items:
        dbname = (it.get("dbname") or "").upper()
        db_display = (it.get("db_display_name") or "").upper()
        if "GO" in dbname or "GO" in db_display or dbname == "GENE_ONTOLOGY":
            go_id = it.get("primary_id") or it.get("id") or it.get("display_id")
            if go_id:
                go_id = str(go_id).strip()
                desc = it.get("description") or it.get("display_id") or ""
                gos.append((go_id, desc))
    return gos

# ----- UniProt helpers -----
def get_uniprot_id_from_ensembl(ensembl_id):
    """Map Ensembl gene ID to UniProt ID using direct search."""
    if not ENABLE_UNIPROT_FALLBACK:
        return None

    # Try multiple search strategies
    search_queries = [
        f"database:ensembl AND {ensembl_id}",
        f"xref:ensembl-{ensembl_id}",
        f"gene:{ensembl_id}"
    ]

    for query in search_queries:
        try:
            url = f"{UNIPROT_REST}/uniprotkb/search"
            params = {
                "query": query,
                "format": "json",
                "size": 1
            }

            time.sleep(SLEEP_BETWEEN_UNIPROT)
            r = retry_get(url, params=params)
            if r and r.status_code == 200:
                data = r.json()
                results = data.get("results", [])
                if results:
                    return results[0].get("primaryAccession")
        except Exception:
            continue
    return None

def get_gene_symbol_from_uniprot(uniprot_id):
    """Get gene symbol from UniProt using UniProt ID."""
    if not ENABLE_UNIPROT_FALLBACK or not uniprot_id:
        return None

    url = f"{UNIPROT_REST}/uniprotkb/{uniprot_id}"
    params = {"fields": "genes"}

    try:
        time.sleep(SLEEP_BETWEEN_UNIPROT)
        r = retry_get(url, params=params)
        if r and r.status_code == 200:
            data = r.json()
            genes = data.get("genes", [])
            if genes:
                gene = genes[0]
                gene_name = gene.get("geneName", {}).get("value")
                if gene_name:
                    return gene_name
    except Exception:
        pass
    return None

def get_go_terms_from_uniprot(uniprot_id):
    """Get GO terms from UniProt using UniProt ID."""
    if not ENABLE_UNIPROT_FALLBACK or not uniprot_id:
        return []

    try:
        url = f"{UNIPROT_REST}/uniprotkb/{uniprot_id}"
        time.sleep(SLEEP_BETWEEN_UNIPROT)
        r = retry_get(url)
        if r and r.status_code == 200:
            data = r.json()
            gos = []
            cross_refs = data.get("uniProtKBCrossReferences", [])
            for xref in cross_refs:
                if xref.get("database") == "GO":
                    go_id = xref.get("id")
                    if go_id and go_id.startswith("GO:"):
                        description = ""
                        properties = xref.get("properties", [])
                        for prop in properties:
                            if prop.get("key") == "GoTerm":
                                description = prop.get("value", "")
                                break
                        gos.append((go_id, description))
            return gos
    except Exception:
        pass

    return []

# ----- NCBI Gene helpers -----
def get_ncbi_gene_id_from_ensembl(ensembl_id):
    """Map Ensembl gene ID to NCBI Gene ID using NCBI E-utilities."""
    if not ENABLE_NCBI_FALLBACK:
        return None

    try:
        search_url = f"{NCBI_EUTILS}/esearch.fcgi"
        params = {
            "db": "gene",
            "term": f"{ensembl_id}[Ensembl]",
            "retmode": "json",
            "retmax": 1
        }

        time.sleep(SLEEP_BETWEEN_NCBI)
        r = retry_get(search_url, params=params)
        if r and r.status_code == 200:
            data = r.json()
            id_list = data.get("esearchresult", {}).get("idlist", [])
            if id_list:
                return id_list[0]
    except Exception:
        pass
    return None

def get_gene_symbol_from_ncbi(ncbi_gene_id):
    """Get gene symbol from NCBI Gene using NCBI Gene ID."""
    if not ENABLE_NCBI_FALLBACK or not ncbi_gene_id:
        return None

    try:
        summary_url = f"{NCBI_EUTILS}/esummary.fcgi"
        params = {
            "db": "gene",
            "id": ncbi_gene_id,
            "retmode": "json"
        }

        time.sleep(SLEEP_BETWEEN_NCBI)
        r = retry_get(summary_url, params=params)
        if r and r.status_code == 200:
            data = r.json()
            result = data.get("result", {}).get(ncbi_gene_id, {})
            return result.get("nomenclature_symbol") or result.get("name")
    except Exception:
        pass
    return None

def get_go_terms_from_ncbi(ncbi_gene_id):
    """Get GO terms from NCBI Gene using NCBI Gene ID."""
    if not ENABLE_NCBI_FALLBACK or not ncbi_gene_id:
        return []

    try:
        summary_url = f"{NCBI_EUTILS}/esummary.fcgi"
        params = {
            "db": "gene",
            "id": ncbi_gene_id,
            "retmode": "json"
        }

        time.sleep(SLEEP_BETWEEN_NCBI)
        r = retry_get(summary_url, params=params)
        if r and r.status_code == 200:
            data = r.json()
            result = data.get("result", {}).get(ncbi_gene_id, {})
            go_terms = []
            for field in ["go_component", "go_function", "go_process"]:
                if field in result:
                    terms = result[field]
                    if isinstance(terms, list):
                        for term in terms:
                            if isinstance(term, dict) and "value" in term:
                                go_id = term["value"]
                                description = term.get("label", "")
                                if go_id.startswith("GO:"):
                                    go_terms.append((go_id, description))
            return go_terms
    except Exception:
        pass
    return []

def is_annotation_complete(annotation):
    """Check if annotation has both gene symbol and GO terms."""
    has_symbol = bool(annotation.get("gene_symbol", "").strip())
    has_go_terms = bool(annotation.get("go_ids", []))
    return has_symbol and has_go_terms

def needs_fallback(annotation):
    """Check if annotation needs UniProt fallback (missing symbol OR GO terms)."""
    has_symbol = bool(annotation.get("gene_symbol", "").strip())
    has_go_terms = bool(annotation.get("go_ids", []))
    return not has_symbol or not has_go_terms

# ----- core processing (no file I/O) -----
def annotate_ensembl_ids(id_list):
    """
    Accept list of Ensembl IDs and return a dict:
      {
        "annotations": [ {ensembl_id, gene_symbol, go_ids, go_terms}, ... ],
        "gene_symbols": [...],
        "go_ids": [...],
        "meta": {...}
      }
    No files are written.
    """
    ids = [str(x).strip() for x in id_list if x and str(x).strip()]
    ids = ids[:MAX_IDS]  # truncate for safety
    annotations = []
    gene_symbols_seen = []
    all_go_ids = set()

    for enid in ids:
        # small polite delay to avoid hammering Ensembl
        time.sleep(SLEEP_BETWEEN)
        symbol = None
        try:
            symbol = get_gene_symbol(enid)
        except Exception:
            symbol = None

        # another short delay before xrefs
        time.sleep(SLEEP_BETWEEN)
        gos = []
        try:
            gos = get_go_xrefs(enid)
        except Exception:
            gos = []

        go_ids = []
        go_terms = []
        for gid, desc in gos:
            if isinstance(gid, str) and gid.upper().startswith("GO:"):
                gid_up = gid.upper()
                go_ids.append(gid_up)
                go_terms.append(desc or "")
                all_go_ids.add(gid_up)

        # Create initial annotation
        annotation = {
            "ensembl_id": enid,
            "gene_symbol": symbol or "",
            "go_ids": go_ids,
            "go_terms": go_terms
        }

        # Check if we need fallback
        if needs_fallback(annotation):
            # Try UniProt first
            if ENABLE_UNIPROT_FALLBACK:
                uniprot_id = None
                try:
                    uniprot_id = get_uniprot_id_from_ensembl(enid)
                except Exception:
                    uniprot_id = None

                if uniprot_id:
                    annotation["_uniprot_fallback_used"] = True

                    # Try to get missing gene symbol from UniProt
                    if not symbol:
                        try:
                            uniprot_symbol = get_gene_symbol_from_uniprot(uniprot_id)
                            if uniprot_symbol:
                                annotation["gene_symbol"] = uniprot_symbol
                                annotation["_uniprot_symbol_added"] = True
                                symbol = uniprot_symbol
                        except Exception:
                            pass

                    # Try to get missing GO terms from UniProt
                    if not go_ids:
                        try:
                            uniprot_gos = get_go_terms_from_uniprot(uniprot_id)
                            if uniprot_gos:
                                annotation["_uniprot_go_added"] = True
                            for gid, desc in uniprot_gos:
                                if isinstance(gid, str) and gid.upper().startswith("GO:"):
                                    gid_up = gid.upper()
                                    go_ids.append(gid_up)
                                    go_terms.append(desc or "")
                                    all_go_ids.add(gid_up)
                            annotation["go_ids"] = go_ids
                            annotation["go_terms"] = go_terms
                        except Exception:
                            pass

            # Try NCBI Gene as second fallback if still incomplete
            if ENABLE_NCBI_FALLBACK and needs_fallback(annotation):
                ncbi_gene_id = None
                try:
                    ncbi_gene_id = get_ncbi_gene_id_from_ensembl(enid)
                except Exception:
                    ncbi_gene_id = None

                if ncbi_gene_id:
                    annotation["_ncbi_fallback_used"] = True

                    # Try to get missing gene symbol from NCBI
                    if not symbol:
                        try:
                            ncbi_symbol = get_gene_symbol_from_ncbi(ncbi_gene_id)
                            if ncbi_symbol:
                                annotation["gene_symbol"] = ncbi_symbol
                                annotation["_ncbi_symbol_added"] = True
                                symbol = ncbi_symbol
                        except Exception:
                            pass

                    # Try to get missing GO terms from NCBI
                    if not go_ids:
                        try:
                            ncbi_gos = get_go_terms_from_ncbi(ncbi_gene_id)
                            if ncbi_gos:
                                annotation["_ncbi_go_added"] = True
                            for gid, desc in ncbi_gos:
                                if isinstance(gid, str) and gid.upper().startswith("GO:"):
                                    gid_up = gid.upper()
                                    go_ids.append(gid_up)
                                    go_terms.append(desc or "")
                                    all_go_ids.add(gid_up)
                            annotation["go_ids"] = go_ids
                            annotation["go_terms"] = go_terms
                        except Exception:
                            pass

        annotations.append(annotation)
        if symbol:
            gene_symbols_seen.append(symbol)

    # Count fallback usage and clean up metadata flags
    uniprot_fallback_used = 0
    ncbi_fallback_used = 0
    uniprot_symbols_added = 0
    uniprot_go_terms_added = 0
    ncbi_symbols_added = 0
    ncbi_go_terms_added = 0

    for annotation in annotations:
        if annotation.get("_uniprot_fallback_used"):
            uniprot_fallback_used += 1
        if annotation.get("_ncbi_fallback_used"):
            ncbi_fallback_used += 1
        if annotation.get("_uniprot_symbol_added"):
            uniprot_symbols_added += 1
        if annotation.get("_uniprot_go_added"):
            uniprot_go_terms_added += 1
        if annotation.get("_ncbi_symbol_added"):
            ncbi_symbols_added += 1
        if annotation.get("_ncbi_go_added"):
            ncbi_go_terms_added += 1
        # Clean up metadata fields from final output
        annotation.pop("_uniprot_fallback_used", None)
        annotation.pop("_ncbi_fallback_used", None)
        annotation.pop("_uniprot_symbol_added", None)
        annotation.pop("_uniprot_go_added", None)
        annotation.pop("_ncbi_symbol_added", None)
        annotation.pop("_ncbi_go_added", None)

    result = {
        "annotations": annotations,
        "gene_symbols": sorted({s for s in gene_symbols_seen if s}),
        "go_ids": sorted(all_go_ids),
        "meta": {
            "version": VERSION,
            "count_input": len(id_list),
            "count_processed": len(ids),
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "uniprot_fallback": {
                "enabled": ENABLE_UNIPROT_FALLBACK,
                "fallback_used_count": uniprot_fallback_used,
                "symbols_added_from_uniprot": uniprot_symbols_added,
                "go_terms_added_from_uniprot": uniprot_go_terms_added
            },
            "ncbi_fallback": {
                "enabled": ENABLE_NCBI_FALLBACK,
                "fallback_used_count": ncbi_fallback_used,
                "symbols_added_from_ncbi": ncbi_symbols_added,
                "go_terms_added_from_ncbi": ncbi_go_terms_added
            }
        }
    }
    return result

# ----- HTTP endpoints -----
@app.route("/health", methods=["GET"])
def health():
    return jsonify({"status": "ok", "time": datetime.utcnow().isoformat() + "Z"})

@app.route("/version", methods=["GET"])
def version():
    return jsonify({"version": VERSION, "started": datetime.utcnow().isoformat() + "Z"})

@app.route("/config", methods=["GET"])
def config():
    return jsonify({
        "uniprot_fallback_enabled": ENABLE_UNIPROT_FALLBACK,
        "max_ids": MAX_IDS,
        "ensembl_rest_url": ENSEMBL_REST,
        "uniprot_rest_url": UNIPROT_REST,
        "version": VERSION
    })

@app.route("/annotate", methods=["POST", "GET"])
def annotate():
    """
    Accepts:
      - POST JSON: { "ids": ["ENSG...","ENSG..."] }  OR  { "id1":"...", "id2":"..." }
      - GET    : /annotate?id1=...&id2=...
    Returns JSON with whatever annotation data is computed.
    """
    ids = []
    if request.method == "POST":
        try:
            data = request.get_json(force=True, silent=True) or {}
        except Exception:
            data = {}
        if isinstance(data, dict):
            if "ids" in data and isinstance(data["ids"], list):
                ids = data["ids"]
            else:
                id1 = data.get("id1") or data.get("ensembl1")
                id2 = data.get("id2") or data.get("ensembl2")
                if id1:
                    ids.append(id1)
                if id2:
                    ids.append(id2)
    # GET fallback
    if not ids:
        id1 = request.args.get("id1") or request.args.get("ensembl1")
        id2 = request.args.get("id2") or request.args.get("ensembl2")
        if id1:
            ids.append(id1)
        if id2:
            ids.append(id2)

    if not ids:
        return jsonify({"error": "No Ensembl IDs provided. Provide JSON {'ids': [...] } or id1/id2 params."}), 400

    if len(ids) > MAX_IDS:
        return jsonify({"error": f"Too many IDs (limit {MAX_IDS}). For large jobs use batch mode."}), 400

    try:
        result = annotate_ensembl_ids(ids)
    except Exception as e:
        return jsonify({"error": "Annotation failed", "detail": str(e)}), 500

    return jsonify(result)

# Serve index.html from project root
@app.route("/", methods=["GET"])
def home():
    root = os.path.abspath(os.path.dirname(__file__) or ".")
    index_path = os.path.join(root, "index.html")
    if os.path.exists(index_path):
        return send_from_directory(root, "index.html")
    return jsonify({"error": "index.html not found"}), 404

@app.route("/index.html", methods=["GET"])
def index_html():
    return home()

# ----- CLI behavior: process file or run server -----
if __name__ == "__main__":
    # legacy CLI mode: `python ensembl.py some_ids.txt`
    if len(sys.argv) == 2 and os.path.exists(sys.argv[1]):
        path = sys.argv[1]
        with open(path, "r", encoding="utf-8") as fh:
            ids = [line.strip() for line in fh if line.strip()]
        print(json.dumps(annotate_ensembl_ids(ids), indent=2))
        sys.exit(0)

    # normal server mode
    print(f"Starting Ensembl annotation server (no file output) on http://127.0.0.1:5000 — version {VERSION}")
    # debug=False by default; use an auto-reload dev loop externally if you want hot reload
    app.run(host="0.0.0.0", port=5000, debug=False)
