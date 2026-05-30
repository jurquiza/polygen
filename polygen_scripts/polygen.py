from flask import Flask, redirect, url_for, render_template, request, Response, send_from_directory, session as flask_session
from io import BytesIO, StringIO
from zipfile import ZipFile
from datetime import date,time
import os, tempfile
import copy
import secrets
import threading

from werkzeug.local import LocalProxy

from engine import *


app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'ee839687a282e6493054c86e00e86925dc04de931acf4aed')
APP_DIR = os.path.dirname(os.path.abspath(__file__))

_SESSION_STORE = {}
_SESSION_LOCK = threading.Lock()


def _session_id():
    if 'polygen_sid' not in flask_session:
        flask_session['polygen_sid'] = secrets.token_urlsafe(32)
    return flask_session['polygen_sid']


def _user_session():
    sid = _session_id()
    with _SESSION_LOCK:
        return _SESSION_STORE.setdefault(sid, {})


session = LocalProxy(_user_session)


def set_session_defaults(defaults):
    for key,value in defaults.items():
        if key not in session:
            session[key] = copy.deepcopy(value)


def set_ptg_defaults():
    set_session_defaults({
        'msg': None,
        'PTG_name': '',
        'oligo_prefix': '',
        'oligo_index': '',
        'bb_linkers': '',
        'ad_linkers': '',
        'PTG_transfer': '',
        'poltype': 'ptg',
        'enzm': 'bbsi',
        'tm_range': [55, 65],
        'staticBorderPrimers': False,
        'noBorderPrimers': False,
        'clr': {'sequence_spacers': '#FFFFFF', 'link': '#FFFFFF', 'poltype_input': '#FFFFFF', 'oligo_index': '#FFFFFF', 'PTG_name': '#FFFFFF'}
    })


def set_peg_defaults():
    set_session_defaults({
        'msg': None,
        'PEG_sequence': '',
        'PEG_edits': '',
        'clr': {'sequence': '#FFFFFF', 'edits': '#FFFFFF'}
    })


def set_tas_defaults():
    set_session_defaults({
        'msg': None,
        'clr': {'exact_spacer': '#FFFFFF', 'scaffold': '#FFFFFF'},
        'TAS_system': 'TasH',
        'TAS_sequence': '',
        'TAS_edge_5': TAS_SYSTEMS['TasH']['edge_5'],
        'TAS_loop': TAS_SYSTEMS['TasH']['loop'],
        'TAS_edge_3': TAS_SYSTEMS['TasH']['edge_3'],
        'TAS_exact_spacer': '',
        'TAS_guides': []
    })

@app.route("/")
def home():
    return render_template("main.html")


@app.route("/learn")
def learn():
    return render_template("learn_more.html")


@app.route("/ptg", methods=["POST","GET"])
def sequence():

    ## Setting up defaults and variables
    set_ptg_defaults()
    session['clr'] = {'sequence_spacers': '#FFFFFF', 'link': '#FFFFFF', 'poltype_input': '#FFFFFF', 'oligo_index': '#FFFFFF', 'PTG_name': '#FFFFFF'}
    enzms={'bsai': ['gaggtctcg', 'cgagacctc'], 'bsmbi': ['tgcgtctca', 'tgagacgca'], 'btgzi': ['ctgcgatggagtatgtta', 'taacatactccatcgcag'], 'bbsi': ['ttgaagactt', 'aagtcttcaa']} #templates found in pUU080 (bsai), pUPD2 (bsmbi), Ortega-Escalante et al. 2018 (btgzi), Kun (bbsi)
    
    if request.method == "POST":
        
        if request.form['submitPTG'] == 'submit':
            
            ## Pulling all input from the page
            session['msg'] = None
            session['poltype'] = request.form["poltype_input"]
            session['enzm'] = request.form['enzm_input']
            session['PTG_name'] = request.form["PTG_name"].replace(" ", "_")
            session['oligo_prefix'] = request.form["oligo_prefix"]
            session['oligo_index'] = request.form["oligo_index"]
            session['PTG_transfer'] = request.form["sequence_spacers"]
            session['bb_linkers'] = request.form["bb_linkers"]
            session['ad_linkers'] = request.form["ad_linkers"]
            session['tm_range'] = [55, 65]
            session['staticBorderPrimers'] = request.form.get('staticBorderPrimers') if request.form.get('staticBorderPrimers') else False
            session['noBorderPrimers'] = request.form.get('noBorderPrimers') if request.form.get('noBorderPrimers') else False
            if session['poltype'] != 'ptg':
                session['staticBorderPrimers'] = False
                session['noBorderPrimers'] = False
            
            if not session['PTG_name']:
                session['PTG_name'] = request.form['poltype_input'].upper()
            if not session['oligo_prefix']:
                session['oligo_prefix'] = 'o'
            if not session['oligo_index']:
                session['oligo_index'] = '0'
            
            
            ## Catching input errors
            if '/' in session['PTG_name']:
                raise InvalidUsage("Polycistron name may not contain a '/'", status_code=400, payload={'pge': 'sequence.html', 'box': 'PTG_name'})
            if session['PTG_transfer'] == '':
                raise InvalidUsage("You must specify a polycistron description", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
            if session['oligo_index'] != '' and re.search(r'^[0-9]*$', session['oligo_index']) is None:
                raise InvalidUsage("Starting index must be a number", status_code=400, payload={'pge': 'sequence.html', 'box': 'oligo_index'})
            
            bbLinkersSplit = session['bb_linkers'].split(';') if session['bb_linkers'] else ['tgcc', 'gttt']
            adLinkersSplit = session['ad_linkers'].split(';') if session['ad_linkers'] else []
            for lnk in bbLinkersSplit+adLinkersSplit:
                if len(lnk) != 4 or re.search(r'^[ACGTacgt]*$', lnk) is None:
                    raise InvalidUsage("Invalid linker input", status_code=400, payload={'pge': 'sequence.html', 'box': 'link'})
            
            
            ## Preprocessing input and catching input errors. The error catching might be duplicated in the individual functions
            ## but catching the errors earlier saves computing time.
            PTG_input = session['PTG_transfer'].replace(" ", "").replace("\r\n", "").split('|')
            PTG_structure = []
            for element in PTG_input:
                e = element.split(';')
                if len(e) != 2:
                    raise InvalidUsage("Invalid input syntax", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                elif e[0] not in ['gRNA', 'pegRNA', 'smRNA', 'crRNA', 'tigRNA']:
                    raise InvalidUsage("Invalid RNA type", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                elif re.search(r'^[ACGTacgt]*$', e[1]) is None:
                    raise InvalidUsage("Invalid sequence input", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                else:
                    PTG_structure.append(e)
            
            
            ## Running inputs through engine functions   
            PTG = PTGbldr(session['PTG_name'], PTG_structure, session['poltype'])
            
            ggArgs = {'poltype': session['poltype'],
                      'enzm': session['enzm'],
                      'tm_range': session['tm_range'],
                      'bb_linkers': bbLinkersSplit,
                      'ad_linkers': adLinkersSplit}
            polycistron = scarless_gg(PTG, **ggArgs)
            
            primerArgs = {'oligo_prefix': session['oligo_prefix'], 
                          'oligo_index': session['oligo_index'], 
                          'staticBorderPrimers': session['staticBorderPrimers'],
                          'noBorderPrimers': session['noBorderPrimers'], 
                          'poltype': session['poltype'], 
                          'enzm': session['enzm'],
                          'bb_linkers': ggArgs['bb_linkers'],
                          'ad_linkers': ggArgs['ad_linkers']}
            session['plcstrn'] = annotatePrimers(polycistron, **primerArgs)
            
            if session['plcstrn'].warning:
                session['msg'] = session['plcstrn'].warning
            
            return render_template("primer_list.html", session=session)
            
            
        ## If the reset button is pressed, erase all input
        elif request.form['submitPTG'] == 'reset':
            session['PTG_name'] = session['oligo_prefix'] = session['oligo_index'] = session['bb_linkers'] = session['ad_linkers'] = session['PTG_transfer'] = ''
            return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None), session=session)
            
    else:
        
        return render_template("sequence.html",PTG_transfer=session.get('PTG_transfer', None), session=session)


@app.route("/peg", methods=["POST","GET"])
def peg_generation():
    
    ## Setting up defaults and variables
    set_peg_defaults()
    session['msg'] = None
    session['clr'] = {'sequence': '#FFFFFF', 'edits': '#FFFFFF'}

    if request.method == "POST":

        if 'file' in request.files:

            file = request.files['file'] # Retrieve file
            fileFormat = file.filename.split(".")[-1] # Retrieve file format
            stringio = StringIO(file.getvalue().decode("utf-8")) # Decode file
            record = SeqIO.read(stringio, fileFormat) # Read file
            session['PEG_sequence'] = str(record.seq) # Extract and save sequence

            return render_template("peg_generation.html", PEG_transfer=session.get('PEG_sequence', None), session=session)

        if request.form['submitPEG'] == 'submit':
        
            ## Pulling all inputs from page
            session['PEG_sequence'] = request.form["sequence"]
            session['PEG_edits'] = request.form["edits"]
            PEG_mode = request.form["mode"]

            ## Preprocessing input
            # Clean up whitespaces and newlines in sequence
            session['PEG_sequence'] = session['PEG_sequence'].replace(' ', '')
            session['PEG_sequence'] = session['PEG_sequence'].replace('\r\n', '')

            # Concatenate edits
            PEG_edits = session['PEG_edits'].split('|')
            pegs_list = []
            for edt in PEG_edits:
                pegs_list.append(edt.split(';'))

            ## Running inputs through engine function
            edits_list,session['msg'] = pegbldr(session['PEG_sequence'], pegs_list, PEG_mode)

            ## Postprocessing output
            session['PTG_transfer'] = ''
            for c,peg in enumerate(edits_list):
                session['PTG_transfer'] += str(peg[1]) + ';' +  str(peg[2])
                if c < len(edits_list)-1:
                    session['PTG_transfer'] += '|'

            return redirect(url_for('sequence'))

        ## If the reset button is pressed, erase all input
        elif request.form['submitPEG'] == 'reset':
            session['PEG_sequence'] = session['PEG_edits'] = ''
            return render_template("peg_generation.html", PEG_transfer=session.get('PEG_sequence', None), session=session)

        else:

            return render_template("main.html")

    else:
        
        return render_template("peg_generation.html", PEG_transfer=session.get('PEG_sequence', None), session=session)


@app.route("/tas", methods=["POST","GET"])
def tas_generation():

    set_tas_defaults()
    session['msg'] = None
    session['clr'] = {'exact_spacer': '#FFFFFF', 'scaffold': '#FFFFFF'}
    session['TAS_guides'] = session.get('TAS_guides', [])

    if request.method == "POST":

        if request.form['submitTAS'] == 'submit':
            session['TAS_system'] = request.form['system']
            defaults = TAS_SYSTEMS[session['TAS_system']]
            session['TAS_sequence'] = ''
            session['TAS_edge_5'] = request.form['edge_5'] or defaults['edge_5']
            session['TAS_loop'] = request.form['loop'] or defaults['loop']
            session['TAS_edge_3'] = request.form['edge_3'] or defaults['edge_3']
            session['TAS_exact_spacer'] = request.form['exact_spacer']

            guides, msg = tas_guide_design('',
                                           session['TAS_system'],
                                           session['TAS_edge_5'],
                                           session['TAS_loop'],
                                           session['TAS_edge_3'],
                                           20,
                                           25,
                                           75,
                                           session['TAS_exact_spacer'])
            session['TAS_guides'] = guides
            session['msg'] = msg
            session['PTG_transfer'] = '|'.join([guide['ptg_part'] for guide in guides])
            return render_template("tas_generation.html", TAS_transfer=session.get('TAS_sequence', None), session=session, tas_systems=TAS_SYSTEMS)

        elif request.form['submitTAS'] == 'send_ptg':
            selected = request.form.getlist('selected_guides')
            guides = [guide for guide in session.get('TAS_guides', []) if guide['name'] in selected]
            if guides == []:
                guides = session.get('TAS_guides', [])
            session['PTG_transfer'] = '|'.join([guide['ptg_part'] for guide in guides])
            session['poltype'] = 'tigRNA'
            return redirect(url_for('sequence'))

        elif request.form['submitTAS'] == 'reset':
            for key in ['TAS_sequence', 'TAS_exact_spacer', 'TAS_guides']:
                session[key] = '' if key == 'TAS_sequence' else []
            session['TAS_exact_spacer'] = ''
            return render_template("tas_generation.html", TAS_transfer=session.get('TAS_sequence', None), session=session, tas_systems=TAS_SYSTEMS)

    return render_template("tas_generation.html", TAS_transfer=session.get('TAS_sequence', None), session=session, tas_systems=TAS_SYSTEMS)
      
        
@app.route("/primer_list", methods=["POST","GET"])
def serve_primers():
    set_ptg_defaults()
    if 'plcstrn' not in session:
        session['msg'] = 'No oligo design is available for download yet.'
        return redirect(url_for('sequence'))
    
    ## Starting to write the csv file
    csv = 'List of oligos:\n'
    csv += 'oligo_ID,sequence\n'
    
    ## Appending the oligo list to the csv file
    for oligo in session['plcstrn'].oligos:
        csv += oligo[0] + ',' + oligo[1] + '\n'
        
    ## Appending the fragment table to the csv file
    csv += '\nTable of fragments:\n'
    csv += 'fragment_id,fragment_type,forward_primer,Tm_forw,reverse_primer,Tm_rev\n'
    for c,fragment in enumerate(session['plcstrn'].parts):
        csv += session['PTG_name']+'_f'+str(c)+','+fragment.type+','+session['plcstrn'].oligos[c*2][0]+','+str(np.round(fragment.primer_forward_tm,1))+','+session['plcstrn'].oligos[c*2+1][0]+','+str(np.round(fragment.primer_reverse_tm, 1))+'\n'
    
    ## Appending the input parameters to the csv file
    csv += '\nInput parameters:\n'
    csv += 'Polycistron name:,' + session['PTG_name'] + '\n'
    csv += 'Oligos prefix:,' + session['oligo_prefix'] + '\n'
    csv += 'Starting index:,' + session['oligo_index'] + '\n'
    csv += 'Border linkers:,' + '+'.join(session['bb_linkers'] if session['bb_linkers'] else ['tgcc', 'gttt']) + '\n'
    csv += 'Addidional linkers:,' + '+'.join(session['ad_linkers'] if session['ad_linkers'] else []) + '\n'
    csv += 'Type of Polycistron:,' + session['poltype'] + '\n'
    csv += 'Restriction enzyme:,' + session['enzm'] + '\n'
    csv += 'Melting temperature range:,' + str(session['tm_range'][0]) + '-' + str(session['tm_range'][1]) + '\n'
    csv += 'Invariable border primers:,' + str(session['staticBorderPrimers']) + '\n'
    csv += 'Omit border primers:,' + str(session['noBorderPrimers']) + '\n'
    
    csv += '\nSequence input:\n'
    PTG_input = session['PTG_transfer'].replace(" ", "").replace("\r\n", "").split("|")
    for part in PTG_input:
        csv += part.replace(";", ",") + '\n'
    
    ## Generating the object for the genbank file
    sr = SeqRecord(seq=Seq(session['plcstrn'].sequence, alphabet=IUPAC.ambiguous_dna), name=session['PTG_name'], annotations={'date': date.today().strftime("%d-%b-%Y").upper(), 'topology': 'linear'})
    for ftr in session['plcstrn'].features:
        sr.features.append(ftr)
    
    ## Jsonifying the raw output
    gb_json = {}
    gb_json['polycistron'] = polyToJson(session['plcstrn'])
    gb_json['msg'] = session['msg']
    
    ## Writing everything to a zip file
    in_memory = BytesIO()
    zf = ZipFile(in_memory, mode='w')
    zf.writestr(session['PTG_name']+"_oligos.csv", csv)
    zf.writestr(session['PTG_name']+".gb", sr.format('genbank'))
    zf.writestr(session['PTG_name']+'_raw.json', json.dumps(gb_json))
    with open(os.path.join(APP_DIR, 'protocol.txt')) as f:
        zf.writestr('protocol.txt', f.read())
    zf.close()
    in_memory.seek(0)
    outpt = in_memory.read()
    
    ## Returning the zip file
    return Response(
        outpt,
        mimetype="application/zip",
        headers={'Content-Disposition': 'attachment;filename=pg_%s_%s.zip' %(session['PTG_name'], str(date.today()))})


@app.route("/impressum")
def impress():
    return render_template("impressum.html")


@app.route("/success")
def success():
    return "Success"


@app.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    set_session_defaults({'clr': {}})
    session['msg'] = error.to_dict()['message']
    session['clr'][error.to_dict()['box']] = '#FA5858'
    return render_template(error.to_dict()['pge'],
                           PEG_transfer=session.get('PEG_sequence'),
                           TAS_transfer=session.get('TAS_sequence'),
                           PTG_transfer=session.get('PTG_transfer', None),
                           session=session,
                           tas_systems=TAS_SYSTEMS)


if __name__=='__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 5000))) ##host 0.0.0.0 is important for docker container
