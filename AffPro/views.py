from django.http import HttpResponse
from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from django.core.files.storage import FileSystemStorage
#from pymongo import MongoClient
from django.conf import settings
import datetime
from django.contrib.auth import authenticate, login, logout
from django.shortcuts import render, redirect
from django.contrib import messages
from django.contrib.auth.hashers import make_password, check_password
from django.conf import settings
from django.http import HttpResponse
#import pymongo
from django.contrib.auth.decorators import login_required
import csv


def aboutUS(request):
    return HttpResponse("<b>Welcome to AffPRo</b>")

def contactDetails(request, id):
    return HttpResponse(id)

def homePage(request):
    return render(request,"index.html")



#from rdkit import Chem
#from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import io
import base64
from rdkit.Chem import Descriptors

def generate_ligand_descriptors(smiles):
    # Convert the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    # Calculate the required descriptors
    descriptors = {
        "MolecularWeight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
        "NumRings": Descriptors.RingCount(mol),
        #"VanDerWaalsVolume": Descriptors.VanDerWaalsVolume(mol),  # Approximation for Van der Waals Volume
        #"BalabanIndex": Descriptors.BalabanIndex(mol),
    }

    return descriptors

from rdkit.Chem import Draw
def visualize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            buffer = io.BytesIO()
            img.save(buffer, format="PNG")
            img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
            buffer.close()
            return f'<img src="data:image/png;base64,{img_str}" alt="SMILES Visualization"/>'
        else:
            return "<p>Invalid SMILES string</p>"
    except Exception as e:
        return f"<p>Error generating visualization: {e}</p>"




def visualize_3D_smiles(smiles):
    try:
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # Add hydrogens
        AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
        AllChem.UFFOptimizeMolecule(mol)  # Optimize structure
        pdb_block = Chem.MolToPDBBlock(mol)  # Convert to PDB format

        # Generate 3D visualization using py3Dmol
        view = py3Dmol.view(width=700, height=300)
        view.addModel(pdb_block, "pdb")  # Add PDB model
        view.setStyle({'stick': {}})  # Use stick model
        view.zoomTo()  # Zoom to fit the molecule

        # Generate HTML for embedding in template
        return view._make_html()
    except Exception as e:
        return "<p>Invalid SMILES string</p>"
        
        #return f"<p style='color: red;'>Error visualizing SMILES: {e}</p>"



    
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def extract_protein_descriptors(sequence):
    try:
        analysis = ProteinAnalysis(sequence)
        descriptors = {
            "Protein Length": len(sequence),
            "Protein Molecular Weight": analysis.molecular_weight(),
            "Protein Aromaticity": analysis.aromaticity(),
            "Protein Instability Index": analysis.instability_index(),
            "Protein Isoelectric Point": analysis.isoelectric_point(),
            "Protein Gravy": analysis.gravy(),
            
        }

        # Amino acid composition
        amino_acid_percent = analysis.get_amino_acids_percent()
        for aa, percent in amino_acid_percent.items():
            descriptors[f"Protein Amino Acid Percent {aa}"] = percent

        return descriptors
    except Exception as e:
        print(f"Error processing protein sequence: {sequence}\n{e}")
        return {}




def PLInteraction(request):
    input1 = ''
    input2 = ''
    result = None
    smiles_visualization = None
    smiles_3D_visualization = None 
    ligand_descriptors= None
    protein_descriptors= None
    active_content_id = 'content-2'

    if request.method == 'POST':
        input1 = request.POST.get('input1', '')
        input2 = request.POST.get('input2', '')
        active_content_id = request.POST.get('activeContentId', 'content-2')
        result = f"Processed: SMILES = {input1}, Protein = {input2}"

        if input1:
            smiles_visualization = visualize_smiles(input1)
            smiles_3D_visualization=visualize_3D_smiles(input1)
            ligand_descriptors=generate_ligand_descriptors(input1)
            

        if input2:
            protein_descriptors=extract_protein_descriptors(input2)    

    return render(request, 'pl_int.html', {
        'input1': input1,
        'input2': input2,
        'result': result,
        'smiles_visualization': smiles_visualization,
        'smiles_3D_visualization':smiles_3D_visualization,
        'descriptors':ligand_descriptors,
        'active_content_id': active_content_id,
        'protein_descriptors':protein_descriptors,
    })


from django.core.files.storage import FileSystemStorage
from django.shortcuts import render





def process_input(request):
    if request.method == 'POST' and 'pdb_file' in request.FILES:
        # Handle the uploaded file
        pdb_file = request.FILES['pdb_file']
        fs = FileSystemStorage()  # Define where files are stored
        filename = fs.save(pdb_file.name, pdb_file)  # Save the file
        file_url = fs.url(filename)  # Get the file's URL if needed
        
        return render(request, 'pl_int.html', {
            'result_1': True,  # Pass this to indicate success
            'file_url': file_url,  # Optional: URL of the uploaded file
            'active_content_id': 'content-5',  # Preserve tab state
        })
    
    return render(request, 'pl_int.html', {
        'active_content_id': 'content-5',  # Preserve tab state
    })


from user_data_app.models import InteractionResult, InteractionResult_PP, InteractionResult_PL
from django.utils.timezone import now

def sidebar_dd(request):

    if not request.user.is_authenticated:  # Check if user is logged in
        messages.error(request, "Please log in first!")
        return redirect(f"{reverse('login')}?next={request.path}")

    

    username = request.user.username  # Get the logged-in username
    input1 = ''
    input2= ''
    input3=None
    result_dd = None
    result_2_dd=None
    csv_data=[]
    top_n=[]
    smiles_visualization_1 = None
    smiles_3D_visualization_1 = None 
    smiles_visualization_2 = None
    smiles_3D_visualization_2 = None 
    ligand_descriptors_1= None
    ligand_descriptors_2= None



    if request.method == 'POST':
        input1 = request.POST.get('input1_dd', '')
        input2 = request.POST.get('input2_dd', '')
        input3 = request.FILES.get('csv_file_dd') 
        

        if input1:
            smiles_visualization_1 = visualize_smiles(input1)
            smiles_3D_visualization_1=visualize_3D_smiles(input1)
            ligand_descriptors_1=generate_ligand_descriptors(input1)

        if input2:
            smiles_visualization_2 = visualize_smiles(input2)
            smiles_3D_visualization_2=visualize_3D_smiles(input2)
            ligand_descriptors_2=generate_ligand_descriptors(input2) 

        if (input1 and input2):
            #result_dd=str(dd_results(input1, input2).item())
            result_dd=str(dd_results(input1, input2))
            result_dd = "Interacting" if result_dd == "1" else "Non Interacting"
            print(result_dd)
            # Save result under the logged-in user
        
            InteractionResult.objects.create(
                user=request.user,  # Automatically associate the result with the logged-in user
                input1=input1,
                input2=input2,
                result=result_dd,
                timestamp=now()
            )

        

        if input3:
                # Process the CSV file and obtain table_data
                
                csv_data = process_csv_file(input3, request)
                
                result_2_dd= prediction_dd_batch(input3)
                top_n=get_top_n_min_predicted_dd(result_2_dd)
                
        else:
                messages.error(request, "CSV file not provided.")

    user_results = InteractionResult.objects.filter(user=request.user).order_by("-timestamp")
    
    data = [
        {
            "input1_user": r.input1,
            "input2_user": r.input2,
            "result_user": r.result,
            "timestamp_user": r.timestamp.strftime("%Y-%m-%d %H:%M:%S")
        } 
        for r in user_results
    ]


           

    return render(request, 'index_sidebar_dd.html', {
        'input1': input1,
        'input2': input2,
        'result_dd': result_dd,
        'result_2_dd':result_2_dd,
        'smiles_visualization_1': smiles_visualization_1,
        'smiles_3D_visualization_1':smiles_3D_visualization_1,
        'smiles_visualization_2': smiles_visualization_2,
        'smiles_3D_visualization_2':smiles_3D_visualization_2,
        'descriptors_1':ligand_descriptors_1,
        'descriptors_2':ligand_descriptors_2,
        'username':username,
        "csv_data":csv_data,
        "top_n":top_n,
        'data':data,
    })

from django.http import JsonResponse
from rdkit import Chem
import json

def validate_smiles(request):
    if request.method == "POST":
        data = json.loads(request.body)
        smiles = data.get('smiles', '')
        valid = Chem.MolFromSmiles(smiles) is not None
        return JsonResponse({'valid': valid})
    return JsonResponse({'valid': False})

# Define valid amino acid single-letter codes
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def validate_protein_sequence(request):
    if request.method == "POST":
        data = json.loads(request.body)
        sequence = data.get('protein', '').upper().strip()

        if not sequence:
            return JsonResponse({'valid': False, 'error': 'No sequence provided.'})

        is_valid = all(aa in VALID_AMINO_ACIDS for aa in sequence)

        return JsonResponse({'valid': is_valid})
    
    return JsonResponse({'valid': False, 'error': 'Invalid request method.'})




def sidebar_pl(request):

    if not request.user.is_authenticated:  # Check if user is logged in
        messages.error(request, "Please log in first!")
        return redirect(f"{reverse('login')}?next={request.path}")
    

    username = request.user.username  # Get the logged-in username

    # Initialize inputs and results
    input1 = ''  # SMILES input (string)
    input2 = ''  # Protein input (possibly file path or string)
    input3 = None
    input4 = None
    csv_data = []
    top_n=[]
 
    result = None
    smiles_visualization = None
    smiles_3D_visualization = None
    ligand_descriptors = None
    protein_descriptors = None
    pdb_visual = "Load the PDB file.."
    pdb_read="Load the PDB file.."

    viewer_style = "height:600px; width:700px; position:absolute; top:225px; left:1150px;"
    viewer_class = "viewer_3Dmoljs"
    data_pdb = ""  # Set your PDB data if available
    data_backgroundcolor = "0xffffff"
    data_style = "ballstick"
    data_ui = "true"
   
    result_2= None
    # model, model_esm, batch_converter =load_models()    
    

    

    if request.method == 'POST':
        
        # Get the hidden field to distinguish the forms.
        form_type = request.POST.get('form_type', 'project')  # default to project if not provided

        # For both forms, get the common fields (if any)
        input1 = request.POST.get('input1_pl', '')
        input2 = request.POST.get('input2_pl', '')
        # Files for project form and pdb form might be in different inputs
        input3 = request.FILES.get('input3_pl')  # From project form (optional)
        input4 = request.FILES.get('csv_file')  

        

        if form_type == 'project':
            # Process the project form
            # SMILES and protein descriptors

            
            if input1:
                smiles_visualization = visualize_smiles(input1)
                smiles_3D_visualization = visualize_3D_smiles(input1)
                ligand_descriptors = generate_ligand_descriptors(input1)
                    
               

            if input2:
                protein_descriptors = extract_protein_descriptors(input2)
                    
                #result= prediction_ki(model,input1, input2, model_esm, batch_converter)
                result= prediction_ki(input1, input2)

                        
        
            # Process PDB file from project form if provided
            if input3:
                pdb_visual = visualize_pdb(input3)
                input2 = str(get_chain_a_sequence(input3))
                protein_descriptors = extract_protein_descriptors(input2)
                pdb_read=parse_atom_records_with_units_html(input3)
                # result= prediction_ki(model,input1, input2, model_esm, batch_converter)
                result= prediction_ki(input1, input2)
            else:
                print("Project form: Not Uploaded file for input3_pl")

            InteractionResult_PL.objects.create(
                    user=request.user,  # Automatically associate the result with the logged-in user
                    input1=input1,
                    input2=input2,
                    result=result,
                    timestamp=now()
                )
           
        

        elif form_type == 'csv_form':
            # Process the CSV form: input4 holds the uploaded CSV file.
            if input4:
                # Process the CSV file and obtain table_data
                
                csv_data = process_csv_file(input4, request)
                
                #result_2= prediction_ki_batch(input4, model,  model_esm,  batch_converter)
                result_2= prediction_ki_batch(input4)

                top_n=get_top_n_min_predicted_ki(result_2)
                
            else:
                messages.error(request, "CSV file not provided.")


        


    user_results_pl = InteractionResult_PL.objects.filter(user=request.user).order_by("-timestamp")
    
    data_pl = [
        {
            "input1_user": r.input1,
            "input2_user": r.input2,
            "result_user": r.result,
            "timestamp_user": r.timestamp.strftime("%Y-%m-%d %H:%M:%S")
        } 
        for r in user_results_pl
    ]

    context = {
        'input1': input1,
        'input2': input2,
        'input4': input4,
        'result': result,
        'smiles_visualization': smiles_visualization,
        'smiles_3D_visualization': smiles_3D_visualization,
        'descriptors': ligand_descriptors,
        'protein_descriptors': protein_descriptors,
        'username': username,
        'pdb_visual': pdb_visual,
        "viewer_style": viewer_style,
        "viewer_class": viewer_class,
        "data_pdb": data_pdb,
        "data_backgroundcolor": data_backgroundcolor,
        "data_style": data_style,
        "data_ui": data_ui,
        "pdb_read":pdb_read,
        "csv_data":csv_data,
        "top_n":top_n,
        "result_2":result_2,
        "data_pl":data_pl,
    }
    
    return render(request, 'index_sidebar_pl.html', context)




def sidebar_pp(request):

    if not request.user.is_authenticated:  # Check if user is logged in
        messages.error(request, "Please log in first!")
        return redirect(f"{reverse('login')}?next={request.path}")

    
    username = request.user.username  # Get the logged-in username
    input1 = ''
    input2 = ''
    input3 = None
    input4 = None
    input5 = None
    csv_data = []
    top_n=[]
    pdb_visual_1 = "Load the PDB-I file.."
    pdb_read_1="Load the PDB-I file.."
    pdb_visual_2 = "Load the PDB-II file.."
    pdb_read_2="Load the PDB-II file.."

    result_pp = None
    result_pp_2=None
    smiles_visualization = None
    smiles_3D_visualization = None 
    protein_descriptors_1= None
    protein_descriptors_2= None

    viewer_style = "height:600px; width:700px; position:absolute; top:225px; left:1150px;"
    viewer_class = "viewer_3Dmoljs"
    data_pdb = ""  # Set your PDB data if available
    data_backgroundcolor = "0xffffff"
    data_style = "ballstick"
    data_ui = "true"
    

    if request.method == 'POST':
        input1 = request.POST.get('input1_pp', '')
        input2 = request.POST.get('input2_pp', '')
        input3 = request.FILES.get('input3_pp')  
        input4 = request.FILES.get('input4_pp')  
        input5 = request.FILES.get('csv_file_pp') 

        if input1:
            
            protein_descriptors_1=extract_protein_descriptors(input1)

        if input2:
            protein_descriptors_2=extract_protein_descriptors(input2) 

        if input3:
            pdb_visual_1 = visualize_pdb_pp(input3)
            input1 = str(get_chain_a_sequence(input3))
            protein_descriptors_1 = extract_protein_descriptors(input1)
            pdb_read_1= parse_atom_records_with_units_html(input3)
        else:
            print("Project form: Not Uploaded file for input3_pl")
            

        if input4:
            pdb_visual_2 = visualize_pdb_pp(input4)
            input2 = str(get_chain_a_sequence(input4))
            protein_descriptors_2 = extract_protein_descriptors(input2)
            pdb_read_2=parse_atom_records_with_units_html(input4)
        else:
            print("Project form: Not Uploaded file for input4_pl")


        if input5:
                # Process the CSV file and obtain table_data
                
                csv_data = process_csv_file(input5, request)
                result_pp_2= prediction_pp_batch(input5)
               
        else:
            messages.error(request, "CSV file not provided.")

        if (input1 and input2):
            result_pp=str(pp_results(input1, input2))
            result_pp = "Interacting" if result_pp == "1.0" else "Non Interacting"
            print(result_pp)

            InteractionResult_PP.objects.create(
                            user=request.user,  # Automatically associate the result with the logged-in user
                            input1=input1,
                            input2=input2,
                            result=result_pp,
                            timestamp=now()
                        )        
            
    user_results_pp = InteractionResult_PP.objects.filter(user=request.user).order_by("-timestamp")
    
    data_pp = [
        {
            "input1_user": r.input1,
            "input2_user": r.input2,
            "result_user": r.result,
            "timestamp_user": r.timestamp.strftime("%Y-%m-%d %H:%M:%S")
        } 
        for r in user_results_pp
    ]



    return render(request, 'index_sidebar_pp.html', {
        'input1': input1,
        'input2': input2,
        'result_pp': result_pp,
        'result_pp_2':result_pp_2,
        'smiles_visualization': smiles_visualization,
        'smiles_3D_visualization':smiles_3D_visualization,
        'protein_descriptors_1':protein_descriptors_1,
        'protein_descriptors_2':protein_descriptors_2,
        'username':username,
        "pdb_read_1":pdb_read_1,
        "pdb_read_2":pdb_read_2,
        "pdb_visual_1":pdb_visual_1,
        "pdb_visual_2":pdb_visual_2,
        "viewer_style": viewer_style,
        "viewer_class": viewer_class,
        "data_pdb": data_pdb,
        "data_backgroundcolor": data_backgroundcolor,
        "data_style": data_style,
        "data_ui": data_ui,
        "csv_data":csv_data,
        "top_n":top_n,
        "data_pp":data_pp,
    })


def dropdown_mol_str(request):
    input1=""
    smiles_visualization_1 = None
    smiles_3D_visualization_1 = None 
    
    if request.method == 'POST':
        input1 = request.POST.get('input1', '')

        if input1:
            smiles_visualization_1 = visualize_smiles_2(input1)
            smiles_3D_visualization_1=visualize_3D_smiles_2(input1)

    return render(request, 'dropdown_visualize_molstr.html', {'input1': input1,'smiles_visualization_1': smiles_visualization_1,
        'smiles_3D_visualization_1':smiles_3D_visualization_1,})
                                                              

def dropdown_prot_str(request):

    input1=""
    pdb_visual="Load the PDB File"
    

    if request.method == 'POST':
        input1 = request.FILES.get('input3_pl')

        if input1:
            pdb_visual = visualize_pdb_protstr(input1)
    
   
    return render(request, 'dropdown_visualize_protstr.html', {'pdb_visual':pdb_visual})

def dropdown_li_desc(request):
    input1=""
    descriptors=None 
   

    if request.method == 'POST':
        input1 = request.POST.get('input1', '')
        descriptors = generate_ligand_descriptors(input1)

    return render(request, 'dropdown_descriptor_li.html',{'input1':input1,'descriptors':descriptors})

def dropdown_pro_desc(request):
    input1=""
    input2=None
    protein_descriptors_1= None


    if request.method == 'POST':
        input1 = request.POST.get('input1', '')
        input2 = request.FILES.get('input3_pl')
        

        if input1:
            protein_descriptors_1 = extract_protein_descriptors(input1)

    
        if input2:
            input1 = str(get_chain_a_sequence(input2))
            protein_descriptors_1 = extract_protein_descriptors(input1)
    
    return render(request, 'dropdown_descriptor_pro.html', {'input1':input1, 'protein_descriptors_1':protein_descriptors_1})

def dropdown_read_pdb(request):
    
    input1=""
    pdb_read_1="Load the PDB file.."

    if request.method == 'POST':
        input1 = request.FILES.get('input3_pl')
    

        if input1:
                pdb_read_1= parse_atom_records_with_units_html(input1)
        else:
            print("Project form: Not Uploaded file for input3_pl")

    return render(request, 'dropdown_read_pdb.html', {'input1':input1, 'pdb_read_1':pdb_read_1})


from django.shortcuts import render, redirect
from django.contrib import messages
from django.contrib.auth.hashers import make_password, check_password



from django.contrib.auth import authenticate, login
from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse

def login_view(request):
    if request.user.is_authenticated:
        return redirect("home")  # Redirect if already logged in

    next_url = request.GET.get("next")  # Get 'next' parameter from URL

    if request.method == "POST":
        username = request.POST.get("username")
        password = request.POST.get("password")

        # Authenticate user using Django's built-in authentication
        user = authenticate(request, username=username, password=password)

        if user is not None:
            login(request, user)
            messages.success(request, "Login successful!")

            # Redirect to the requested page or home
            return redirect(next_url) if next_url else redirect("home")

        else:
            messages.error(request, "Invalid username or password!")
            return redirect(reverse("login") + f"?next={next_url}" if next_url else reverse("login"))

    return render(request, "login.html", {"next": next_url})


# Logout View
"""
def logout_view(request):
    if "username" in request.session:
        del request.session["username"]  # Clear session
        messages.success(request, "Logged out successfully!")
    return redirect("home")
"""

from django.contrib.auth import logout
from django.shortcuts import redirect
from django.contrib import messages

def logout_view(request):
    logout(request)  # Logs out the user
    messages.success(request, "Logged out successfully!")
    return redirect("home")  # Redirect to home page after logout


def visualize_pdb(pdbfile):
    """
    # Read the PDB file content
    with open(pdbfile, "r") as pdb_file:
        pdb_data = pdb_file.read()
    """

    pdb_data = pdbfile.read().decode('utf-8')

    # Create a py3Dmol viewer
    view = py3Dmol.view(width=800, height=600)

    # Load the PDB structure
    view.addModel(pdb_data, "pdb")  

    # Apply a cartoon representation
    view.setStyle({"cartoon": {"color": "spectrum"}})

    # Zoom to fit the structure
    view.zoomTo()

    embed_code = view._make_html()  

    return embed_code


def visualize_pdb_pp(pdbfile):
    """
    # Read the PDB file content
    with open(pdbfile, "r") as pdb_file:
        pdb_data = pdb_file.read()
    """

    pdb_data = pdbfile.read().decode('utf-8')

    # Create a py3Dmol viewer
    view = py3Dmol.view(width=800, height=300)

    # Load the PDB structure
    view.addModel(pdb_data, "pdb")  

    # Apply a cartoon representation
    view.setStyle({"cartoon": {"color": "spectrum"}})

    # Zoom to fit the structure
    view.zoomTo()

    embed_code = view._make_html()  

    return embed_code


def visualize_pdb_protstr(pdbfile):
    """
    # Read the PDB file content
    with open(pdbfile, "r") as pdb_file:
        pdb_data = pdb_file.read()
    """

    pdb_data = pdbfile.read().decode('utf-8')

    # Create a py3Dmol viewer
    view = py3Dmol.view(width=900, height=500)

    # Load the PDB structure
    view.addModel(pdb_data, "pdb")  

    # Apply a cartoon representation
    view.setStyle({"cartoon": {"color": "spectrum"}})

    # Zoom to fit the structure
    view.zoomTo()

    embed_code = view._make_html()  

    return embed_code



def get_chain_a_sequence(pdb_input):
    """
    Extracts the primary sequence for Chain A from a PDB file using SEQRES records.
    
    Args:
        pdb_input (str or file-like): Either a file path to the PDB file or a file-like object.
        
    Returns:
        str: The one-letter FASTA sequence for chain A. Returns an empty string if not found.
    """
    # Mapping from three-letter to one-letter amino acid codes
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    residues_chain_a = []
    
    # If pdb_input is a string (file path), open it; otherwise assume it's a file-like object.
    if isinstance(pdb_input, (str, bytes)):
        # Open using the file path
        with open(pdb_input, "r") as pdb_file:
            lines = pdb_file.readlines()
    else:
        # pdb_input is assumed to be a file-like object (e.g., InMemoryUploadedFile)
        pdb_input.seek(0)  # ensure we're at the beginning
        # Read and decode the content (assuming it's UTF-8 encoded)
        content = pdb_input.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8")
        lines = content.splitlines()
    
    for line in lines:
        if line.startswith("SEQRES"):
            # Chain identifier is in column 12 (index 11)
            chain_id = line[11].strip()
            if chain_id == 'A':
                # Residues start at column 19 (index 18) and are space-separated
                residues = line[18:].split()
                residues_chain_a.extend(residues)
    
    # Convert three-letter codes to one-letter codes
    chain_a_sequence = "".join(three_to_one.get(res.upper(), 'X') for res in residues_chain_a)
    return chain_a_sequence


from tabulate import tabulate
import os

def parse_atom_records_with_units_html(pdb_input):
    """
    Parses ATOM records from a PDB file and returns a structured HTML table,
    including the dimensions/units in the headers.
    
    Args:
        pdb_input (str or file-like object): A path to the PDB file or a file-like object.
    
    Returns:
        str: A formatted HTML table of ATOM records.
    """
    atom_data = []
    
    # Determine if pdb_input is a file path or a file-like object.
    if isinstance(pdb_input, (str, os.PathLike)):
        with open(pdb_input, "r") as pdb_file:
            lines = pdb_file.readlines()
    else:
        # Assume file-like object, ensure reading from the beginning.
        pdb_input.seek(0)
        content = pdb_input.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8")
        lines = content.splitlines()
    
    for line in lines:
        if line.startswith("ATOM"):
            # Using fixed-width slicing based on PDB format specifications.
            record_name  = line[0:6].strip()
            serial       = line[6:11].strip()
            atom_name    = line[12:16].strip()
            alt_loc      = line[16].strip()
            residue_name = line[17:20].strip()
            chain_id     = line[21].strip()
            residue_seq  = line[22:26].strip()
            insertion    = line[26].strip()
            x            = line[30:38].strip()  # Coordinates in Angstroms (Å)
            y            = line[38:46].strip()  # Coordinates in Angstroms (Å)
            z            = line[46:54].strip()  # Coordinates in Angstroms (Å)
            occupancy    = line[54:60].strip()
            temp_factor  = line[60:66].strip()  # Temperature factor in Å²
            element      = line[76:78].strip()
            charge       = line[78:80].strip()
            
            atom_data.append([
                record_name, serial, atom_name, residue_name,
                chain_id, residue_seq, x, y, z, occupancy, temp_factor,
                element, charge
            ])
    
    # Include units/dimensions in the headers.
    headers = [
        "Rec", "Serial", "Atom", "Res", "Chain", "ResSeq", 
        "X (Å)", "Y (Å)", "Z (Å)", "Occupancy", "TempFactor (Å²)",
        "Elem", "Charge"
    ]
    # Return the table formatted as HTML
    return tabulate(atom_data, headers=headers, tablefmt="html")




  












def dd_results(ligand1, ligand2, threshold=0.5):
    time.sleep(2)
    return random.randint(0, 1)
    

def prediction_dd_batch(csv_file):
    if not csv_file:
        print("⚠ No file uploaded!")
        return []
    
    csv_file.seek(0)  # Reset file pointer
    decoded_file = csv_file.read().decode("utf-8").strip()
    
    if not decoded_file:
        print("⚠ Error: File is empty!")
        return []

    print("🔹 Raw File Content:\n", repr(decoded_file))  # Debugging
    
    # Detect delimiter (comma vs. tab)
    delimiter = ',' if decoded_file.count(',') > decoded_file.count('\t') else '\t'
    print(f"📝 Detected delimiter: {repr(delimiter)}")

    reader = csv.reader(decoded_file.splitlines(), delimiter=delimiter)
    
    # Read the header
    header = next(reader, None)
    print("📝 Header:", header)

    if not header or len(header) < 2:
        print("⚠ Error: Header not detected or missing columns!")
        return []

    results = []
    
    for row in reader:
        print("🔸 Row Data:", row)  # Debugging Step
        
        if len(row) >= 2:
            val1 = row[0].strip()  # Ligand
            val2 = row[1].strip()  # Protein
            
            # Call the prediction function
            #result = str(dd_results(val1, val2).item())
            result = str(dd_results(val1, val2))
            result = "Interacting" if result == "1" else "Non Interacting"
            results.append({'ligand': val1, 'protein': val2, 'prediction': result})
    
    print("✅ Final Results:", results)
    return results


def get_top_n_min_predicted_dd(results):
    pass








def visualize_3D_smiles_2(smiles):
    try:
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # Add hydrogens
        AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
        AllChem.UFFOptimizeMolecule(mol)  # Optimize structure
        pdb_block = Chem.MolToPDBBlock(mol)  # Convert to PDB format

        # Generate 3D visualization using py3Dmol
        view = py3Dmol.view(width=900, height=300)
        view.addModel(pdb_block, "pdb")  # Add PDB model
        view.setStyle({'stick': {}})  # Use stick model
        view.zoomTo()  # Zoom to fit the molecule

        # Generate HTML for embedding in template
        return view._make_html()
    except Exception as e:
        return "<p>Invalid SMILES string</p>"
        #return f"<p style='color: red;'>Error visualizing SMILES: {e}</p>"
    



def visualize_smiles_2(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Specify desired image size (e.g., 300x300 pixels)
            img = Draw.MolToImage(mol, size=(300, 300))
            
            buffer = io.BytesIO()
            img.save(buffer, format="PNG")
            img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
            buffer.close()
            
            # You can also control the displayed size via HTML if you like
            return f'<img src="data:image/png;base64,{img_str}" alt="SMILES Visualization" width="500" height="350" />'
        else:
            return "<p>Invalid SMILES string</p>"
    except Exception as e:
        return f"<p>Error generating visualization: {e}</p>"
    


def process_csv_file(csv_file, request=None):
    """
    Processes an uploaded CSV file and returns the data as a list of rows.
    Optionally adds messages to the request if provided.
    """
    table_data = []  # This will hold rows from the CSV file
    if not csv_file:
        if request:
            messages.error(request, "No CSV file uploaded.")
        return table_data

    try:
        # Read and decode the CSV file into a list of lines
        decoded_file = csv_file.read().decode("utf-8").splitlines()
        reader = csv.reader(decoded_file)
        table_data = list(reader)
        if request:
            messages.success(request, "CSV file processed successfully!")
    except Exception as e:
        if request:
            messages.error(request, f"Error processing CSV file: {e}")
    return table_data
            


import pandas as pd

def get_top_n_min_predicted_ki(result_2, n=3):
    """
    Get the top N entries with the minimum 'prediction' values.
    
    Args:
    result_2 (list of dict): List containing ligand, protein, and prediction.
    n (int): Number of minimum prediction values to return.

    Returns:
    list of dict: Top N entries sorted by minimum prediction.
    """
    # Sort the list based on 'prediction' key
    sorted_results = sorted(result_2, key=lambda x: x["prediction"])
    
    # Return the top N minimum predictions
    return sorted_results[:n]




from rdkit import Chem
from rdkit.Chem import rdchem

import random

import time
def prediction_ki(smiles, sequence):
    
    time.sleep(2)
    return random.uniform(1, 1000)

import csv
import io

#def prediction_ki_batch(csv_file, model, model_esm, batch_converter):
def prediction_ki_batch(csv_file):
    if not csv_file:
        print("⚠ No file uploaded!")
        return []
    
    csv_file.seek(0)  # Reset file pointer
    decoded_file = csv_file.read().decode("utf-8").strip()
    
    if not decoded_file:
        print("⚠ Error: File is empty!")
        return []

    print("🔹 Raw File Content:\n", repr(decoded_file))  # Debugging
    
    # Detect delimiter (comma vs. tab)
    delimiter = ',' if decoded_file.count(',') > decoded_file.count('\t') else '\t'
    print(f"📝 Detected delimiter: {repr(delimiter)}")

    reader = csv.reader(decoded_file.splitlines(), delimiter=delimiter)
    
    # Read the header
    header = next(reader, None)
    print("📝 Header:", header)

    if not header or len(header) < 2:
        print("⚠ Error: Header not detected or missing columns!")
        return []

    results = []
    
    for row in reader:
        print("🔸 Row Data:", row)  # Debugging Step
        
        if len(row) >= 2:
            val1 = row[0].strip()  # Ligand
            val2 = row[1].strip()  # Protein
            
            # Call the prediction function
            
            result = prediction_ki(val1, val2)
            results.append({'ligand': val1, 'protein': val2, 'prediction': result})
    
    print("✅ Final Results:", results)
    return results


from django.http import JsonResponse
from rdkit import Chem
import re

def validate_inputs(request):
    input1 = request.GET.get("input1", "").strip()  # SMILES
    input2 = request.GET.get("input2", "").strip()  # Protein sequence
    input3 = request.FILES.get('input3_pl')

    def is_valid_smiles(smiles):
        """Check if a given SMILES string is valid using RDKit."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except:
            return False

    def is_valid_protein_sequence(sequence):
        """Check if the given protein sequence contains only valid amino acids."""
        valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # Standard 20 amino acids
        return bool(re.fullmatch(f"[{valid_amino_acids}]+", sequence))  # Match only valid characters

    # Validate SMILES (input1)
    if not is_valid_smiles(input1):
        return JsonResponse({"valid": False, "message": "Invalid SMILES string."})

    # Validate Protein Sequence (input2)
    if input2:
        if not is_valid_protein_sequence(input2):
            return JsonResponse({"valid": False, "message": "Invalid protein sequence. Only valid amino acids (ACDEFGHIKLMNPQRSTVWY) are allowed."})

    return JsonResponse({"valid": True})





# PP Workflow



def pp_results(sequence1, sequence2):
    
    time.sleep(2)
    pred= random.randint(0, 1)
    return pred

def prediction_pp_batch(csv_file):
    if not csv_file:
        print("⚠ No file uploaded!")
        return []
    
    csv_file.seek(0)  # Reset file pointer
    decoded_file = csv_file.read().decode("utf-8").strip()
    
    if not decoded_file:
        print("⚠ Error: File is empty!")
        return []

    print("🔹 Raw File Content:\n", repr(decoded_file))  # Debugging
    
    # Detect delimiter (comma vs. tab)
    delimiter = ',' if decoded_file.count(',') > decoded_file.count('\t') else '\t'
    print(f"📝 Detected delimiter: {repr(delimiter)}")

    reader = csv.reader(decoded_file.splitlines(), delimiter=delimiter)
    
    # Read the header
    header = next(reader, None)
    print("📝 Header:", header)

    if not header or len(header) < 2:
        print("⚠ Error: Header not detected or missing columns!")
        return []

    results = []
    
    for row in reader:
        print("🔸 Row Data:", row)  # Debugging Step
        
        if len(row) >= 2:
            val1 = row[0].strip()  # Protein_1
            val2 = row[1].strip()  # Protein_2
            
            # Call the prediction function
            result = str(pp_results(val1, val2))
            result = "Interacting" if result == "1.0" else "Non Interacting"
            results.append({'Protein_1': val1, 'Protein_2': val2, 'Interaction': result})
    
    print("✅ Final Results:", results)
    return results
