
import streamlit as st
st.set_page_config(layout="wide")  # Set the page to wide layout

import os
import pandas as pd
import streamlit as st
import numpy as np
from Bio import SeqIO
import io
import tempfile
import glob
import re
import subprocess
import streamlit as st
import zipfile


# Ensure the folder path is valid and correct for your system
def normalize_folder_path(path):
    """Normalize path to work on different OS platforms."""
    if os.name == 'nt':
        return os.path.normpath(path)  # For Windows
    else:
        return os.path.abspath(path)  # For Unix/Linux/Mac

def validate_input(output_folder, uploaded_pdb, color1, color2, color3):
    """Validates that necessary inputs are provided."""
    if not os.path.exists(output_folder):
        st.error(f"The directory does not exist: {output_folder}")
        return False
    
    if uploaded_pdb is None:
        st.error("Please upload a PDB file.")
        return False
    
    # Ensure colors are provided
    if not any([color1, color2, color3]):
        st.error("Please select at least one color.")
        return False
    
    return True

def hex_to_rgb(hex_color):
    """Converts hexadecimal color code to RGB tuple."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def setup_pymol_color_palette(color1, color2, color3):
    """Sets up the color gradient in PyMOL based on the given colors."""
    
    
    # Define the PyMOL color palette from the provided colors
    colors = [color1, color2, color3]
    
    
    pymol_palette = []
    pymol_commands = []
    
    for i, color in enumerate(colors):
        if color:
            pymol_palette.append(f"color{i+1}")  # Add color name to the palette
            rgb = hex_to_rgb(color)
            pymol_commands.append(f"set_color color{i+1}, [{rgb[0]/255.0}, {rgb[1]/255.0}, {rgb[2]/255.0}]")
            
        else:
            st.write(f"Skipped invalid color: {color}")
    
    # Join the palette with underscores for PyMOL formatting
    palette_str = "_".join(pymol_palette)
    
    
    return pymol_commands, palette_str


# Function to read the FASTA file from an uploaded file object
def read_fasta(uploaded_fasta):
    """Reads the content of the uploaded FASTA file and returns DataFrames for sequence, delta, and weights."""
    try:
        # Decode the content of the uploaded file and parse it as a FASTA sequence
        fasta_content = uploaded_fasta.read().decode("utf-8")
        fasta_sequences = list(SeqIO.parse(io.StringIO(fasta_content), "fasta"))
        
        # Check if any sequences were found in the FASTA file
        if len(fasta_sequences) == 0:
            st.error("The FASTA file does not contain any valid sequences.")
            return None, None, None
        
        # Extract the sequence from the first entry in the FASTA file
        sequence = str(fasta_sequences[0].seq)
        amino_acids = list(range(1, len(sequence) + 1))
        
        # Create DataFrames for the sequence, delta, and weights
        fasta_df = pd.DataFrame({'sequence': list(sequence)})
        data = {'seq': list(sequence), 'AminoAcid': amino_acids}
        delta = pd.DataFrame(data)
        weights = pd.DataFrame(data)
        
        # Adjust indices to start from 1 instead of 0 for easier processing
        delta.index += 1
        weights.index += 1

        return fasta_df, delta, weights
    
    except Exception as e:
        # Handle any errors that occur during FASTA file reading and parsing
        st.error(f"Error reading FASTA file: {str(e)}")
        return None, None, None

# Function to process uptake data for each state

def process_data_per_state(df_state_data, state, output_folder):
    """Processes the uptake data for a given state and saves it to CSV, adjusting exposure times as needed."""
    # Filter the data for the specified state
    df_state = df_state_data[df_state_data['State'] == state]
    
    
    
    # Calculate the sum of uptake per peptide
    sum_per_peptide = df_state.groupby('Sequence')[['Uptake']].sum()
    final_rows = []
    
    # Process each peptide group
    for peptide, group in df_state.groupby('Sequence', sort=False):
        final_rows.append(group)
        
        # Retrieve summed values for uptake and define a row with these summary values
        sums = sum_per_peptide.loc[peptide]
        sum_row_values = [peptide, group['Start'].iloc[0], group['End'].iloc[0], state, '100000000', sums['Uptake']]
        sum_row = pd.DataFrame([sum_row_values], columns=df_state.columns)
        final_rows.append(sum_row)

    # Concatenate final rows into a DataFrame for output
    Uptake_per_pep = pd.concat(final_rows, ignore_index=True)
    state_folder = os.path.join(output_folder, state)
    if not os.path.exists(state_folder):
        os.makedirs(state_folder)
    
    # Export to CSV
    output_filename = f'Uptake_data_per_pep_{state}.csv'
    output_path = os.path.join(state_folder, output_filename)
    #Uptake_per_pep.to_csv(output_path, header=True, index=False)
    
    return output_path, state_folder, Uptake_per_pep


# Adjusted process_data function to return a list of CSV file paths
def process_data(diff_data, state_folder, df, delta, weights,state):
    """Processes uptake data and returns max, min, and list of CSV file paths."""
    unique_time = sorted(diff_data['Exposure'].astype(float).unique())
    max_flattened_data = min_flattened_data = None
    csv_file_paths = []  # To store paths of generated CSV files

    for e in unique_time:
        if e == 0:
            continue
        
        # Processing logic remains the same as before
        diff_d_per_time = diff_data[diff_data['Exposure'] == e][['Sequence', 'Start', 'End', 'State', 'Uptake']].reset_index(drop=True)
        for j, row in diff_d_per_time.iterrows():
            n1, n2 = int(row['Start']), int(row['End'])
            length = list(range(n1 + 1, n2 + 1))
            sequence = row['Sequence']
            P_count = sequence.count('P')
            n_P = len(sequence) - P_count - 1
            a = row['Uptake'] / n_P
            delta.loc[length, str(j)] = a
            weights.loc[length, str(j)] = n_P

        delta_time_v2 = delta.iloc[:, 2:].apply(pd.to_numeric, errors='coerce')
        weights_v2 = weights.iloc[:, 2:].apply(pd.to_numeric, errors='coerce')
        weights_sum = weights_v2.sum(axis=1).replace(0, np.nan)
        weighted_average = delta_time_v2 * weights_v2
        delta_time_v2['Flattened_data'] = weighted_average.sum(axis=1) / weights_sum
        
        if e == unique_time[-1]:
            max_flattened_data = delta_time_v2['Flattened_data'].max()
            min_flattened_data = delta_time_v2['Flattened_data'].min()

        data = pd.DataFrame({
            'Flattened_data': delta_time_v2['Flattened_data'],
            'AminoAcid': delta['AminoAcid'],
            'seq': delta['seq']
        })
        data.loc[data['seq'] == 'P', 'Flattened_data'] = np.nan
        
        # Save each CSV file
        output_filename = f'Uptake_avr_{state}_{e}.csv'
        csv_path = os.path.join(state_folder, output_filename)
        data.to_csv(csv_path, index=True, header=True)
        csv_file_paths.append(csv_path)  # Add path to list


    return max_flattened_data, min_flattened_data, csv_file_paths




def export_to_pymol_aa_diff(folder_path, min_val, max_val, pdb_file, chain, color1, color2, color3, state):
    """Exports amino acid uptake data to PyMOL and applies a custom color gradient palette for each state."""
    
    # Check if the folder exists
    if not os.path.exists(folder_path):
        st.error(f"Directory does not exist: {folder_path}")
        return
    
    # Find CSV files matching the specific state
    file_list = glob.glob(os.path.join(folder_path, f'Uptake_avr_{state}_*.csv'))
    all_files = glob.glob(os.path.join(folder_path, '*.csv'))
    

    if not file_list:
        st.error(f"No CSV files found for state '{state}' in directory.")
        return
    
    # Sort the CSV files by numeric values
    file_list.sort(key=lambda x: float(re.search(r"(\d+\.?\d*)", os.path.split(x)[-1]).group(1)))
    
    # Temporary directory for PyMOL processing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Save the uploaded PDB file
        uploaded_pdb_path = os.path.join(temp_dir, "uploaded_structure.pdb")
        with open(uploaded_pdb_path, 'wb') as f:
            f.write(pdb_file.getvalue())
        
        # Generate PyMOL commands
        pymol_commands, color_palette_str = setup_pymol_color_palette(color1, color2, color3)
        
        # Create the PML file for this state
        pml_file_path = os.path.join(temp_dir, f"pymol_script_{state}.pml")
        with open(pml_file_path, 'w') as pml_file:
            # Write PyMOL setup commands
            pml_file.write("set grid_mode, on\n")
            for command in pymol_commands:
                pml_file.write(f"{command}\n")
            
            # Load PDB structure
            pml_file.write(f"load {uploaded_pdb_path}, structure\n")
            
            for csv_file_path in file_list:
                csv_filename = os.path.basename(csv_file_path)
                csv_filename = csv_filename.replace(" ", "_").replace("-", "_")  # Add more replacements as needed

                pml_file.write(f"create {csv_filename}, structure\n")
                PyClusters = pd.read_csv(csv_file_path, header=0)
                PyClusters = PyClusters[['AminoAcid', 'Flattened_data']].dropna()
                PyClusters['Flattened_data'].replace(0, np.nan, inplace=True)
                
                pml_file.write(f"alter {csv_filename}, b=0\n")
                pml_file.write(f"color white, {csv_filename}\n")
                
                for _, row in PyClusters.iterrows():
                    residue_number = int(row['AminoAcid'])
                    selection = f"{csv_filename} and resi {residue_number} and chain {chain}"
                    b_value = row['Flattened_data']
                    pml_file.write(f"alter {selection}, b={b_value}\n")
                
                pml_file.write(f"spectrum b, {color_palette_str}, selection={csv_filename}, minimum={min_val}, maximum={max_val}, byres=1\n")
           
            # Add grey coloring for residues with missing or zero uptake data
            for csv_file_path in file_list:
                csv_filename = os.path.basename(csv_file_path)
                csv_filename = csv_filename.replace(" ", "_").replace("-", "_")  # Add more replacements as needed

                PyClusters = pd.read_csv(csv_file_path, header=0)
                empty_data = PyClusters[PyClusters['Flattened_data'].isna() | 
                        (PyClusters['Flattened_data'] == 0) ]
                
                for _, row in empty_data.iterrows():
                    residue_number = int(row['AminoAcid'])
                    selection = f"{csv_filename} and resi {residue_number} and chain {chain}"
                    pml_file.write(f"color grey, {selection}\n")
                    
                    
            pml_file.write(f"delete structure\n")
            
            
            # Additional PyMOL visualization settings
            pml_file.write("set transparency_mode, 3\n")
            pml_file.write("set ray_trace_mode, 0\n")
            pml_file.write("set volume_layers, 1000\n")
            pml_file.write("set solvent_radius, 2.5\n")
            pml_file.write("set depth_cue, 1\n")
            pml_file.write("set fog_start, 0.4\n")
            pml_file.write("set cartoon_fancy_helices, 1\n")
            pml_file.write("set cartoon_transparency, 0\n")
            pml_file.write("set ambient, 0.35\n")
            pml_file.write("set direct, 0.7\n")
            pml_file.write("set reflect, 0.3\n")
            pml_file.write("set spec_power, 1250\n")
            pml_file.write("set spec_reflect, 50\n")
            pml_file.write("set ambient_occlusion_scale, 27\n")
            pml_file.write("set ambient_occlusion_smooth, 5\n")
            
            
            
            # Define PyMOL session file for this state
            pse_session_file = os.path.join(folder_path, f"UD_per_aa_{state}.pse")
            pml_file.write(f"save {pse_session_file}\n")
            st.success(f"PyMOL session for {state} is ready for download!")
            st.write(f"Download the session from: {pse_session_file}")
            
            
        # Run PyMOL script for the current state
        subprocess.run(["pymol", "-cq", pml_file_path])
        
        return pse_session_file




def load_pymol_session(state):
    if f"pse_session_{state}" in st.session_state:
        pse_session_file = st.session_state[f"pse_session_{state}"]
        # Provide the download link for the saved session file
        st.download_button(label="Download PyMOL Session", data=pse_session_file, file_name=f"UD_per_aa_{state}.pse", mime="application/octet-stream")
    else:
        st.warning(f"No session available for {state}. Please generate it first.")






# Custom CSS for styling
st.markdown("""
    <style>
        /* Smaller font size for the entire app */
        body {
            font-size: 12px;
        }

        /* Smaller font size for the title */
        h1, h2, h3, h4, h5, h6 {
            font-size: 20px;
        }

        /* Smaller font size for subheaders */
        .streamlit-expanderHeader {
            font-size: 14px;
        }

        /* Adjusting button size */
        .css-1emrehy.edgvbvh3 {
            font-size: 12px;  /* Decrease font size */
            padding: 5px 10px;  /* Reduce button padding */
        }

        /* Styling all columns to have white background and white outline */
        div[data-testid="stHorizontalBlock"] > div {
            background-color: white;  /* Set all columns to white */
            border: 1px solid #cccccc;  /* Change the border color to white */
            padding: 10px;  /* Add some padding inside the columns */
        }
    </style>
""", unsafe_allow_html=True)

# ------------------------------- User Input Section -------------------------------
st.title("Amino Acid Uptake Data Processing and Visualization")
st.markdown("Process amino acid uptake data and export visualizations to PyMOL.")

# ------------------------------- Process Button Section -------------------------------

# Create the layout with 3 main columns
col1, col2, col3 = st.columns([1, 1, 1])  # Equal-width columns

# Column 1: Folder path input
with col1:
    st.markdown('<div class="custom-section">', unsafe_allow_html=True)
    st.markdown("### Output Folder Path")
    st.markdown('</div>', unsafe_allow_html=True)

    output_folder = st.text_input(
        "Enter the output folder path:", 
        "C:\\Users\\mk678\\OneDrive - University of Exeter\\Python_project_Monika\\NEw_scripts\\output_tests", 
        key="unique_output_folder_path"
    )


    # File uploads for CSV, FASTA, and PDB files
    st.subheader("File Uploads")
    uploaded_csv = st.file_uploader("Upload CSV File:", type=["csv"])
    uploaded_fasta = st.file_uploader("Upload FASTA File:", type=["fasta"])
    uploaded_pdb = st.file_uploader("Upload PDB File:", type=["pdb"])

# Column 2: Input fields
with col2:
    st.subheader("Protein Chain Selection")
    chain_input = st.text_input("Enter the chain (e.g., A, B, C):", "A")

    st.subheader("Color Gradient Settings")
    st.write("Choose colors for the amino acid uptake gradient:")
    color1 = st.color_picker("Color for Low Uptake (e.g., Blue):", "#0000FF")
    color2 = st.color_picker("Color for Medium Uptake (e.g., White):", "#FFFFFF")
    color3 = st.color_picker("Color for High Uptake (e.g., Red):", "#FF0000")

# Column 3: Button to process data
with col3:
    if st.button("Process Data and Export to PyMOL"):
        if validate_input(output_folder, uploaded_pdb, color1, color2, color3):
            if uploaded_csv and uploaded_fasta and uploaded_pdb:
                df_state_data = pd.read_csv(uploaded_csv)
                required_columns = ['Sequence', 'Start', 'End', 'State', 'Exposure', 'Uptake']
                df_state_data = df_state_data[required_columns]
                df, delta, weights = read_fasta(uploaded_fasta)

                if df is not None and delta is not None and weights is not None:
                    # Create a ZIP file for PyMOL sessions
                    pse_zip_buffer = io.BytesIO()
                    csv_zip_buffer = io.BytesIO()
                    with zipfile.ZipFile(pse_zip_buffer, mode='w', compression=zipfile.ZIP_DEFLATED) as pse_zip_file, \
                         zipfile.ZipFile(csv_zip_buffer, mode='w', compression=zipfile.ZIP_DEFLATED) as csv_zip_file:
                        
                        for state in df_state_data['State'].unique():
                            
                            output_path, state_folder, Uptake_per_pep = process_data_per_state(
                                df_state_data, state, output_folder
                            )
                            max_val, min_val, csv_file_paths = process_data(Uptake_per_pep, state_folder, df, delta, weights, state)

                            st.markdown(f"### Results for State: **{state}**")
                            st.write(f"**Max Flattened Data:** {round(max_val, 3)}")
                            st.write(f"**Min Flattened Data:** {round(min_val, 3)}")

                            pse_session_file = export_to_pymol_aa_diff(
                                    state_folder, min_val, max_val, uploaded_pdb, chain_input, color1, color2, color3, state
                                        )
   
                            # Add PyMOL session file to ZIP
                            with open(pse_session_file, "rb") as session_file:
                                pse_zip_file.writestr(f"UD_per_aa_{state}.pse", session_file.read())

                            # Add CSV files to ZIP
                            for csv_file in csv_file_paths:
                                csv_zip_file.write(csv_file, os.path.basename(csv_file))

                    # Seek back to the start of the buffers
                    pse_zip_buffer.seek(0)
                    csv_zip_buffer.seek(0)

                    # Save the ZIP files in session state
                    st.session_state["pse_zip"] = pse_zip_buffer.getvalue()
                    st.session_state["csv_zip"] = csv_zip_buffer.getvalue()
                    st.session_state["zip_ready"] = True
                else:
                    st.error("Error processing FASTA file. Please verify the file format and content.")
            else:
                st.error("All files (CSV, FASTA, PDB) must be uploaded to proceed.")
        else:
            st.error("Invalid input. Please check the output folder path and color gradient settings.")

# Display download buttons if ZIP files are ready
if st.session_state.get("zip_ready"):
    st.download_button(
        label="Download All PyMOL Sessions (Zipped)",
        data=st.session_state["pse_zip"],
        file_name="PyMOL_Sessions.zip",
        mime="application/zip"
    )
    st.download_button(
        label="Download All Processed CSVs (Zipped)",
        data=st.session_state["csv_zip"],
        file_name="Processed_CSVs.zip",
        mime="application/zip"
    )
