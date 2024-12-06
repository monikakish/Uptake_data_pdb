import streamlit as st
import subprocess
import os

# Streamlit app title
st.title("PyMOL Integration Test")

# Upload a PDB file
uploaded_file = st.file_uploader("Upload a PDB file", type=["pdb"])

if uploaded_file:
    # Save uploaded file temporarily
    temp_dir = "temp"
    os.makedirs(temp_dir, exist_ok=True)
    pdb_path = os.path.join(temp_dir, uploaded_file.name)
    with open(pdb_path, "wb") as f:
        f.write(uploaded_file.read())

    # Run PyMOL using subprocess
    st.write("## PyMOL Visualization")
    output_image = os.path.join(temp_dir, "pymol_output.png")
    pymol_script = os.path.join(temp_dir, "script.pml")

    try:
        # Create a PyMOL script to automate visualization
        with open(pymol_script, "w") as script:
            script.write(f"""
load {pdb_path}
bg_color white
show cartoon
color cyan
png {output_image}, width=800, height=600, dpi=300
quit
            """)

        # Run PyMOL
        subprocess.run(["pymol", "-c", pymol_script], check=True)

        # Display the generated image
        if os.path.exists(output_image):
            st.image(output_image, caption="PyMOL Visualization", use_column_width=True)
        else:
            st.error("PyMOL did not generate an image.")
    except Exception as e:
        st.error(f"An error occurred while running PyMOL: {e}")
    finally:
        # Cleanup files
        if st.button("Clean Up"):
            for file in [pdb_path, output_image, pymol_script]:
                if os.path.exists(file):
                    os.remove(file)
            os.rmdir(temp_dir)
            st.success("Temporary files cleaned up.")
