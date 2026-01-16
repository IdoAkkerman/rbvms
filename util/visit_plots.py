import sys
import os
import glob

OUTPUT_DIR = "plots"
RESOLUTION = (1024, 1024)

if not os.path.exists(OUTPUT_DIR):
    try:
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")
    except OSError as e:
        print(f"Error creating directory: {e}")
        sys.exit()

# Find all directories starting with "solution"
solution_dirs = glob.glob("solution*")
solution_dirs.sort()

if not solution_dirs:
    print("No 'solution*' directories found.")
    sys.exit()

# Setup global SaveWindow attributes
s = SaveWindowAttributes()
s.format = s.PNG
s.width, s.height = RESOLUTION
s.screenCapture = 0 
s.family = 0 
SetSaveWindowAttributes(s)

def save_variable(variable_name, output_path):
    """
    Clears old plots, adds the new variable, draws, and saves.
    """
    # clear previous plots ensures we don't overlay error on top of phi
    DeleteAllPlots() 

    # I disabled the mesh, because for high refinement, it's just pure black
    # AddPlot("Mesh", "main")
    
    # Add the Pseudocolor plot
    AddPlot("Pseudocolor", variable_name)
    
    DrawPlots()

    # Set filename
    s.fileName = output_path
    SetSaveWindowAttributes(s)
    SaveWindow()
    print(f"Saved: {output_path}.png")

# 4. Main Loop
for d in solution_dirs:
    # Find ALL .mfem_root files inside the solution folder
    root_files = glob.glob(os.path.join(d, "*.mfem_root"))
    
    # Sort them so we process step_00, step_01 in order
    root_files.sort()
    
    if not root_files:
        print(f"Skipping {d}: No .mfem_root files found.")
        continue
        
    print(f"--- Processing directory: {d} ---")

    # Inner loop: Iterate over every time step (every .mfem_root file)
    for db_path in root_files:
        
        # Extract step name (e.g. "step_000000")
        filename_only = os.path.basename(db_path)
        step_name = os.path.splitext(filename_only)[0]

        # Skip step 0 (boundary conditions only)
        if step_name.endswith("_000000"):
            continue

        # Open the specific time step
        if OpenDatabase(db_path, 0):
            
            if step_name.startswith("step_"):
                # Plot phi
                full_path_phi = os.path.join(OUTPUT_DIR, f"{d}_{step_name}_phi")
                try:
                    save_variable("phi", full_path_phi)
                except Exception as e:
                    print(f"Error plotting phi for {db_path}: {e}")

                # Plot error (skip for reference solution)
                if d != "solution":
                    full_path_err = os.path.join(OUTPUT_DIR, f"{d}_{step_name}_error")
                    try:
                        save_variable("error", full_path_err)
                    except Exception as e:
                        print(f"Error plotting error for {db_path}: {e}")
            
            elif step_name.startswith("tau_"):
                # Plot tau
                full_path_tau = os.path.join(OUTPUT_DIR, f"{d}_{step_name}_tau")
                try:
                    save_variable("tau", full_path_tau)
                except Exception as e:
                    print(f"Error plotting tau for {db_path}: {e}")

            # Delete plots *before* closing the database
            DeleteAllPlots()
            
            CloseDatabase(db_path)
        else:
            print(f"Failed to open database: {db_path}")

print("Batch processing complete.")
sys.exit()
