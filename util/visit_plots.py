import sys
import os
import glob

OUTPUT_DIR = "plots"
RESOLUTION = (1024, 1024)
# We will set the Multires resolution to 5 to ensure smoothness
SMOOTHING_LEVEL = 5

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

def save_variable(variable_name, output_path, ref_level):
    """
    Clears old plots, adds the new variable, applies Multires smoothing,
    conditionally adds mesh, draws, and saves.
    """
    DeleteAllPlots()

    # 1. Add the Pseudocolor plot
    AddPlot("Pseudocolor", variable_name)

    # 2. Apply Multires Control for smoothing
    AddOperator("MultiresControl", 0)
    m_atts = MultiresControlAttributes()
    m_atts.resolution = SMOOTHING_LEVEL
    SetOperatorOptions(m_atts)

    # 3. Conditionally Add Mesh Plot (if refinement level <= 3)
    if ref_level <= 3:
        AddPlot("Mesh", "main")
        # Optional: Set mesh to be slightly transparent or a specific color
        # so it doesn't obscure the Pseudocolor data
        m_plot_atts = MeshAttributes()
        m_plot_atts.legendFlag = 0 # Hide mesh legend to keep it clean
        SetPlotOptions(m_plot_atts)

    DrawPlots()

    # Set filename and save
    s.fileName = output_path
    SetSaveWindowAttributes(s)
    SaveWindow()
    print(f"Saved: {output_path}.png (Ref: {ref_level}, Smooth: {SMOOTHING_LEVEL})")

# 4. Main Loop
for d in solution_dirs:
    # Logic to extract refinement level from folder name (e.g., solution_r4 -> 4)
    # If it's just 'solution', we treat it as a high-ref reference
    try:
        if "_r" in d:
            current_ref_level = int(d.split('_r')[-1])
        else:
            current_ref_level = 99 # Reference solution
    except ValueError:
        current_ref_level = 99

    root_files = glob.glob(os.path.join(d, "*.mfem_root"))
    root_files.sort()

    if not root_files:
        continue

    print(f"--- Processing directory: {d} (Level {current_ref_level}) ---")

    for db_path in root_files:
        filename_only = os.path.basename(db_path)
        step_name = os.path.splitext(filename_only)[0]

        if step_name.endswith("_000000"):
            continue

        if OpenDatabase(db_path, 0):
            if step_name.startswith("step_"):
                # Plot phi
                full_path_phi = os.path.join(OUTPUT_DIR, f"{d}_{step_name}_phi")
                try:
                    save_variable("phi", full_path_phi, current_ref_level)
                except Exception as e:
                    print(f"Error plotting phi for {db_path}: {e}")

                # Plot error
                if d != "solution":
                    full_path_err = os.path.join(OUTPUT_DIR, f"{d}_{step_name}_error")
                    try:
                        save_variable("error", full_path_err, current_ref_level)
                    except Exception as e:
                        print(f"Error plotting error for {db_path}: {e}")

            elif step_name.startswith("tau_"):
                full_path_tau = os.path.join(OUTPUT_DIR, f"{d}_{step_name}_tau")
                try:
                    save_variable("tau", full_path_tau, current_ref_level)
                except Exception as e:
                    print(f"Error plotting tau for {db_path}: {e}")

            DeleteAllPlots()
            CloseDatabase(db_path)

print("Batch processing complete.")
sys.exit()