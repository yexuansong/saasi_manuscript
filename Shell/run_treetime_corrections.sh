#!/bin/bash

# Loop through all 1000 simulations
for i in {1..1000}; do
    echo "Processing simulation $i..."
    
    SIM_DIR="sim_${i}"
    
    # TreeTime with multiple sampling bias corrections
    cd "$SIM_DIR"
    
    # Define correction factors
    CORRECTION_FACTORS=(1.0 2.0 5.0 10.0)
    
    # Loop through each correction factor
    for correction in "${CORRECTION_FACTORS[@]}"; do
        # Create output directory name (replace . with _ for folder name)
        correction_suffix=$(echo $correction | tr '.' '_')
        output_dir="treetime_results_${correction_suffix}"
        
        mkdir -p "$output_dir"
        
        # Run TreeTime with sampling bias correction
        treetime_start=$(date +%s)
        treetime mugration \
            --tree "tree_${i}.nwk" \
            --states "metadata_${i}.csv" \
            --attribute numeric_state \
            --confidence \
            --sampling-bias-correction $correction \
            --outdir "$output_dir" \
            --verbose 0 > /dev/null 2>&1
        treetime_success=$?
        treetime_time=$(($(date +%s) - treetime_start))
        
        # Determine success status
        if [ $treetime_success -eq 0 ]; then
            success_status="TRUE"
        else
            success_status="FALSE"
        fi
        
        # Append to external_timing.csv
        echo "${i},treetime_${correction},${treetime_time},${success_status}" >> "external_timing.csv"
    done
    
    # Go back to base directory
    cd ..
    
    echo "Completed simulation $i"
done

echo "All simulations completed!"
