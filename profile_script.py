import os
import sys
import pstats
import pandas as pd
import subprocess
from datetime import datetime
import cProfile

def profile_script():
    """Profiles the main.py script and saves the results in a CSV file with git commit and a timestamp."""
    # Profile the main.py script
    cProfile.run('exec(open("main.py").read())', 'profile_results.prof')

    # Get current timestamp and Git commit hash
    timestamp = datetime.now().strftime('%d-%m-%Y_%H-%M-%S')
    try:
        commit_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip().decode('utf-8')
    except subprocess.CalledProcessError:
        commit_hash = "no_git_repo"

    # Create a directory to save the profiling results (if it doesn't exist)
    profiling_directory = 'profiling_directory'
    if not os.path.exists(profiling_directory):
        os.makedirs(profiling_directory)

    # Define a filename based on timestamp and Git commit hash
    filename_prefix = f"profiling_{timestamp}_{commit_hash}"

    # Load profiling results from the file
    p = pstats.Stats('profile_results.prof')

    # Sort by cumulative time
    p.sort_stats('cumulative')

    # Create an empty list to hold the profiling data
    profile_data = []

    # Loop through the profile statistics and extract relevant statistics
    for func in p.stats.items():
        func_name = func[0]
        calls = func[1][0]  # Number of calls
        total_time = func[1][2]  # Total time spent in the function
        per_call_time = total_time / calls if calls > 0 else 0
        cum_time = func[1][3]  # Cumulative time spent in the function

        # Append each function's data as a dictionary
        profile_data.append({
            'Function': func_name,
            'Calls': calls,
            'Total Time (s)': total_time,
            'Per Call Time (s)': per_call_time,
            'Cumulative Time (s)': cum_time
        })

    # Create a Pandas DataFrame from the profiling data
    df = pd.DataFrame(profile_data)

    # Set the index names explicitly for a MultiIndex
    df.index.names = ['Index']
    df = df.sort_values(by='Cumulative Time (s)', ascending=False)

    # Save the DataFrame as a CSV
    csv_filename = os.path.join(profiling_directory, f"{filename_prefix}.csv")
    df.to_csv(csv_filename, index=False)

    # Print completion message
    print(f"Profiling data and plot saved in '{profiling_directory}' directory with filename '{filename_prefix}'.")

def main():
    """Main entry point for the script."""
    if len(sys.argv) != 2 or sys.argv[1] != 'run':
        print("Usage: python main.py run")
        sys.exit(1)

    print("Running profiling...")
    profile_script()

if __name__ == "__main__":
    main()