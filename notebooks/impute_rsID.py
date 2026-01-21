import pandas as pd
import sqlite3
import os
import glob
from tqdm import tqdm

if __name__ == "__main__":
    # ----------------------------
    # Configuration and Setup
    # ----------------------------

    # Path to your database
    db_path = './variants.db'

    # Directory containing your CSV files
    csv_directory = './CancerPGS/Annotated' 

    # Path to save the updated CSV files
    output_directory = './CancerPGS/Annotated/Clean_ReAnnoted'

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # SQLite maximum variables per query (default is 999)
    SQLITE_MAX_VARIABLE_NUMBER = 999

    # Maximum number of positions per chunk (since each position uses two variables)
    max_positions_per_chunk = SQLITE_MAX_VARIABLE_NUMBER // 2

    # ----------------------------
    # Function Definitions
    # ----------------------------

    # Function to check if hm_rsID is valid
    def is_valid_rsID(rsID):
        return isinstance(rsID, str) and rsID.startswith('rs')

    # Function to query a chunk of positions
    def query_positions(chunk, cursor):
        placeholders = ', '.join(['(?, ?)'] * len(chunk))
        query = f'SELECT chrom, pos, id FROM variants WHERE (chrom, pos) IN ({placeholders})'
        # Flatten the list of tuples for the query parameters
        query_params = [item for position in chunk for item in position]
        cursor.execute(query, query_params)
        results = cursor.fetchall()
        # Return a dictionary mapping (chrom, pos) to rsID
        return { (str(row[0]), int(row[1])): row[2] for row in results }

    # ----------------------------
    # Connect to the Database
    # ----------------------------

    # Check if the database exists
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database not found at path: {db_path}")

    # Connect to the database
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Ensure that the variants table has an index on (chrom, pos)
    c.execute('CREATE INDEX IF NOT EXISTS idx_chrom_pos ON variants (chrom, pos)')

    # ----------------------------
    # Process Each CSV File
    # ----------------------------

    # Get a list of all CSV files in the directory
    csv_files = glob.glob(os.path.join(csv_directory, '*.csv'))

    print(f"Found {len(csv_files)} CSV files in directory: {csv_directory}")

    for csv_file_path in csv_files:
        print(f"\nProcessing file: {os.path.basename(csv_file_path)}")

        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file_path)

        # Ensure required columns exist
        required_columns = {'hm_rsID', 'hm_chr', 'hm_pos'}
        if not required_columns.issubset(df.columns):
            print(f"Error: Required columns missing in {os.path.basename(csv_file_path)}")
            continue  # Skip this file

        # Ensure hm_chr and hm_pos are the correct data types
        df['hm_chr'] = df['hm_chr'].astype(str)
        df['hm_pos'] = df['hm_pos'].astype(int)

        # Identify invalid hm_rsIDs
        invalid_rsID_mask = ~df['hm_rsID'].apply(is_valid_rsID)
        invalid_entries = df[invalid_rsID_mask]

        num_invalid_before = invalid_rsID_mask.sum()
        print(f"Number of invalid hm_rsID before imputation: {num_invalid_before}")

        if num_invalid_before == 0:
            print("No invalid hm_rsID found. Skipping imputation.")
            # Save a copy of the original file to the output directory
            output_file_path = os.path.join(output_directory, os.path.basename(csv_file_path))
            df.to_csv(output_file_path, index=False)
            continue  # Move to the next file

        # Get unique (hm_chr, hm_pos) pairs from invalid entries
        unique_positions = invalid_entries[['hm_chr', 'hm_pos']].drop_duplicates()

        # Initialize an empty dictionary to store the mappings
        position_to_rsID = {}

        # Prepare positions for querying
        positions = list(unique_positions.itertuples(index=False, name=None))

        # Split the positions into chunks and query
        for i in tqdm(range(0, len(positions), max_positions_per_chunk), desc="Querying database"):
            chunk = positions[i:i+max_positions_per_chunk]
            mapping = query_positions(chunk, c)
            position_to_rsID.update(mapping)

        # Function to impute rsID using the mapping
        def impute_rsID(row):
            if not is_valid_rsID(row['hm_rsID']):
                key = (str(row['hm_chr']), int(row['hm_pos']))
                return position_to_rsID.get(key, row['hm_rsID'])  # If not found, keep original hm_rsID
            else:
                return row['hm_rsID']

        # Apply the imputation
        df['hm_rsID'] = df.apply(impute_rsID, axis=1)

        # Count the number of invalid hm_rsID after imputation
        num_invalid_after = (~df['hm_rsID'].apply(is_valid_rsID)).sum()
        print(f"Number of invalid hm_rsID after imputation: {num_invalid_after}")

        # Save the updated DataFrame to the output directory
        output_file_name = os.path.basename(csv_file_path).replace('.csv', '_updated.csv')
        output_file_path = os.path.join(output_directory, output_file_name)
        df.to_csv(output_file_path, index=False)

        print(f"Updated file saved to: {output_file_path}")

    # ----------------------------
    # Cleanup
    # ----------------------------

    # Close the database connection
    conn.close()

    print("\nProcessing of all CSV files completed.")
